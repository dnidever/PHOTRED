#!/usr/bin/env python

import os
import time
import numpy as np
import time
import shutil
from datetime import datetime
import subprocess
import tempfile
import logging
import re
from glob import glob
import warnings
import traceback
from astropy.io import fits,ascii
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from dlnpyutils import utils as dln
from . import utils,io,iraf,combine as comb,mkopt,getpsf

# Filtering warnings
warnings.simplefilter('ignore', category=AstropyWarning)
#warnings.filterwarnings(action="ignore", message=r'FITSFixedWarning:*')

def allfprep(filename,xoff=0.0,yoff=0.0,maxiter=1,scriptsdir=None,
             detectprog='sextractor',logfile=None,maskfile=None):
    """
    IDL version of Jamie's allfprep.cl that runs SExtractor 
    iteratively on stacked images and then runs ALLSTAR 
    to get PSF sources. 
 
    The final ALS file will be called file+'_allf.als' 
    A file with all of the sources that SExtractor found 
    that can be matched to the ALS file using IDs is called 
    file+'_allf.sex' 
 
    Parameters
    ----------
    filename : str
       Filename of the stacked images 
    xoff : float, optional
       The offset in X between the original and shifted images 
         xorig = xshift + xoff 
    yoff : float, optional
       The offset in Y between the original and shifted images 
         yorig = yshift + yoff 
    maxiter : int, optional
       Maximum number of times to iterate. By default, maxiter=1
    scriptsdir : str, optional
       Directory that contains the scripts. 
    detectprog : str, optional
       The program to use for source detection.  The options 
         are 'SEXTRACTOR' (the default) or 'DAOPHOT'. 
    logfile : str, optional
       A logfile to print output to. 
    maskfile : str, optional
       The name of a mask/weight file to use for the combined image. 

    Returns
    -------
    als : table,
       Final ALLSTAR table.
 
    Example
    -------
    als = allfprep('ccd1001_comb.fits',scriptsdir=scriptsdir,logfile=logfile)
 
    By D. Nidever    February 2008 (copied from Jamie's allfprep.cl) 
    Translated to Python by D. Nidever,  April 2022
    """
 
    # Logfile
    if logfile is None:
        logger = dln.basiclogger()
     
    # Make sure file exists 
    if os.path.exists(filename)==False:
        raise ValueError(filename+' NOT FOUND')
     
    # What program are we using for source detection? 
    #------------------------------------------------ 
    detectprog = detectprog.lower()
    if detectprog == 'sex': 
        detectprog = 'sextractor' 
    if detectprog=='dao': 
        detectprog = 'daophot' 
    if detectprog != 'sextractor' and detectprog != 'daophot': 
        raise ValueError('DETECTPROG = '+detectprog+' NOT AN OPTION.  Use sextractor or daophot')
     
    # Copy the scripts 
    #--------------------- 
    # No scriptsdir
    if scriptsdir is None:
        raise ValueError('scriptsdir NOT INPUT')
    # Check if the scripts exist in the current directory 
    scripts = ['default.sex','default.param','default.nnw','default.conv'] 
    nscripts = len(scripts) 
    # Loop through the scripts 
    for i in range(nscripts):
        exists = os.path.exists(scriptsdir+'/'+scripts[i])
        if exists:
            size = os.stat(scripts[i])
        else:
            size = 0
        curexists = os.path.exists(scripts[i])
        if curexists:
            cursize = os.stat(scripts[i])
        else:
            cursize = 0
        # No file
        if exists==False or size==0:
            raise ValueError(scriptsdir+'/'+scripts[i]+' NOT FOUND or EMPTY')
        # Check if the two files are the same size, if not copy it
        if size!=cursize:
            if curexists:
                os.remove(scripts[i])
            shutil.copyfile(scriptsdir+'/'+scripts[i],scripts[i])

    # Check that the SEXTRACTOR program exists
    out = subprocess.check_output(['which','sex'],shell=False)
    if type(out) is bytes: out = out.decode()
    out = out.strip()
    if os.path.exists(out)==False:
        raise ValueError('No sextractor program found')
    # Check that the ALLSTAR program exists
    out = subprocess.check_output(['which','allstar'],shell=False)
    if type(out) is bytes: out = out.decode()
    out = out.strip()
    if os.path.exists(out)==False:
        raise ValueError('No allstar program found')    
    # Check that the DAOPHOT program exists
    out = subprocess.check_output(['which','daophot'],shell=False)
    if type(out) is bytes: out = out.decode()
    out = out.strip()
    if os.path.exists(out)==False:
        raise ValueError('No daophot program found') 
     
    base,ext = os.path.splitext(os.path.basename(filename))

    # Read the .opt file
    opt = io.readopt(base+'.opt')
    satlevel = opt['HI']
    gain = opt['GA']
    fwhm = opt['FW']
     
    # Get pixel scale
    scale = io.getpixscale(filename)
    if scale is  None: # default 
        scale = 0.5 
     
    # Filenames 
    origfile = filename
    subbase = base+'_sub' 
    subfile = subbase+'.fits' 
    catfile = subbase+'.cat' 
    sexfile = base+'_allf.sex'
    if os.path.exists(subfile): os.remove(subfile)
    shutil.copyfile(filename,subfile)
     
    # Make customized SEXTRACTOR file 
    #-------------------------------- 
    if (detectprog == 'sextractor'): 
        sexconfigfile = base+'.sex'
        if os.path.exists(sexconfigfile): os.remove(sexconfigfile)
        shutil.copyfile('default.sex',sexconfigfile)
        sexlines = dln.readlines(sexconfigfile)
        sexlines2 = np.char.array(sexlines).copy()
        for i,line in enumerate(sexlines):
            # CATALOG_NAME
            if line.find('CATALOG_NAME')>-1:
                sexlines2[i] = 'CATALOG_NAME    '+catfile+'  # name of the output catalog' 
            # SATUR_LEVEL
            elif line.find('SATUR_LEVEL')>-1:
                sexlines2[i] = 'SATUR_LEVEL     '+str(satlevel)+'  # level (in ADUs) at which arises saturation' 
            # GAIN
            elif line.find('GAIN')>-1:
                sexlines2[i] = 'GAIN            '+str(gain)+'  # detector gain in e-/ADU.' 
            # PIXEL_SCALE
            elif line.find('PIXEL_SCALE')>-1:
                sexlines2[i] = 'PIXEL_SCALE     '+str(scale)+' # size of pixel in arcsec (0=use FITS WCS info).' 
            # SEEING_FWHM
            elif line.find('SEEING_FWHM')>-1:
                fwhmas = float(fwhm)*float(scale) 
                sexlines2[i] = 'SEEING_FWHM     '+str(fwhmas)+'  # stellar FWHM in arcsec' 
        # MASK/WEIGHT file 
        if maskfile is not None:
            sexlines2 = np.append(sexlines2,'WEIGHT_IMAGE  '+maskfile)
            sexlines2 = np.append(sexlines2,'WEIGHT_TYPE   MAP_WEIGHT')
        # Write the file
        dln.writelines(sexconfigfile,sexlines2)
          
    logger.info('Using '+detectprog.upper()+' for Source Detection')
     
    # Loop to find all the objects the first time, bright star problems 
    flag = 0 
    nals = 0 
    als = []
    count = 1 
    allsex = None
    while (flag==0):
        logger.info('--Iteration '+str(count)+'--')
        nlastals = nals 
                  
        ####################### 
        # Detect New Sources 
        ####################### 
         
        #---------------------------- 
        # SExtractor Source Detection 
        #---------------------------- 
        if (detectprog == 'sextractor'):
             
            #------------------------------------- 
            # Initial run of SExtractor -- output sex.cat 
            #!sex allf.fits -c default.sex 
            # Run sextractor on star subtracted file 
            logger.info('Running SExtractor')
            if os.path.exists(catfile): os.remove(catfile) # delete sextractor catalog file 
            out = subprocess.check_output(['sex',subfile,'-c',sexconfigfile],shell=False,stderr=subprocess.STDOUT)
            if os.path.exists(catfile)==False or os.path.getsize(catfile)==0:  # no output file 
                raise ValueError('Error when running SExtractor')
             
            #------------------------------------- 
            # Load sextractor output file 
            # default.param specifies the output columns 
            if utils.file_isfits(catfile)==False:
                #fields = dln.readlines('default.param')
                #fields = np.char.array(fields)
                #gd, = np.where((fields.startswith('#')==False) & (fields.strip() != ''))
                #fields = fields[gd] 
                #fields = fields.replace('(1)','')
                #fields = ['ID','X','Y','MAG','ERR','FLAGS','STAR']
                sex = Table.read(catfile,format='ascii')
            else:
                sex = Table.read(catfile,format='fits')
                if len(sex.colnames)==0:
                    sex = Table.read(catfile,2,format='fits')
            nsex = len(sex)
            sex['ndetiter'] = 0    # add the detection iteration 
            sex['ndetiter'] = count 
            logger.info('SExtractor found '+str(nsex)+' sources')

            #------------------------------------- 
            # Make coordinate input file for ALLSTAR 
             
            # Get ALS header
            with open(base+'.als','r') as f:
                line1 = f.readline().replace("\n","")
                line2 = f.readline().replace("\n","")                        
            head = [line1,line2]
             
            # Concatenate with the ALS file to make a combined 
            # list of sources 
            coofile = base+'_all.coo'
            if os.path.exists(coofile): os.remove(coofile) # delete final coordinate file 
            nals = len(als) 
             
            # ALS file exists, concatenate 
            if nals > 0: 
                maxid = np.max(als['id'])# The last ALS ID                  
                nsex = len(sex) 
                # Making ALS structure for new SEX sources 
                dt = [('id',int),('x',float),('y',float),('mag',float),('err',float),
                      ('sky',float),('niter',int),('chi',float),('sharp',float)]
                sex2 = np.zeros(nsex,dtype=np.dtype(dt))
                sex2['id'] = sex['number'] + maxid   # offset the IDs, sex IDs start at 1 
                sex2['x'] = sex['x_image']
                sex2['y'] = sex['y_image']
                sex2['mag'] = sex['mag_aper']
                sex2['err'] = sex['magerr_aper']
                sex2['sky'] = np.median(als['sky']) 
                sex2['niter'] = 1 
                sex2['chi'] = 1.0 
                sex2['sharp'] = 0.0 
                 
                # Add to final sextractor file 
                # This will be a file that has all sources SExtractor found 
                # with their IDs properly offset so they can be matched to the 
                # final ALS file 
                sex['number'] += maxid    # offset the IDs 
                allsex = np.hstack((allsex,sex))
                 
                # Concatenate 
                allals = np.hstack((als,sex2))
                nall = len(allals) 
                 
                # Write to file 
                io.writeals(coofile,allals,alshead)
                 
            # First time, no ALS file yet 
            else:  
                # Making ALS structure for new SEX sources 
                dt = [('id',int),('x',float),('y',float),('mag',float),('err',float),
                      ('sky',float),('niter',int),('chi',float),('sharp',float)]
                sex2 = np.zeros(nsex,dtype=np.dtype(dt))
                sex2['id'] = sex['NUMBER']
                sex2['x'] = sex['X_IMAGE']
                sex2['y'] = sex['Y_IMAGE']
                sex2['mag'] = sex['MAG_APER'] 
                sex2['err'] = sex['MAGERR_APER']
                sex2['sky'] = 0.0 
                sex2['niter'] = 1 
                sex2['chi'] = 1.0 
                sex2['sharp'] = 0.0 
                 
                # Write to file 
                io.writeals(coofile,sex2,head)
                 
                # Initialize the structure of all SExtractor detections 
                allsex = np.copy(sex)
             
             
        #------------------------- 
        # DAOPHOT Source Detection 
        #------------------------- 
        else: 
             
            # Remove the SExtractor file.  If it exists it's an old one.
            if os.path.exists(sexfile): os.remove(sexfile)
             
            #------------------------------------- 
            # Run DAOPHOT/FIND and PHOT on star subtracted file 
            logger.info('Running DAOPHOT')
             
            # Copy the OPT file from _comb.opt to _comb_sub.opt 
            if os.path.exists(subbase+'.opt')==False:
                if os.path.exists(subbase+'.opt'): os.remove(subbase+'.opt')
                shutil.copyfile(base+'.opt',subbase+'.opt')
             
            # Copy the .opt file daophot.opt 
            if os.path.exists('daophot.opt')==False:
                if os.path.exists('daophot.opt'): os.remove('daophot.opt')
                shutil.copyfile(subbase+'.opt','daophot.opt')
             
            # Make temporary photo.opt file
            tid,tphotofile = tempfile.mkstemp(prefix="photo",dir='.')
            tphotofile = os.path.basename(tphotofile)
            photlines = []
            photlines += ['A1 = %7.4f' % (fwhm*3)]
            photlines += ['IS = 45.0000']
            photlines += ['OS = 50.0000']
            dln.writelines(tphotofile,photlines)
                        
            # Make a temporary script to run FIND 
            lines = []
            lines += ['#!/bin/sh']
            lines += ['export image=${1}']
            lines += ['rm ${image}.temp.log      >& /dev/null']
            lines += ['rm ${image}.temp.coo      >& /dev/null']
            lines += ['rm ${image}.temp.ap       >& /dev/null']
            lines += ['daophot << END_DAOPHOT >> ${image}.temp.log']
            lines += ['OPTIONS']
            lines += ['${image}.opt']
            lines += ['']
            lines += ['ATTACH ${image}.fits']
            lines += ['FIND'] 
            lines += ['1,1']
            lines += ['${image}.temp.coo']
            lines += ['y']
            lines += ['PHOTOMETRY']
            lines += [os.path.basename(tphotofile)]
            lines += ['']
            lines += ['${image}.temp.coo'] 
            lines += ['${image}.temp.ap']
            lines += ['EXIT']
            lines += ['END_DAOPHOT']
            tid2,tmpfile = tempfile.mkstemp(prefix="dao",dir='.')
            tmpfile = os.path.basename(tmpfile)
            dln.writelines(tmpfile,lines)
            os.chmod(tmpfile,0o755)
             
            # Run the program
            out = subprocess.check_output(tmpfile+' '+subbase,shell=False,stderr=subprocess.STDOUT)
            if os.path.exists(tmpfile): os.remove(tmpfile) # delete the temporary script 
             
            #------------------------------------- 
            # Load DAOPHOT FIND and PHOT output files
            io.readcoo(subbase+'.temp.coo',coo,coohead)
            io.readaper(subbase+'.temp.ap',aper,aperhead)
            for f in [subbase+'.temp.coo',subbase+'.temp.ap']:
                if os.path.exists(f): os.remove(f)
            ncoo = len(coo) 
            logger.info('DAOPHOT found '+str(ncoo)+' sources')
                          
            # Get ALS header
            with open(base+'.ls','r') as f:
                line1 = f.readline().replace("\n","")
                line2 = f.readline().replace("\n","")                        
            head = [line1,line2,''] 
             
            # Concatenate with the ALS file to make a combined 
            # list of sources 
            coofile = base+'_all.coo'
            if os.path.exists(coofile): os.remove(coofile) # delete final coordinate file 
            nals = len(als) 
            # ALS file exists, concatenate 
            if (nals > 0): 
                maxid = np.max(als['id'])  # The last ALS ID 
                 
                # Make a fake ALS structure 
                #-------------------------- 
                ndao = len(coo) 
                # Making ALS structure for new DAOPHOT sources
                dt = [('id',int),('x',float),('y',float),('mag',float),('err',float),
                      ('sky',float),('niter',int),('chi',float),('sharp',float)]
                dao = np.zeros(ndao,dtype=np.dtype(dt))
                dao['id'] = coo['id'] + maxid  # offset the IDs, sex IDs start at 1 
                dao['x'] = coo['x'] 
                dao['y'] = coo['y'] 
                dao['mag'] = aper['mag'][0] 
                dao['err'] = aper['err'][0] 
                dao['sky'] = aper['sky']
                dao['niter'] = 1 
                dao['chi'] = 1.0 
                dao['sharp'] = 0.0 
                 
                # Concatenate 
                allals = np.hstack((als,dao))
                nall = len(allals) 
                 
                # Write to file 
                io.writeals(coofile,allals,alshead)

                 
            # First time, no ALS file yet 
            else: 
                                  
                # Make a fake ALS structure 
                #-------------------------- 
                ndao = len(coo) 
                # Making ALS structure for new DAOPHOT sources 
                dt = [('id',int),('x',float),('y',float),('mag',float),('err',float),
                      ('sky',float),('niter',int),('chi',float),('sharp',float)]
                dao = np.zeros(ndao,dtype=np.dtype(dt))
                dao['id'] = coo['id'] 
                dao['x'] = coo['x'] 
                dao['y'] = coo['y'] 
                dao['mag'] = aper['mag'][0] 
                dao['err'] = aper['err'][0] 
                dao['sky'] = aper['sky'] 
                dao['niter'] = 1 
                dao['chi'] = 1.0 
                dao['sharp'] = 0.0 
                 
                # Use the DAOPHOT FIND/PHOT data 
                io.writeals(coofile,dao,head)
         
        #----------------------------------------------------- 
        # Run ALLSTAR on all sources found so far 
        # on original frame 
        logger.info('Running ALLSTAR')
        if os.path.exists(subfile): os.remove(subfile)  # delete fits subfile
        if os.path.exists(subbase+'.als'): os.remove(subbase+'.als')  # delete als output file 
        # Sometimes the filenames get too long for daophot/allstar, 
        # use temporary files and symlinks
        tid3,tbase = tempfile.mkstemp(prefix="base",dir='.')  # create base, leave so other processes won't take it
        tbase = os.path.basename(tbase)
        tid4,tsub = tempfile.mkstemp(prefix="sub",dir='.')    # create sub base, leave so other processes won't take it
        tsub = os.path.basename(tsub)
        tbasefits = tbase+'.fits'
        if os.path.exists(tbasefits): os.remove(tbasefits)
        os.symlink(base+'.fits',tbasefits)
        tbasepsf = tbase+'.psf'
        if os.path.exists(tbasepsf): os.remove(tbasepsf)
        os.symlink(base+'.psf',tbasepsf) 
        tid5,tcoofile = tempfile.mkstemp(prefix="coo",suffix='.coo',dir='.')
        tcoofile = os.path.basename(tcoofile)
        if os.path.exists(tcoofile): os.remove(tcoofile)
        os.symlink(coofile,tcoofile)
        tsubals = tsub+'.als'
        if os.path.exists(tsubals): os.remove(tsubals)
        tsubfits = tsub+'.fits'
        if os.path.exists(tsubfits): os.remove(tsubfits)
        # Make Input file 
        cmd = []
        alsoptlines = dln.readlines(base+'.als.opt') # make sure we use the right als options 
        cmd += alsoptlines
        cmd += [tbasefits] # image file (symlink) 
        cmd += [tbasepsf]  # psf file (symlink)
        cmd += [tcoofile]  # coordinate file (symlink) 
        cmd += [tsubals]   # new als file (temporary file) 
        cmd += [tsubfits]  # subtracted fits file (temporary file) 
        cmdid,cmdfile = tempfile.mkstemp(prefix='temp',dir='.')
        cmdfile = os.path.basename(cmdfile)
        dln.writelines(cmdfile,cmd)
        out = subprocess.check_output('allstar < '+cmdfile,shell=True,stderr=subprocess.STDOUT)
        out = out.decode().split('\n')
        for o in out: print(o)
        # Copy and clean up
        if os.path.exists(subbase+'.als'): os.remove(subbase+'.als')
        shutil.move(tsubals,subbase+'.als')
        if os.path.exists(subbase+'.fits'): os.remove(subbase+'.fits')
        shutil.move(tsubfits,subbase+'.fits')
        for f in [cmdfile,tbase,tsub,tbasefits,tbasepsf,tcoofile]:
            if os.path.exists(f): os.remove(f)
         
        # Load ALS file
        als,alshead = io.readals(subbase+'.als')
        nals = len(als) 
        logger.info('ALLSTAR found '+str(nals)+' sources')
         
        # How many new stars 
        nnew = nals-nlastals 
        logger.info(str(nnew)+' new stars found')
        if nnew < 10: 
            flag = 1 
        if nnew < int(np.round(nals*0.01)):  # more than 1% of total 
            flag = 1 
        if count >= maxiter: 
            flag = 1 
         
        # Increment counter 
        count += 1 
         
     
    # Write out the SExtractor catalog 
    Table(allsex).write(sexfile,format='fits',overwrite=True)
     
    # The X/Y pixel coordinates should already be in the reference image 
    # coordinate system 
     
    # Need to offset the final X/Y pixel coordinates for XOFF/YOFF 
    logger.info('Applying offsets: Xoff=%.2f Yoff=%.2f' % (xoff,yoff))
    als,alshead = io.readals(subbase+'.als')
    als['x'] += xoff 
    als['y'] += yoff
    io.writeals(base+'_allf.als',als,alshead)
    logger.info('Final ALS file = '+base+'_allf.als')
    if os.path.exists(sexfile):
        logger.info('Final SExtractor file = '+sexfile)

    return als


def makemag(tfrfile,outfile,nowrite=False):
    """
    This combines the ALLFRAME alf photometry output files 

    Parameters
    ----------
    tfrfile : str
       Filename of the TFR file.
    outfile : str, optional
       The output filename to write the final combined photometry to.
    nowrite : boolean, optional
       Do not write final table to a file.  Default is False.

    Returns
    -------
    phot : table
      The final table of the combined alf photometry.

    Example
    -------
    phot = makemag(tfrfile)

    """

    if os.path.exists(tfrfile)==False:
        raise ValueError(tfrfile+' NOT FOUND')
     
    # GETTING the number of files and filenames 
    #------------------------------------------ 
    nfiles = 0
    lines = []
    line = ''
    with open(tfrfile,'r') as f:
        line = f.readline()
        while line.strip().startswith('===')==False:
            lines += [line]
            nfiles += 1
            line = f.readline()
     
    # No files 
    if nfiles < 1: 
        print('No files in '+tfrfile)
        return None
     
    # Getting filenames 
    lines = np.char.array(lines).strip()
    files = [l.split()[0] for l in lines]
     
    # Loading the entire TFR file 
    #---------------------------- 
    nlines = dln.numlines(tfrfile)
    with open(tfrfile,'r') as f:
        for i in range(nfiles+1):
            f.readline()
        nrow = nlines-nfiles-1
        ncol = 3+nfiles

        # Only 9-18 per line, the rest wraps on multiple lines 
        # but ID/X/Y stay the same 
        ids = np.zeros(nrow,int) 
        x = np.zeros(nrow,float) 
        y = np.zeros(nrow,float) 
        num = np.zeros((nrow,nfiles),int)-1 
        flag = 0 
        count = 0
        istar = 0 
        line = f.readline()
        while line!='':
            arr = line.split()
            id1 = int(arr[0]) 
            x1 = float(arr[1]) 
            y1 = float(arr[2]) 
            num1 = np.array(arr[3:]).astype(int)
            nnum1 = len(num1) 
         
            # First line 
            if count == 0: 
                ids[istar] = id1 
                x[istar] = x1 
                y[istar] = y1 
                num[istar,0:nnum1] = num1 
                numcount = nnum1 
            # Second or later lines 
            else: 
                # Same star, wrapped line 
                if id1 == ids[istar]: 
                    num[istar,numcount:numcount+nnum1] = num1 
                    numcount += nnum1 
                # New star 
                else: 
                    istar += 1
                    ids[istar] = id1 
                    x[istar] = x1 
                    y[istar] = y1 
                    num[istar,0:nnum1] = num1 
                    numcount = nnum1 

            line = f.readline()         
            count += 1

    # Trim excess lines 
    nstars = istar+1 
    ids = ids[0:nstars] 
    x = x[0:nstars] 
    y = y[0:nstars] 
    num = num[0:nstars,:] 
     
     
    magarr = np.zeros((nstars,nfiles),float)+99.9999 
    magerrarr = np.zeros((nstars,nfiles),float)+9.9999 
    skyarr = np.zeros((nstars,nfiles),float)
    iterarr = np.zeros((nstars,nfiles),float)
    chiarr = np.zeros((nstars,nfiles),float)+np.nan  # NANs are ignored by MEDIAN 
    sharparr = np.zeros((nstars,nfiles),float)+np.nan 
    countarr = np.zeros(nstars,int)  # how many good mags for this stars 
              
    # Loop through the ALF files 
    #--------------------------- 
    for i in range(nfiles):          
        if os.path.exists(files[i]) == False:
            print(files[i]+' NOT FOUND')
            return 

        # Load the ALF file 
        alf,alfhead = io.readals(files[i])
        nalf = len(alf)
         
        ind = num[:,i] 
        gd, = np.where(ind > 0) 
        bd, = np.where(ind == 0) 
         
        alfind = ind[gd]-1   # python indices 
         
        if len(gd) > 0: 
            magarr[gd,i] = alf['mag'][alfind]
            magerrarr[gd,i] = alf['err'][alfind]
            skyarr[gd,i] = alf['sky'][alfind] 
            iterarr[gd,i] = alf['niter'][alfind]
            chiarr[gd,i] = alf['chi'][alfind]
            sharparr[gd,i] = alf['sharp'][alfind]
            countarr[gd] += 1
     
    # Calculating the Median CHI and SHARP 
    if nfiles > 1: 
        chi = np.nanmedian(chiarr,axis=1) 
        sharp = np.nanmedian(sharparr,axis=1)
        # The fortran makemag code used the mean 
        chimean = np.nansum(chiarr,axis=1)/countarr 
        sharpmean = np.nansum(sharparr,axis=1)/countarr 
    # only 1 file 
    else: 
        chi = chiarr 
        sharp = sharparr 
        chimean = chiarr 
        sharpmean = sharparr 
     
    # Creating mag/magerr output array 
    magoutarr = np.zeros((nstars,nfiles*2),float)
    magoutarr[:,np.arange(nfiles)*2] = magarr 
    magoutarr[:,np.arange(nfiles)*2+1] = magerrarr 
     
    # Create PHOT structure 
    dt = [('id',int),('x',float),('y',float)]
    for i in range(nfiles): 
        dt += [('mag'+str(i+1),float),('mag'+str(i+1)+'err',float)]
    dt += [('chiarr',(float,nfiles)),('sharparr',(float,nfiles)),('chi',float),('sharp',float)]
    phot = np.zeros(nstars,dtype=np.dtype(dt))
    phot['id'] = ids
    phot['x'] = x 
    phot['y'] = y 
    for i in range(nfiles): 
        phot['mag'+str(i+1)] = magarr[:,i]
        phot['mag'+str(i+1)+'err'] = magerrarr[:,i]
    phot['chiarr'] = chiarr
    phot['sharparr'] = sharparr
    phot['chi'] = chimean 
    phot['sharp'] = sharpmean 
          
    # Print the output 
    if not nowrite:
        print('Writing output to '+outfile)
        with open(outfile,'w') as f:
            for i in range(nstars): 
                form = '%1s%8d%9.3f%9.3f'+(nfiles*2+2)*'%9.4f'+'\n'
                data = ('',ids[i],x[i],y[i],*magoutarr[i,:],chimean[i],sharpmean[i])
                f.write(form % data)

    return phot



def cleanup(mchbase,files,fpack,mchdir,workdir,tempdir):
    """ Clean up """
    
    # Delete temporarily funpacked files 
    bdfpack, = np.where(fpack)
    if len(bdfpack)>0:
        for f in np.array(files)[bdfpack]:
            base1,ext = os.path.splitext(os.path.basename(f))
            if os.path.exists(base1+'.fits'): os.remove(base1+'.fits')
         
    #------------------------------------------------ 
    # Working in temporary directory, copy files back 
    #------------------------------------------------ 
    if len(workdir) > 0: 
        logger.info('Copying files back to original directory')
        # Make sure all files are writeable
        files0 = glob(tempdir+'/*')
        #files0 = file_search(tempdir+'/*',count=nfiles0,/match_initial_dot)
        nfiles0 = len(files0)
        if nfiles0 > 0:
            for f in files0:
                os.chmod(f,0o755)
        # Delete some files
        for e in ['.mch','.raw.','.tfr']:
            if os.path.exists(tempdir+'/'+mchbase+e): os.remove(tempdir+'/'+mchbase+e)
        for f in files:
            base1,ext = os.path.splitext(os.path.basename(f))
            for e in ['fits','fits.head','opt','als.opt','ap','als','log','psf']:
                if os.path.exists(tempdir+'/'+f+'.'+e): os.remove(tempdir+'/'+f+'.'+e)
            if os.path.exists(tempdir+'/.'+f+'.fits'): os.remove(tempdir+'/.'+f+'.fits')
        if fake:
            for e in ['.weights','.scale','.zero','_comb.psf','_comb.mch']:
                if os.path.exists(tempdir+'/'+mchbase+e): os.remove(tempdir+'/'+mchbase+e)
        # Copy files back
        allfiles = glob(tempdir+'/*')
        nallfiles = len(allfiles)
        for f in allfiles:
            base1 = os.path.basename(f)
            if os.path.exists(mchdir+base1): os.remove(mchdir+base1)
            shutil.copyfile(f,mchdir)
            os.remove(f)  # Delete all temporary files 
        # CD back
        os.chdir(curdir)
        # Delete temporary directory 
        #  leave the base working directory 
        os.remove(tempdir)
 


def allframe(infile,tile=None,setupdir=None,setup=None,scriptsdir=None,
             detectprog='sextractor',logfile=None,finditer=2,irafdir=None,
             satlevel=6e4,nocmbimscale=False,trimcomb=False,usecmn=False,fake=False,
             catformat='ASCII',imager=None,workdir=None,geocoef=None):
    """
    This runs ALLFRAME on images 

    You need to have run daophot, allstar, daomatch and daomaster 
    already.  There needs to be a fits, opt, als.opt, ap and als 
    file for each file in the MCH file.  There also needs to be an 
    associated RAW file for the MCH file. 
    
    Parameters
    ---------- 
    infile : str
       The MCH filename 
    tile : struct, optional
       Information on the tiling to use for the combined 
         image.  If not set then the "original" method is 
         used. 
    finditer : int, optional
       The maximum number of iterations to use when finding 
         sources with SExtractor/ALLSTAR.  The default is 1, 
         and the maximum allowed it 10. 
    detectprog : str, optional
       The program to use to detect sources.  Either 
         'sextractor' or 'daophot'.  'sextractor' is the 
         default.  'daophot' is better for VERY crowded 
         regions.  Default is 'sextractor'.
    nocmbimscale : boolean, optional
       Don't scale the images when combining them.  Not 
         recommended, but the old way of doing it.  Bright 
        stars can be missed this way. 
    combtrim : boolean, optional
       Trim the combined images to the overlapping region. 
         This used to be the default, but now the default 
         is to keep the entire original region. 
    setup : dict, optional
       The information contained in the photred.setup file.
    setupdir : str, optional
       The original base directory which contains photred.setup. 
    scriptsdir : str, optional
       The directory that contains all of the necessary scripts. 
    irafdir : str, optional
       The IRAF home directory. 
    logfile : str, optional
        A logfile to print to output to. 
    usecmn : boolean, optional
        Use the common sources file of the reference image. 
    fake : boolean, optional
        Run for artificial star tests. 
    catformat : str, optional
        Catalog format to use: FITS or ASCII.  Default is ASCII.
    imager : struct
        Imager structure with basic information. 
    workdir : str, optional
        Use a temporary working directory with this as the base. 
    geocoef : int, optional
        The number of geometric coefficients to allow for 
          fitting in ALLFRAME. 

    Returns
    -------
    The final allframe output file name is filebase+'.mag' 
    
    Example
    -------
    
    allframe('ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/',finditer=2)
    
    
    By D.Nidever   February 2008 
    Automation of steps and scripts by J.Ostheimer and Rachael Beaton 
    Translated to Python by D. Nidever   April 2022
    """

    # Set up logging to screen and logfile
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    logger = logging.getLogger() 
    while logger.hasHandlers(): # some existing loggers, remove them   
        logger.removeHandler(logger.handlers[0]) 
    logger = logging.getLogger()
    logtime = datetime.now().strftime("%Y%m%d%H%M%S") 
    if logfile is not None:
        if os.path.exists(logfile): os.remove(logfile)
        fileHandler = logging.FileHandler(logfile)
        fileHandler.setFormatter(logFormatter)
        logger.addHandler(fileHandler)
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)
        logger.setLevel(logging.NOTSET)
    else:
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)
        logger.setLevel(logging.NOTSET)      

    # Get the setup information
    if setup is None:
        setup = io.readsetup(setupdir)
        if setup is None:
            raise ValueError('No setup file found')

    # FIND iterations 
    finditer = np.minimum(finditer, 10)  # maximum 10. 
         
    # Getting scripts directory and iraf directory 
    scriptsdir = setup['scriptsdir']
    irafdir = setup['irafdir']
         
    # No irafdir
    if scriptsdir is None:
        raise ValueError('scriptsdir must be input')
         
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if iraf.check(irafdir=irafdir)==False:
        raise ValueError('IRAF TEST FAILED.')

    # Check if the scripts exist in the current directory 
    scripts = ['getpsfnofind.sh','allstar.sh','photo.opt','apcor.opt','lstfilter.py',
               'goodpsf.pro','allframe.opt','default.sex','default.param',
               'default.nnw','default.conv'] 
    nscripts = len(scripts) 
    # Loop through the scripts 
    for i in range(nscripts): 
        exists = os.path.exists(scriptsdir+'/'+scripts[i])
        if exists:
            size = os.stat(scriptsdir+'/'+scripts[i])
        else:
            size = 0
        curexists = os.path.exists(scripts[i])
        if curexists:
            cursize = os.stat(scripts[i])
        else:
            cursize = 0
        # No file
        if exists==False or size==0:
            raise ValueError(scriptsdir+'/'+scripts[i]+' NOT FOUND or EMPTY')
        # Check if the two files are the same size, if not copy it
        if size!=cursize:
            if curexists:
                os.remove(scripts[i])
            shutil.copyfile(scriptsdir+'/'+scripts[i],scripts[i])
     
    # Check that the ALLFRAME program exists
    out = subprocess.check_output(['which','allframe'],shell=False,stderr=subprocess.STDOUT)
    if type(out) is bytes: out = out.decode()
    out = out.strip()
    if os.path.exists(out)==False:
        raise ValueError('No allframe found')
     
    # Combination settings and inputs
    cmborig = False
    if tile is None:
        cmborig = True
    else:
        if type(tile) is not dict and type(tile):
            raise ValueError('tile must be a dictionary')
        if tile.get('type') is None:
            raise ValueError('tile must have type column')
        if tile['type']=='ORIG': 
            cmborig = True
     
    logger.info('')
    logger.info('')
    logger.info('=====================================')
    logger.info('RUNNING ALLFRAME on '+str(infile))
    logger.info('=====================================')
    logger.info('')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                        
     
    # FILENAME 
    mchfile = os.path.basename(infile) 
    mchdir = os.path.dirname(infile) 
    mchbase = os.path.splitext(os.path.basename(infile))[0]

    # CD to the directory
    curdir = os.getcwd()
    os.chdir(mchdir)

    # Check that the mch, als, and opt files exist
    for f in [mchfile,mchbase+'.raw']:
        if os.path.exists(f)==False:
            raise ValuError(f+' NOT FOUND')

     
    ############################################ 
    # CHECK NECESSARY FILES 
     
    # Load the MCH file
    files,trans,magoff = io.readfile(mchfile)
     
    # Check that the fits, als, opt, and psf files exist 
    nfiles = len(files) 
    fpack = np.zeros(nfiles,bool)
    for i in range(nfiles): 
        base,ext = os.path.splitext(os.path.basename(files[i]))
        # Checking FITS file 
        if os.path.exists(base+'.fits')==False and os.path.exists(base+'.fits.fz')==False:
            raise ValueError(base+'.fits/.fits.fz NOT FOUND')
        # Uncompress FPACK FITS files if necessary, temporarily 
        if os.path.exists(base+'.fits')==False and os.path.exists(base+'.fits.fz'):
            fpack[i] = True
            logger.info('Temporarily uncompressing '+base+'.fits.fz')
            out = subprocess.run(['funpack',base+'.fits.fz'],shell=False,stderr=subprocess.STDOUT)
        # Checking OPT file
        for e in ['opt','als.opt','ap','als','log','psf']:
            if os.path.exists(base+'.'+e)==False:
                raise ValueError(base+'.'+e+' NOT FOUND')
        # REMOVE ALF if it exists 
        if os.path.exists(base+'.alf'): 
            os.remove(base+'.alf')
     

    # REMOVE the .mag file if it exists 
    if os.path.exists(mchbase+'.mag'): 
        os.remove(mchbase+'.mag')
 
 
    # FAKE, check that we have all the files that we need
    if fake:
        # In early versions of ALLFRAME _comb.mch was called _shift.mch 
        #  Use new name with link 
        if os.path.exists(mchbase+'_comb.mch')==False and os.path.exists(mchbase+'_shift.mch'): 
            logger.info('Linking '+mchbase+'_comb.mch to '+mchbase+'_shift.mch')
            os.symlink(mchbase+'_shift.mch',mchbase+'_comb.mch')
        # weights, scale, zero, comb_psf, _comb.mch 
        chkfiles = [mchbase+e for e in ['.weights','.scale','.zero','_comb.psf','_comb.mch']]
        chkexists = [os.path.exists(f) for f in chkfiles]
        bdfiles, = np.where(np.array(chkexists)==False)
        if len(bdfiles)>0:
            raise ValueError('FAKE.  Some necessary files not found. '+' '.join(chkfiles[bdfiles]))
 
 
    #------------------------------------ 
    # Using a temporary working directory 
    #------------------------------------ 
    if workdir is not None:
        # Create a temporary directory in WORKDIR 
        if os.path.exists(workdir)==False:
            os.makedirs(workdir)
        tid2,tfile2 = tempfile.mkstemp(prefix="alf",dir=workdir)
        os.chmod(tempdir,0o755)
        logger.info('Working in temporary directory '+tempdir)
        # Copy over the files that we need 
        #  this will copy the contents of symlinks 
        shutil.copyfile(mchbase+['.mch','.raw','.tfr'],tempdir)
        for i in range(nfiles):
            base1,ext = os.path.splitext(os.path.basename(files[i]))
            for e in ['fits','opt','als.opt','ap','als','log','psf']:
                shutil.copyfile(base1+'.'+e,tempdir)
            # Copy resource files and headers if they exist 
            base1 = os.path.basename(files[i],'.als') 
            if os.path.exists('.'+base1+'.fits'): 
                shutil.copyfile('.'+base1+'.fits',tempdir)
            if os.path.exists(base1+'.fits.head'): 
                shutil.copyfile(base1+'.fits.head',tempdir) 
        # Copy files for FAKE 
        if fake:
            for e in ['.weights','.scale','.zero','_comb.psf','_comb.mch']:
                shutil.copyfile(mchbase+e,tempdir)
        # Copy the scripts
        for f in scripts:
            shutil.copyfile(f,tempdir)
        # Go there
        os.chdir(tempdir)

 
    ############################################ 
    # STEP 1: COMBINE the images  
    logger.info('--------------------------')
    logger.info('STEP 1: Combine the images')
    logger.info('--------------------------')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
    # Use the original combine code
    if cmborig:
        maskdatalevel,xoff,yoff = comb.combine_orig(infile,fake=fake,scriptsdir=scriptsdir,logger=logger,
                                                    setup=setup,irafdir=irafdir,satlevel=satlevel,
                                                    nocmbimscale=nocmbimscale,trimcomb=trimcomb)
                                                    
        combmch = mchbase+'.mch' 
    # New combine code 
    else:
        maskdatalevel,fileinfo = comb.combine(infile,tile=tile,fake=fake,scriptsdir=scriptsdir,logger=logger,
                                              setup=setup,irafdir=irafdir,satlevel=satlevel,
                                              nocmbimscale=nocmbimscale,imager=imager)
        xoff = 0.0 
        yoff = 0.0 
        combmch = mchbase+'_comb.mch' 
        combfile = mchbase+'_comb.fits' 
        combweightfile = mchbase+'_comb.mask.fits'
 
    ############################################ 
    # STEP 2: Get PSF for Combined Image 
    logger.info('----------------------------------------')
    logger.info('STEP 2: Getting PSF for Combined Image')
    logger.info('----------------------------------------') 
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
    combbase = os.path.splitext(os.path.basename(combfile))[0]
    if fake==False:
        # Make .opt files, set saturation just below the mask data level
        mkopt.mkopt(combfile,va=1,hilimit=maskdatalevel-1000)
        #MKOPT,combfile,satlevel=maskdatalevel-1000 
        # THIS IS NOW DONE IN ALLFRAME_COMBINE/ALLFRAME_COMBINE_ORIG.PRO ABOVE 
        # Using CMN.LST of reference frame if it exists 
        #if file_test(mchbase+'.cmn.lst') and keyword_set(usecmn) then begin 
        #  print,'Using reference image COMMON SOURCE file' 
        #  shutil.copyfile(mchbase+'.cmn.lst',mchbase+'_comb.cmn.lst',/over,/allow 
        #endif
        try:
            pfile = getpsf.getpsf(combbase)
        except:
            traceback.print_exc()
            cleanup(mchbase,files,fpack,mchdir,workdir,tempdir)  # cleanup
            return            
 
    # FAKE, use existing comb.psf file 
    else: 
        logger.info('Using existing '+combbase+'.psf file and running ALLSTAR.')
        os.chmod('allstar.sh',0o755)
        out = subprocess.run('./allstar.sh '+combbase,shell=True,stderr=subprocess.STDOUT)
        logger.info(' ')
 
 
    ############################################ 
    # STEP 3: Run allframe prep 
    #  This iteratively runs SExtractor on the combined image 
    #  This can take a while. 
    logger.info('--------------------------------')
    logger.info('STEP 3: Running allframe prep')
    logger.info('--------------------------------') 
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
    # Make sure we have an allstar.opt file
    if os.path.exists('allstar.opt')==False:
        shutil.copyfile(base[0]+'.als.opt','allstar.opt')
    try:
        als = allfprep(combfile,xoff,yoff,logfile=logfile,
                       detectprog=detectprog,scriptsdir=scriptsdir,maxiter=finditer,
                       maskfile=combweightfile)
    except:
        traceback.print_exc()
        cleanup(mchbase,files,fpack,mchdir,workdir,tempdir)  # cleanup
        return        

 
    ############################################ 
    # STEP 4: Running ALLFRAME 
    logger.info('----------------------------')
    logger.info('STEP 4: Running ALLFRAME')
    logger.info('----------------------------')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
                
    # What we need 
    # allf.mag     List of coordinates made by allfprep 
    # allf.mch     List of transformations 
    # allframe.opt 
    # obj????.psf 
    # obj????.als 
    # obj????.fits 
 
    # Delete any temporary ALLFRAME files from possible 
    # previous runs of allframe.  Otherwise ALLFRAME 
    # will start from where it left off. 
    # For each ALLFRAME run there is: 
    #  mchbasename+'.bck' 
    #  mchbasename+'.nmg' 
    #  mchbasename+'.tfr' 
    # For each file in the MCH file there are: 
    #  filebasename+'.alf' 
    #  filebasename+'j.fits' 
    #  filebasename+'k.fits' 
    if os.path.exists(mchbase+'.tfr'):  # copy original
        if os.path.exists(mchbase+'.tfr.orig'): os.remove(mchbase+'.tfr.orig')
        shutil.copyfile(mchbase+'.tfr',mchbase+'.tfr.orig')
    if os.path.exists(mchbase+'.bck'): os.remove(mchbase+'.bck')
    for e in ['.bck','.nmg','.tfr','j.fits','k.fits','.alf']:
        if os.path.exists(mchbase+e): os.remove(mchbase+e)
     
    # Sometimes the filenames get too long for allframe, 
    # use temporary files and symlinks
    tid,tbase = tempfile.mkstemp(prefix="allf",dir='.') # create base, leave so other processes won't take it         
    tmch = tbase+'.mch'
    if os.path.exists(tmch): os.remove(tmch)
    os.symlink(combmch,tmch)
    tals = tbase+'.als'
    if os.path.exists(tals): os.remove(tals)
    os.symlink(mchbase+'_comb_allf.als',tals)
     
    # Geometric coefficients 
    if geocoef is not None:
        logger.info('Modifying ALLFRAME Geometric Coefficients to '+str(geocoef))
        if os.path.exists('allframe.opt'):
            optlines = dln.readlines('allframe.opt')
            opt = io.readopt('allframe.opt')
            g = dln.grep(optlines,'GE = ',index=True)
            if len(g) > 0: 
                optlines[g[0]] = 'GE = '+str(geocoef)
            else:
                optlines.append('GE = '+str(geocoef))
            dln.writelines('allframe.opt',optlines)
         
    # Make input file
    cmd = []
    if geocoef is not None: # geometric coefficients
        cmd += ['ge='+str(geocoef)+'\n']
    cmd += ['    \n']
    cmd += [tmch+'\n']
    cmd += [tals+'\n']
    cmd += ['    \n']
    cmdid,cmdfile = tempfile.mkstemp(prefix='temp',dir='.')
    dln.writelines(cmdfile,cmd)
         
    out = subprocess.run('allframe < '+cmdfile,shell=True,stderr=subprocess.STDOUT)
         
    # Rename tfr and nmg files
    if os.path.exists(combbase+'.tfr'): os.remove(combbase+'.tfr')
    shutil.move(tbase+'.tfr',combbase+'.tfr')
    if os.path.exists(combbase+'.nmg'): os.remove(combbase+'.nmg')
    shutil.move(tbase+'.nmg',combbase+'.nmg')
    if os.path.exists(cmdfile): os.remove(cmdfile)
    for f in files:
        base1,ext = os.path.splitext(os.path.basename(f))
        if os.path.exists(base1+'j.fits'): os.remove(base1+'.fits') # delete subtracted images 
        if os.path.exists(tbase): os.remove(tbase) # delete temporary files and links 
        if os.path.exists(tmch): os.remove(tmch)
        if os.path.exists(tals): os.remove(tals)

    import pdb; pdb.set_trace()
         
    ############################################ 
    # STEP 5: Combine photometry with MAKEMAG 
    # This combines the photometry in the N alf files 
    # and averages chi and sharp 
    logger.info('--------------------------') 
    logger.info('STEP 5: Running MAKEMAG')
    logger.info('--------------------------')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                                            
                                      
    #FILE_COPY,scriptsdir+'makemag','.',/overwrite 
    if os.path.exists(mchbase+'.makemag'): os.remove(mchbase+'.makemag')
         
    # Make input file 
    #magfile = mchbase+'.makemag' 
    #undefine,cmd 
    #push,cmd,mchbase+'.tfr'           ; final tfr file 
    #push,cmd,strtrim(nfiles,2)+',0'   ; nfiles, offset 
    #push,cmd,magfile                  ; final output file 
    #push,cmd,'2'                      ; do not renumber 
    #cmdfile = maketemp('temp','.inp') 
    #cmdfile = MKTEMP('temp') 
    #WRITELINE,cmdfile,cmd 
    #SPAWN,'./makemag < '+cmdfile 
    #FILE_DELETE,cmdfile,/allow        ; delete temporary input file 
         
    magfile = mchbase+'.makemag' 
    # With the new combined files we are using _comb.mch 
    #   and _comb.tfr, 10/23/16 
    # combmch = FILEBASE.mch        ORIG 
    # combmch = FILEBASE_comb.mch   NEW 
    # The tfr file will have the same name but with .tfr 
    try:
        makemag(combbase+'.tfr',magfile)
    except:
        traceback.print_exc()
        cleanup(mchbase,files,fpack,mchdir,workdir,tempdir)  # cleanup
        return
         
    # Prepend the ALF header to the makemag file
    with open(os.path.splitext(os.path.basename(files[0]))[0]+'.alf','rb') as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()            
    head = [line1,line2,line3]
    # Prepend header
    with open(magfile, 'r+') as f:
        content = f.read()
        f.seek(0)
        f.write(head + content)
         
    import pdb; pdb.set_trace()
         
    ####################################################### 
    # STEP 6: Adding SExtractor information 
    logger.info('----------------------------------------') 
    logger.info('STEP 6: Adding SExtractor information')
    logger.info('----------------------------------------')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                                            
                                      
    # combfile_allf.sex can be matched to the makemag file using IDs 
    # Load the SExtractor file 
    sexfile = combbase+'_allf.sex' 
    if os.path.exists(sexfile): 
             
        #------------------------------------- 
        # Load sextractor output file 
        # default.param specifies the output columns 
        if file_isfits(sexfile)==False:
            fields = dln.readlines('default.param')
            #gd , = np.where(strmid(fields,0,1) != '#' and strtrim(fields,2) ne '',ngd) 
            #fields = fields[gd]
            sex = Table.read(sexfile)
            #sex = IMPORTASCII(sexfile,fieldnames=fields,/noprint) 
        else:
            sex = Table.read(sexfile,1)
            #sex = MRDFITS(sexfile,1,/silent) 
            if n_tags(sex) == 1: 
                sex = Table.read(sexfile,2) 
        nsex = len(sex)
             
        # Load the MAKEMAG file 
        mag,alfhead = io.readraw(mchbase+'.makemag')
        nmag = len(mag) 
             
        # Match them with IDs
        ind1,ind2 = dln.match(mag['id'],sex['number'])
        count = len(ind1)
             
        # Add SExtractor information to mag file 
        sextags = sex.colnames
        # New columns 
        newcols = ['FLAGS','CLASS_STAR','MAG_AUTO','MAGERR_AUTO','BACKGROUND','THRESHOLD','ISOAREA_WORLD',
                   'A_WORLD','B_WORLD','THETA_WORLD','ELLIPTICITY','FWHM_WORLD'] 
        newname = ['FLAG','PROB','MAG_AUTO','MAGERR_AUTO','BACKGROUND','THRESHOLD','ISOAREA',
                   'ASEMI','BSEMI','THETA','ELLIPTICITY','FWHM'] 
        for k in range(len(newcols)):
            if newcols[k] in sex.colnames:
            #colind, = np.where(sextags == newcols[k],ncolind) 
            #if len(colind) > 0:
                mag[newname[k]][ind1] = sex[newcols[k]][ind2]
                #add_tag,mag,newname[k],fix('',type=size(sex[0].(colind),/type)),mag 
                #mag[ind1].(n_tags(mag)-1) = sex[ind2].(colind) 
                # Convert to arcsec 
                if newcols[k]=='A_WORLD' or newcols[k]=='B_WORLD' or newcols[k]=='FWHM_WORLD' : 
                    mag[newname[k]][ind1] *= 3600 
             
        if nind < nmag: 
            logger.info('DID NOT MATCH ALL THE STARS!')
             
             
        # Write the final output file 
        #---------------------------- 
        finalfile = mchbase+'.mag' 
        # FITS has a limit 999 columns/fields for binary tables, use ASCII 
        # if over that limit 
        if catformat == 'FITS' and n_tags(mag) > 999 : 
            logger.info('Cannot use FITS output because number of columns>999.  Using ASCII instead')
        if (catformat == 'FITS') and (len(mag.colnames)<1000):
            if os.path.exists(finalfile+'.fits'): os.remove(finalfile+'.fits')
            mag.writeto(finalfile+'.fits',overwrite=True)
            if os.path.exists(finalfile): os.remove(finalfile)
            shutil.move(finalfile+'.fits',finalfiel)
        else:  # ASCII
            mag.writeto(finalfile,overwrite=True)                
            #PRINTSTR,mag,finalfile,/silent 
             
    # DAOPHOT 
    else: 
        # No SExtractor information, just copy .makemag to .mag 
        finalfile = mchbase+'.mag' 
        if catformat == 'FITS' and n_tags(mag) > 999 : 
            logger.info('Cannot use FITS output because number of columns>999.  Using ASCII instead')
        if (catformat == 'FITS') and (n_tags(mag) < 1000):
            mag,alfhead = io.readraw(mchbase+'.makemag')
            if os.path.exists(finalfile+'.fits'): os.remove(finalfile+'.fits')
            mag.writeto(finalfile+'.fits',overwrite=True)
            if os.path.exists(finalfile): os.remove(finalfile)
            shutil.move(finalfile+'.fits',finalfiel)               
        else:  # ASCII
            if os.path.exists(finalfile): os.remove(finalfile)
            shutil.copyfile(mchbase+'.makemag',finalfile)
         
    logger.info('FINAL ALLFRAME file = '+finalfile)
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))

    # Clean up
    cleanup(mchbase,files,fpack,mchdir,workdir,tempdir)
         

    
