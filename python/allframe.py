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
from astropy.io import fits,ascii
from astropy.table import Table
from dlnpyutils import utils as dln

def loadsetup():
    pass

def loadmch():
    pass

def loadmch(mchfile):
    """
    This loads a DAOMATCH/DAOMASTER mch file 
     
    Parameters
    ----------
    mchfile  The MCH filename 
     
    Returns
    -------
    files    The list of files in the MCH file 
    trans    The transformation equations in a Nfilesx6 array. 
    magoff   The magnitude offsets and errors array. 
     
    Example
    -------

    files,trans,magoff = loadmch('ccd1001.mch')
     
    By D.Nidever   February 2008 
    Translated to Python by D. Nidever,  April 2022
    """ 
     
    count = 0 
     
    # Test the file 
    if os.path.exists(mchfile)==False:
        raise ValueError(mchfile+' NOT FOUND')

    # Read in the file
    lines = dln.readlines(mchfile)

    # Creating the trans array
    lines2 = [re.sub("'","",l) for l in lines]
    #lines2 = repstr(lines,"'",'') 
    nlines = len(lines)

    arr = lines2[0].split()
    ntrans = len(arr)-3  # first line is the samee, last two are mag offsets 
 
    # Getting the file names
    arr2 = [l.split() for l in lines2]
    files = [a[0] for a in arr2]
 
    # Initializing the array
    trans = np.zeros((nlines,ntrans),float)
    magoff = np.zeros((nlines,2),float)
    # Filling the aray 
    for i in range(nlines):
        arr = lines2[i].split()
        trans[i,:] = arr[1:ntrans+1]
        magoff[i,:] = arr[ntrans+1:]
 
    return files,trans,magoff


def readpar(array,keyword):
    """
    This allows you to get a parameter value from 
    a 2xN array of keyword/value pairs. 
    This is similar to getting keyword values from 
    headers with SXPAR.PRO. 
    
    Parameters
    ----------
    array    A 2xN array of keyword-value pairs 
    keyword  A keyword string for which to return the value. 
            Case insensitive. 
 
    Returns
    -------
    value    The value corresponding to the input keyword 
              is output.  If none is found then '0' is returned. 
              If the keyword exists in array but has not value 
              then an empty string '' is returned. 
    count   The number of parameters found by READPAR.  count=0 
              if the parameter was not found. 

    Example
    -------

    value = readpar(setup,'MOSAIC') 
 
    By D. Nidever    Oct. 2007 
    Translated to Python by D. Nidever,  April 2022
    """
    
    #sz = size(array)
    #if sz[0] != 2 or sz[1] != 2: # must be 2xN 
    #    return None
    #if size(array,/type) != 7: # must be string 
    #    return None 
     
    # Looking for keyword 
    keyword2 = strlowcase(str(keyword[0],2)) 
    keys = reform(array[0,:]) 
    values = reform(array[1,:])
    gd, = dln.grep(keys==keyword2)
    if len(gd)==0:
        return None
    value = str(values[gd[0]])  # returning the first value 
    return value

def file_isfits(filename):
    """
    Check if this file is a FITS file or not. 
 
    Parameters
    ----------
    filename   The name of the file to check. 
 
    Returns
    -------
    return     1 if the file is a FITS file and 0 otherwise. 
 
    Example
    ------
  
    test = file_isfits(filename) 
 
    By D. Nidever, Jan 2019 
    Based partially from is_fits.pro by Dave Bazell 
    Translated to python by D. Nidever, April 2022
    """

    # Does the file exist 
    if os.path.exists(filename)==False:
        return False
     
    # Four possible possibilities: 
    # 1) Regular FITS file (this includes fpacked FITS files) 
    # 2) Gzipped FITS file 
    # 3) ASCII file 
    # 4) Some other binary file 
     
    # Try to open the file normally
    # astropy will catch lots of problems
    #  empty file, corrupted file, not a FITS file
    try:
        hdu = fits.open(filename)
    except:
        return False

    hdu.close()
    return True


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
 


def allframe(infile,tile=None,setupdir=None,scriptsdir=None,detectprog='sextractor',
             logfile=None,finditer=2,irafdir=None,satlevel=6e4,
             nocmbimscale=False,trimcomb=False,usecmn=False,fake=False,
             catformat='ASCII',imager=None,workdir=None,geocoef=None):
    """
    This runs ALLFRAME on images 

    You need to have run daophot, allstar, daomatch and daomaster 
    already.  There needs to be a fits, opt, als.opt, ap and als 
    file for each file in the MCH file.  There also needs to be an 
    associated RAW file for the MCH file. 
    
    Parameters
    ---------- 
    file           The MCH filename 
    tile          Information on the tiling to use for the combined 
                  image.  If not set then the "original" method is 
                  used. 
    finditer      The maximum number of iterations to use when finding 
                  sources with SExtractor/ALLSTAR.  The default is 1, 
                  and the maximum allowed it 10. 
    detectprog    The program to use to detect sources.  Either 
                  'sextractor' or 'daophot'.  'sextractor' is the 
                  default.  'daophot' is better for VERY crowded 
                  regions. 
    nocmbimscale  Don't scale the images when combining them.  Not 
                  recommended, but the old way of doing it.  Bright 
                  stars can be missed this way. 
    combtrim      Trim the combined images to the overlapping region. 
                  This used to be the default, but now the default 
                  is to keep the entire original region. 
    setupdir      The original base directory which contains photred.setup. 
    scriptsdir    The directory that contains all of the necessary scripts. 
    irafdir       The IRAF home directory. 
    logfile       A logfile to print to output to. 
    usecmn        Use the common sources file of the reference image. 
    fake          Run for artificial star tests. 
    catformat     Catalog format to use: FITS or ASCII.  Default is ASCII. 
    imager        Imager structure with basic information. 
    workdir       Use a temporary working directory with this as the base. 
    geocoef       The number of geometric coefficients to allow for 
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
     
    global setup 

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
        rootLogger.addHandler(fileHandler)
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
        setup = loadsetup(setupdir)
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
    check_iraf(iraftest,irafdir=irafdir)
    if iraftest == 0: 
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
     
    # Check that the ALLFRAME program exists
    out = subprocess.check_output(['which','allframe'],shell=False)
    if type(out) is bytes: out = out.decode()
    out = out.strip()
    if os.path.exists(out)==False:
        raise ValueError('No allframe found')
     
    # Combination settings and inputs
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
    logger.info('RUNNING ALLFRAME on ',infile)
    logger.info('=====================================')
    logger.info('')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                        
     
    # FILENAME 
    mchfile = os.path.basename(infile) 
    mchdir = os.path.dirname(infile) 
    mchbase = os.path.basename(infile,'.mch') 

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
    files,trans = loadmch(mchfile)
     
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
            out = subprocess.run(['funpack',base+'.fits.fz'],shell=False)     
        # Checking OPT file
        for e in ['opt','als.opt','ap','als','log','psf']:
            if os.path.exists(base+'.'+e):
                raise ValueError(base+'.'+e+' NOT FOUND')
        # REMOVE ALF if it exists 
        if os.path.exists(base+'.alf'): 
            os.remove(base+'.alf')
     

    # REMOVE the .mag file if it exists 
    if os.path.exists(mchbase+'.mag'): 
        os.remove(mchbase+'.mag')
 
 
    # FAKE, check that we have all the files that we need
    if fakse:
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
    if len(workdir) > 0: 
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
        comb.combine_orig(infile,fake=fake,scriptsdir=scriptsdir,error=error,logfile=logfile,
                           irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,trimcomb=trimcomb,
                           maskdatalevel=maskdatalevel,xoff=xoff,yoff=yoff)
        combmch = mchbase+'.mch' 
    # New combine code 
    else:
        comb.combine(infile,tile=tile,fake=fake,scriptsdir=scriptsdir,error=error,logfile=logfile,
                     irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,
                     maskdatalevel=maskdatalevel,filestr=filestr,imager=imager)
        xoff = 0.0 
        yoff = 0.0 
        combmch = mchbase+'_comb.mch' 
        #  There was an error in combination 
        if len(error) > 0: 
            logger.info(error)
            cleanup(mchbase,files,fpack,mchdir,workdir,tempdir)  # cleanup
            return            
        combfile = mchbase+'_comb.fits' 
        combweightfile = mchbase+'_comb.mask.fits'
 
 
    ############################################ 
    # STEP 2: Get PSF for Combined Image 
    logger.info('----------------------------------------')
    logger.info('STEP 2: Getting PSF for Combined Image')
    logger.info('----------------------------------------') 
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
    combbase = os.path.basename(combfile,'.fits')
    if fake:
        # Make .opt files, set saturation just below the mask data level
        mkopt(combfile,va=1,hilimit=maskdatalevel-1000,error=opterror)
        if len(opterror) > 0: 
            logger.info(opterror)
            cleanup(mchbase,files,fpack,mchdir,workdir,tempdir)  # cleanup
            return
        #MKOPT,combfile,satlevel=maskdatalevel-1000 
        # THIS IS NOW DONE IN ALLFRAME_COMBINE/ALLFRAME_COMBINE_ORIG.PRO ABOVE 
        # Using CMN.LST of reference frame if it exists 
        #if file_test(mchbase+'.cmn.lst') and keyword_set(usecmn) then begin 
        #  print,'Using reference image COMMON SOURCE file' 
        #  shutil.copyfile(mchbase+'.cmn.lst',mchbase+'_comb.cmn.lst',/over,/allow 
        #endif
        getpsf(combbase)
        if len(error)>0:
            cleanup(mchbase,files,fpack,mchdir,workdir,tempdir)  # cleanup
            return            
 
    # FAKE, use existing comb.psf file 
    else: 
        logger.info('Using existing '+combbase+'.psf file and running ALLSTAR.')
        out = subprocess.run('./allstar.sh '+combbase)
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
 
    allfprep(combfile,als,xoff,yoff,logfile=logfile,error=error,
             detectprog=detectprog,scriptsdir=scriptsdir,maxiter=finditer,
             maskfile=combweightfile)
    if len(error) > 0:
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
        if os.path.exists(base+e): os.remove(base+e)
     
    # Sometimes the filenames get too long for allframe, 
    # use temporary files and symlinks
    tid,tbase = tempfile.mkstemp(prefix="allf") # create base, leave so other processes won't take it         
    tmch = tbase+'.mch'
    if os.path.exists(tmch): os.remove(tmch)
    os.symlink(combmch,tmch)
    tals = tbase+'.als'
    if os.path.exists(tals): os.remove(tals)
    os.symlink(mchbase+'_comb_allf.als',tals)
     
    # Geometric coefficients 
    if len(geocoef) > 0: 
        logger.info('Modifying ALLFRAME Geometric Coefficients to '+str(geocoef))
        if os.path.exists('allframe.opt'):
            optlines = dln.readlines('allframe.opt')
            g = dln.grep(optlines,'GE = ')
            if ng > 0: 
                optlines[g[0]] = 'GE = '+str(geocoef)
            else:
                optlines.append('GE = '+str(geocoef))
            dln.writeline('allframe.opt',optlines)
         
        # Make input file
        cmd = ''
        if len(geocoef) > 0: # geometric coefficients
            cmd += 'ge='+str(geocoef)+'\n'
        cmd += '    \n'
        cmd += tmch+'\n'
        cmd += tals+'n'
        cmd += '    \n' 
        cmdfile = MKTEMP('temp')
        dln.writelines(cmdfile,cmd)
         
        out = subprocess.run('allframe < '+cmdfile,shell=True)
         
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
        makemag(combbase+'.tfr',magfile,error=magerror)
        if len(magerror)>0:
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
            mag,alfhead = loadmakemag(mchbase+'.makemag')
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
                mag,alfhead = loadmakemag(mchbase+'.makemag')
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
         

    
