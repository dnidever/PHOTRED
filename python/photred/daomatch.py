#!/usr/bin/env python

import os
import time
import numpy as np
import tempfile
from astropy.wcs import WCS
from astropy.table import Table
from scipy.optimize import curve_fit  
from dlnpyutils import utils as dln,coords
import subprocess
import shutil
from glob import glob
from . import io,utils

def test_trans(trans):
    """
    This function tests if a transformation equation from 
    daomatch is good or not.  The scale should be nearly 1.0 
    and the rotation should be near 0. 
    
    Return value: 
     1  Good 
     0  Bad 
    -1  Fail 
    """

    test = -1 
     
    # The test mainly looks at the rotation/scale values 
    # and not the xoff/yoff values. 
     
    if trans.ndim != 2 or trans.shape[0] != 2 or trans.shape[1] != 6: 
        return -1 
     
     
    xoff = trans[1,0][0] 
    yoff = trans[1,1][0]
    c = trans[1,2][0] 
    e = trans[1,3][0] 
    d = trans[1,4][0]
    f = trans[1,5][0]
    # from ccdpck.txt 
    #              x(1) = A + C*x(n) + E*y(n) 
    #              y(1) = B + D*x(n) + F*y(n) 
     
    # C=F~1 and D=E~0 
    test = 1 
    if abs(c-f) > 0.1: 
        test = 0 
    if abs(d-e) > 0.1: 
        test = 0 
    if abs(c-1.0) > 0.1: 
        test = 0 
    if abs(e) > 0.1: 
        test = 0 
     
    return test 

def daomaster(mchbase):
    """
    Running DAOMASTER 
    """
         
    # DAOMASTER has problems with files that have extra dots in them 
    # (i.e. F1.obj1123_1.mch). 
    # Do everything with a temporary file, then rename the output files 
    # at the end. 
    #tempbase = MAKETEMP('temp','')
    tid,tempbase = tempfile.mkstemp(prefix='temp')
    dln.remove(tempbase,allow=True) # remove empty file
    tempbase = tempbase.replace('.','')  # remove the do
    tempmch += '.mch'
    dln.remove(tempmch,allow=True)
    shutil.copyfile(mchbase+'.mch',tempmch)
         
    # Make the DAOMASTER script 
    #--------------------------
    cmdlines = []
    cmdlines += ['#!/bin/csh']
    cmdlines += ['set input=${1}']
    cmdlines += ['daomaster <<DONE']
    cmdlines += ['${input}.mch']
    cmdlines += ['1,1,1']
    cmdlines += ['99.']
    cmdlines += ['6']
    cmdlines += ['10'] 
    cmdlines += ['5']
    cmdlines += ['4']
    cmdlines += ['3']
    cmdlines += ['2']
    cmdlines += ['1']
    cmdlines += ['1']
    cmdlines += ['1']
    cmdlines += ['1']
    cmdlines += ['0']
    cmdlines += ['y']
    cmdlines += ['n']
    cmdlines += ['n']
    cmdlines += ['y']
    cmdlines += ['']
    cmdlines += ['y'] 
    cmdlines += ['']
    cmdlines += ['']
    cmdlines += ['y'] 
    cmdlines += ['']
    cmdlines += ['n'] 
    cmdlines += ['n']
    cmdlines += ['DONE']
    tid,tempscript = tempfile.mkstemp(prefix='daomaster')  # absolute filename
    dln.writelines(tempscript,cmdlines)
    os.chmod(tempscript,0o755)
         
    # Run DAOMASTER 
    #--------------- 
    cmd2 = tempscript+' '+tempbase
    out2 = subprocess.check_output(cmd2,shell=False)
    #SPAWN,cmd2,out2,errout2 
         
    # Remove temporary DAOMASTER script 
    #----------------------------------- 
    dln.remove(tempscript,allow=True)
                  
    # Rename the outputs 
    #------------------- 
    # MCH file 
    mchfile = glob(tempbase+'.mch') 
    if (len(mchfile) > 0):
        dln.remove(files2[0]+'.mch',allow=True)
        shutil.copyfile(mchfile[0],files2[0]+'.mch')
        dln.remove(mchfile,allow=True)
    else: 
        raise ValueError('NO FINAL MCH FILE')
    # TFR file 
    tfrfile = glob(tempbase+'.tfr') 
    if (len(tfrfile) > 0):
        dln.remove(files2[0]+'.tfr',allow=True)
        shutil.copyfile(tfrfile[0],files2[0]+'.tfr')
        dln.remove(tfrfile,allow=True)
    else: 
        raise ValueError('NO FINAL TFR FILE')
    # RAW file 
    rawfile = glob(tempbase+'.raw')
    if (len(rawfile) > 0):
        dln.remove(files2[0]+'.raw',allow=True)
        shutil.copyfile(rawfile[0],files2[0]+'.raw')
        dln.remove(rawfile,allow=True)
    else: 
        raise ValueError('NO FINAL RAW FILE')

 
#--------------------------------------------------------------- 
 
def daomatch(files,usewcs=True,verbose=True,logfile=None,
             maxshift=5000,fake=False):
    """
    This matches stars and finds the transformations using 
    MATCHSTARS.PRO (originally DAOMATCH was used) and then 
    combining them with DAOMASTER.  INputs need to be ALS files. 
 
    Parameters
    ----------
    files : list
       Array of ALS files,  The first file will be used 
         as the reference.  If /fake set then the first ALS file 
         in "files" should already have an associated MCH file. 
    maxshift : int, optional
       Constraints on the initial X/Y shifts. 
    usewcs : bool, optional
       Use the WCS for initial matching.  This is the default.
    fake : bool, optional
       Run for artificial stars.  The MCH file should be input 
         and daomaster run to create raw/tfr files. 
    verbose : bool, optional
       Verbose output.  Default is True.
 
    Returns
    -------
    An MCH, TFR and RAW file with basename of the first file. 
 
    Example
    -------

    daomatch(['obj1034_1.als','obj1035_1.als','obj1036_1.als'])
 
    Add options to freeze the scale and/or rotation. 
 
    By D. Nidever   December 2006  
    Translated to Python by D. Nidever, April 2022
    """

    t0 = time.time() 

    nfiles = dln.size(files)

    # Logfile 
    if logfile is not None:
        logf = logfile 
    else: 
        logf = -1 
     
    # Only one file, can't match 
    if nfiles == 1: 
        raise ValueError('ONLY ONE FILE INPUT.  No matching, simply creating .mch and .raw file')

    import pdb; pdb.set_trace()
     
    # Current directory
    curdir = os.getcwd()
     
    fdir = os.path.abspath(os.path.dirname(files[0]))
    os.path.chdir(fdir)
     
    files2 = [os.path.splitext(os.path.basename(f))[0] for f in files]
    mchbase = files2[0]
     
    # FAKE, running for artificial star tests 
    if fake:
        # Check that MCH file exists 
        if os.path.exists(mchbase+'.mch') == False:
            raise ValueError('/fake set but '+mchbase+'.mch NOT FOUND')
         
        # Keep backup of original mch file
        dln.remove(mchbase+'.mch.orig',allow=True)
        shutil.copyfile(mchbase+'.mch',mchbase+'.mch.orig')
         
        # Remove the output files
        dln.remove([mchbase+'.raw',mchbase+'.tfr'],allow=True)
                   
        #goto,rundaomaster 
     
    # Remove the output files 
    dln.remove(mchbase+'.mch',allow=True)
    dln.remove(mchbase+'.raw',allow=True)
    dln.remove(mchbase+'.tfr',allow=True)
     
    mchfinal = []
     
    # Check that the reference file exists 
    if os.path.exists(files[0])==False:
        raise ValueError('REFERENCE FILE '+files[0]+' NOT FOUND')
     
    # Load the reference data
    refals,alshead = io.readals(files[0])
     
    # Use WCS 
    if usewcs:
        # Checking WCS of first file 
        fitsfile1 = os.path.splitext(os.path.basename(files[0]))[0]+'.fits' 
        if os.path.exists(fitsfile1) == False: 
            fitsfile1 = os.path.splitext(os.path.basename(files[0]))[0]+'.fits.fz' 
        if os.path.exists(fitsfile1) == False: 
            raise ValueError(fitsfile1+' NOT FOUND. Cannot use WCS for matching')
        else:
            if fitsfile1[-7:] =='fits.fz':
                head1 = io.readfile(fitsfile1,exten=1,header=True)
                # Fix the NAXIS1/2 in the header
                head1['NAXIS1'] = head1['ZNAXIS1']
                head1['NAXIS2'] = head1['ZNAXIS2']                   
            else:
                head1 = io.readfile(fitsfile1,header=True)
            wcs1 = WCS(head1)
            if wcs1.has_celestial==False:
                print(fitsfile1+' has NO WCS.  Cannot use WCS for matching')
     
     
    fmt = '%2s%-30s%1s%10.2f%10.2f%10.5f%10.5f%10.5f%10.5f%10.3f%10.3f'
    data = "'",files[0],"'",0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0
    newline = fmt % data
    mchfinal += [newline]
     
    if verbose: 
        print('Initial Transformations:')
     
    # Printing the first line 
    if verbose: 
        print('%-20s%10.4f%10.4f%12.8f%12.8f%12.8f%12.8f' % (files[0], 0.0, 0.0, 1.0, 0.0, 0.0, 1.0))
     
     
    # Run DAOMATCH for each pair (N-1 times) 
    for i in range(nfiles): 
        # Check that the file exists 
        if os.path.exists(files[i])==False:
            raise ValueError('FILE '+files[i]+' NOT FOUND')
                   
        # Load the current data
        als,alshead = io.readals(files[i])
         
        # Getting FRAD
        headarr = alshead[1].split()
        fread = float(headarr[-1])
         
        # Make CHI, SHARP and ERR cuts here 
         
        # CUTS for REFALS 
        gdref, = np.where((np.abs(refals['sharp'])<1.0) & (refals['chi']<2.0) & (refals['mag']<50.0) & (refals['err']<0.2))
        if (len(gdref) < 100): 
            gdref, = np.where((np.abs(refals['sharp'])<1.5) & (refals['chi']<3.0) & (refals['mag']<50.0) & (refals['err']<0.5))
        if (len(gdref) < 100):
            gdref, = np.where((np.abs(refals['sharp'])<1.5) & (refals['chi']<3.0) & (refals['mag']<50.0) & (refals['err']<1.0))
        if (len(gdref) == 0): 
            raise ValueError('NO good reference stars '+files[0])
        # Cuts for ALS 
        gdals, = np.where((np.abs(als['sharp'])<1.0) & (als['chi']<2.0) & (als['mag']<50.0) & (als['err']<0.2))
        if (len(gdals) < 100): 
            gdals, = np.where((np.abs(als['sharp'])<1.5) & (als['chi']<3.0) & (als['mag']<50.0) & (als['err']<0.5))
        if (len(gdals) < 100): 
            gdals, = np.where((np.abs(als['sharp'])<1.5) & (als['chi']<3.0) & (als['mag']<50.0) & (als['err']<1.0))
        if (len(gdals) == 0): 
            raise ValueError('NO good stars for '+files[i])
         
        # --- Use WCS --- 
        if usewcs and noparams1 >= 1: 
            # Checking WCS of second file 
            fitsfile2 = os.path.splitext(os.path.basename(files[i]))[0]+'.fits' 
            if os.path.exists(fitsfile2) == 0 : 
                fitsfile2 = os.path.splitext(os.path.basename(files[i]))[0]+'.fits.fz' 
            if os.path.exists(fitsfile2) == 0: 
                print(' NOT FOUND. Cannot use WCS for matching')
                continue
            if fitsfile2[-7:]=='fits.fz':
                head2 = io.readfile(fitsfile2,exten=1,header=True)
                # Fix the NAXIS1/2 in the header
                head2['NAXIS1'] = head2['ZNAXIS']
                head2['NAXIS2'] = head2['ZNAXIS2']                              
            else:
                head2 = io.readfile(fitsfile2,header=True)
            wcs2 = WCS(head2)
            if wcs2.has_celestial==False:
                print(fitsfile2+' has NO WCS.  Cannot use WCS for matching')
                continue
             
            # Get coordinates for the stars
            coo1 = wcs1.pixel_to_world(refals['x']-1,refals['y']-1)
            coo2 = wcs2.pixel_to_world(als['x']-1,als['y']-1)                 
            #head_xyad,head1,refals.x-1,refals.y-1,a1,d1,/deg 
            #head_xyad,head2,als.x-1,als.y-1,a2,d2,/deg 

            ind1,ind2,dist = coords.xmatch(a1[gdref],d1[gdref],a2[gdals],d2[gdals],1.0,sphere=False)
            count = len(ind1)
            #SRCMATCH,a1[gdref],d1[gdref],a2[gdals],d2[gdals],1.0,ind1,ind2,count=count,/sph
             
            # If no matches, try with looser cuts 
            if count<3: 
                print('No matches, trying looser cuts')
                gdref, = np.where((refals['mag']<50.0) & (refals['err']<1.0))
                gdals, = np.where((als['mag']<50.0) & (als['err']<1.0))
                ind1,ind2,dist = coords.xmatch(a1[gdref],d1[gdref],a2[gdals],d2[gdals],1.0,sphere=False)
                count = len(ind1)
                #SRCMATCH,a1[gdref],d1[gdref],a2[gdals],d2[gdals],1.0,ind1,ind2,count=count,/sph 
                if count < 3:
                    ind1,ind2,dist = coords.xmatch(a1[gdref],d1[gdref],a2[gdals],d2[gdals],3.0,sphere=False)
                    count = len(ind1)                                  
                    #SRCMATCH,a1[gdref],d1[gdref],a2[gdals],d2[gdals],3.0,ind1,ind2,count=count,/sph 
             
            if count > 0: 
                xdiff = refals['x'][gdref[ind1]]-als['x'][gdals[ind2]]
                ydiff = refals['y'][gdref[ind1]]-als['y'][gdals[ind2]]
                xmed = np.median(xdiff) 
                ymed = np.median(ydiff)
                # Fit rotation with linear fits if enough points 
                if count>1:
                    slp1 = dln.mediqrslope(als['y'][gdals[ind2]],xdiff)  # fit rotation term
                    slp1rms = dln.mad(xdiff-als['y'][gdals[ind2]]*slp1)
                    #coef1 = robust_poly_fitq(als[gdals[ind2]].y,xdiff,1)# fit rotation term 
                    #coef1b = dln_poly_fit(als[gdals[ind2]].y,xdiff,1,measure_errors=xdiff*0+0.1,sigma=coef1err,/bootstrap)
                    slp2 = dln.mediqrsope(als['x'][gdals[ind2]],ydiff)
                    slp2rms = dln.mad(ydiff-als['x'][gdals[ind2]]*slp2)                    
                    #coef2 = robust_poly_fitq(als[gdals[ind2]].x,ydiff,1)# fit rotation term 
                    #coef2b = dln_poly_fit(als[gdals[ind2]].x,ydiff,1,measure_errors=ydiff*0+0.1,sigma=coef2err,/bootstrap) 
                    #theta = mean([-coef1[1],coef2[1]]) 
                    #WMEANERR,[-coef1[1],coef2[1]],[coef1err[1],coef2err[1]],theta,thetaerr 
                    theta = dln.wtmean([-slp1,slp2],[slp1rms,slp2rms])
                    
                    # [xoff, yoff, cos(th), sin(th), -sin(th), cos(th)] 
                    trans = np.array([xmed, ymed, 1.0-theta**2, theta, -theta, 1.0-theta**2])
                    # Adjust Xoff, Yoff with this transformation 
                    xout,yout = utils.trans_coo(als['x'][gdals[ind2]],als['y'][gdals[ind2]],trans) 
                    trans[0] += np.median(refals['x'][gdref[ind1]]-xout) 
                    trans[1] += np.median(refals['y'][gdref[ind1]]-yout) 
                else:
                    trans = np.array([xmed, ymed, 1.0, 0.0, 0.0, 1.0])
                     
                # Fit full six parameters if there are enough stars 
                if count > 10: 
                    fa = {'x1':refals['x'][gdref[ind1]],'y1':refals['y'][gdref[ind1]],
                          'x2':als['x'][gdals[ind2]],'y2':als['y'][gdals[ind2]]}
                    initpar = np.copy(trans)
                    fpar,cov = curve_fit(trans_coo_dev,xdata,ydata,initpar)
                    trans = fpar
             
        # Match stars with X/Y coordinates 
        if (count < 1): 
            #MATCHSTARS,refals.x,refals.y,als.x,als.y,ind1,ind2,trans,count=count,/silent
            MATCHSTARS(refals[gdref].x,refals[gdref].y,als[gdals].x,als[gdals].y,ind1,ind2,trans,count=count)
             
        # No good matches, try srcmatch with "small" shifts 
        if (count < 1):
            ind1a,ind2a,dist = coords.xmatch(refals['x'][gdref],refals['y'][gdref],als['x'][gdals],als['y'][gdals],100,sphere=False)
            count = len(ind1)                                                            
            #SRCMATCH,refals[gdref].x,refals[gdref].y,als[gdals].x,als[gdals].y,100,ind1a,ind2a,count=count1 
            if count1 > 0: 
                xdiff1 = refals[gdref[ind1a]].x-als[gdals[ind2a]].x 
                ydiff1 = refals[gdref[ind1a]].y-als[gdals[ind2a]].y 
                xmed1 = np.median(xdiff1) 
                ymed1 = np.median(ydiff1) 
                # redo the search
                ind1,ind2,dist = coords.xmatch(refals['x'][gdref],refals['y'][gdref],als['x'][gdals]+xmed1,als['y'][gdals]+ymed1,20,sphere=False)
                count = len(ind1) 
                #SRCMATCH,refals[gdref].x,refals[gdref].y,als[gdals].x+xmed1,als[gdals].y+ymed1,20,ind1,ind2,count=count 
                if count == 0:
                    ind1,ind2,dist = coords.xmatch(refals['x'][gdref],refals['y'][gdref],als['x'][gdals]+xmed1,als['y'][gdals]+ymed1,100,sphere=False)
                    count = len(ind1) 
                    #SRCMATCH,refals[gdref].x,refals[gdref].y,als[gdals].x+xmed1,als[gdals].y+ymed1,100,ind1,ind2,count=count 
                xdiff = refals['x'][gdref[ind1]]-als['x'][gdals[ind2]]
                ydiff = refals['y'][gdref[ind1]]-als['y'][gdals[ind2]]
                xmed = np.median(xdiff) 
                ymed = np.median(ydiff) 
                trans = np.array([xmed, ymed, 1.0, 0.0, 0.0, 1.0])
             
        # No good match 
        if (count < 1): 
            print('NO MATCHES.  Using XSHIFT=YSHIFT=ROTATION=0')
            trans = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0])
             
        # Shift too large 
        if maxshift is not None:
            if np.max(np.abs(trans[0:1]))>maxshift: 
                print('SHIFTS TOO LARGE. ',str(trans[0:1]),' > ',str(maxshift),' Using XSHIFT=YSHIFT=ROTATION=0')
                trans = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0])
             
        # The output is: 
        # filename, xshift, yshift, 4 trans, FRAD (from als file), 0.0 
        fmt = '%2s%-30s%1s%10.2f%10.2f%10.5f%10.5f%10.5f%10.5f%10.3f%10.3f'
        data = "'",files[i],"'",*trans, frad, 0.0
        newline = fmt % data
        mchfinal += [newline]
             
        # Printing the transformation 
        if verbose: 
            print('%-20s%10.4f%10.4f%12.8f%12.8f%12.8f%12.8f' % (files[i],*trans))
         
    # Writing the final mchfile
    dln.writelines(mchbase+'.mch',mchfinal)

                    
    # Running DAOMASTER 
    #------------------
    daomaster(mchbase)
         
    # FAKE, copy back the original MCH file 
    if fake:
        dln.remove(mchbase+'.mch.daomaster',allow=True)
        shutil.copyfile(mchbase+'.mch',mchbase+'.mch.daomaster')
        dln.remove(mchbase+'.mch',allow=True)
        shutil.move(mchbase+'.mch.orig',mchbase+'.mch')
          
    # Print out the final transformations 
    if verbose:
        files,trans,magoff = io.readmch(mchbase+'.mch')
        nfiles = len(files) 
        print('')
        print('Final DAOMASTER Transformations:')
        for i in range(nfiles):
            print('%-20s%10.4f%10.4f%12.8f%12.8f%12.8f%12.8f' % (files[i],*list(trans[i,:])))
                  
    # Back to the original directory
    os.chdir(curdir)
    

def daomaster_tile(mchbase,info,tile,group,verbose=False):
    """
    Perform cross-matching of the various ALS catalogs
    and write out the .mch, tfr and .raw files that
    daomaster would create.
 
    Parameters
    ----------
    mchbase : str
       The base name of the main mch file.
    info : catalog
       Catalog of information about the fits and als files.
    tile : dict
       The tiling scheme information.
    groups : dict
       Grouping information.
    verbose : bool, optional
       Verbose output.  Default is True.
 
    Returns
    -------
    An MCH, TFR and RAW file with basename of the first file. 
 
    Example
    -------

    daomaster_tile(mchbase,info,tile,group)
  
    By D. Nidever   December 2006  
    Translated to Python by D. Nidever, April 2022

    """

    ###################### 
    # Create the RAW file 
    ###################### 
    # I can't use daomaster because it won't work for non-overlapping 
    # images, and it always makes the first image the reference. 

    nfiles = len(info)
    bases = np.char.array(info['base'])
     
    # Create the Schema
    dt = [('id',int),('x',float),('y',float)]
    for i in range(nfiles):
        dt += [('mag'+str(i+1),float),('err'+str(i+1),float)]
    dt += [('chi',float),('sharp',float),('nobs',int)]
    #schema = {id:0L,x:0.0,y:0.0} 
    #for i in range(nfiles): 
    #    schema=create_struct(schema,'MAG'+str(i+1,2),99.99,'ERR'+str(i+1,2),9.99) 
    #schema = create_struct(schema,'chi',99.99,'sharp',99.99,'nobs',0) 
    #rawtags = tag_names(schema) 
     
    # Number of sources in each als 
    nalsarr = np.zeros(nfiles,int)
    for i in range(nfiles): 
        nalsarr[i] = dln.numlines(bases[i]+'.als')-3
     
    # Initialize the RAW structure 
    #raw = REPLICATE(schema,100000L>nalsarr[0])
    raw = np.zeros(np.maximum(100000,nalsarr[0]),dtype=np.dtype(dt))
    raw = Table(raw)
    for i in range(nfiles):
        raw['mag'+str(i+1)] = 99.99
        raw['err'+str(i+1)] = 9.99
    raw['chi'] = 99.99
    raw['sharp'] = 99.99    
    rawcols = raw.colnames
    
    # Loop over the images 
    #---------------------
    tfr = np.zeros((len(raw),nfiles),float)  # tfr array
    count = 0
    for i in range(nfiles): 
        if verbose:
            print(str(i+1)+' '+bases[i])
         
        # Get the header
        if info['file'][i][-7:]=='fits.fz':
            fhead = io.readfile(info['file'][i],exten=1,header=True)
            # Fix the NAXIS1/2 in the header
            fhead['NAXIS1'] = fhead['ZNAXIS1']
            fhead['NAXIS2'] = fhead['ZNAXIS2']            
        else:
            fhead = io.readfile(info['file'][i],header=True)

        # Load the ALS file
        als,alshead = io.readals(bases[i]+'.als')
        nals = len(als)
        alsind = np.arange(nals)+1
        if i==0: 
            rawhead = alshead.copy()
         
        # Convert to tile coordinates
        fwcs = WCS(fhead)
        coo = fwcs.pixel_to_world(als['x']-1,als['y']-1)
        ra = coo.ra.deg
        dec = coo.dec.deg
        wcs = WCS(tile['head'])
        xref,yref = wcs.world_to_pixel(coo)
        xref += 1 # convert 0-indexed to 1-indexed 
        yref += 1 
        # Convert to tile group X/Y values 
        xref -= group['x0']
        yref -= group['y0']
         
        # Get mag/err columns
        magcol = 'mag'+str(i+1)
        errcol = 'err'+str(i+1)                  
        #magind, = np.where(rawtags == 'MAG'+str(i+1,2),nmagind) 
        #errind, = np.where(rawtags == 'ERR'+str(i+1,2),nerrind)
         
        # First image 
        if i==0: 
            raw['id'][0:nals] = np.arange(nals)+1 
            raw['x'][0:nals] = xref 
            raw['y'][0:nals] = yref 
            raw[magcol][0:nals] = als['mag']
            raw[errcol][0:nals] = als['err']
            raw['chi'][0:nals] = als['chi']
            raw['sharp'][0:nals] = als['sharp']
            raw['nobs'][0:nals] += 1
            # TFR 
            tfr[0:nals,i] = alsind 
            count += nals 
             
        # Second and later images, crossmatch 
        else: 
            # Cross-match
            ind1,ind2,dist = coords.xmatch(raw['x'][0:count],raw['y'][0:count],xref,yref,2,sphere=False)
            nmatch = len(ind1)
            #SRCMATCH,raw[0:count-1].x,raw[0:count-1].y,xref,yref,2,ind1,ind2,count=nmatch
            if verbose:
                print(str(nmatch)+' matches')

            # Some matches, add data to existing records for these sources 
            if nmatch > 0: 
                raw[magcol][ind1] = als['mag'][ind2]
                raw[errcol][ind1] = als['err'][ind2]
                raw['chi'][ind1] += als['chi'][ind2]     # cumulative sum 
                raw['sharp'][ind1] += als['sharp'][ind2]  # cumulative sum 
                raw['nobs'][ind1] += 1
                # TFR 
                tfr[ind1,i] = ind2 
                # Remove stars 
                if nmatch < nals:
                    als = np.delete(als,ind2)
                    xref = np.delete(xref,ind2)
                    yref = np.delete(yref,ind2)
                    alsind = np.delete(alsind,ind2)                  
                    #REMOVE,ind2,als,xref,yref,alsind 
                else: 
                    als = None
                nals = len(als)
             
            # Add leftover sources 
            if nals > 0: 
                # Add more elements to RAW and TFR 
                if count+nals > len(raw): 
                    if verbose:
                        print('Adding more elements')
                    raw = dln.add_elements(np.array(raw),np.maximum(100000,nals))
                    raw = Table(raw)
                    oldtfr = tfr 
                    del tfr
                    tfr = np.zeros((len(raw),nfiles),int)
                    tfr[0:len(oldtfr),:] = oldtfr
                    del oldtfr
                if verbose:
                    print('Adding '+str(nals)+' leftover sources')
                raw['id'][count:count+nals] = np.arange(nals)+1+count 
                raw['x'][count:count+nals] = xref 
                raw['y'][count:count+nals] = yref 
                raw[magcol][count:count+nals] = als['mag']
                raw[errcol][count:count+nals] = als['err']
                raw['chi'][count:count+nals] += als['chi']
                raw['sharp'][count:count+nals] += als['sharp']
                raw['nobs'][count:count+nals] += 1
                # TFR 
                tfr[count:count+nals,i] = alsind 
                count += nals
    # Trim extra elements 
    raw = raw[0:count] 
    tfr = tfr[0:count,:]
    # Calculate average chi/sharp 
    raw['chi'] /= np.maximum(raw['nobs'],1)
    raw['sharp'] /= np.maximum(raw['nobs'],1)
    nraw = len(raw)
                       
     
    # Write out the RAW file 
    #-----------------------
    with open(mchbase+'.raw','w') as f:
                  # Header                   
        f.write(rawhead[0]+'\n')
        f.write(rawhead[1]+'\n')
        f.write('\n')                  
        # Create MAG/ERR array 
        magarr = np.zeros((nraw,nfiles*2),float)
        for i in range(nfiles):
            magcol = 'mag'+str(i+1)
            errcol = 'err'+str(i+1)
            #magind , = np.where(rawtags == 'MAG'+str(i+1,2),nmagind) 
            #errind , = np.where(rawtags == 'ERR'+str(i+1,2),nerrind) 
            magarr[:,i*2] = raw[magcol]
            magarr[:,i*2+1] = raw[errcol]
        # Only 12 mag/err/chi/sharp columns per row
        nrows = int(np.ceil((nfiles*2+2)/12.))
        # Loop over raw elements 
        for i in range(nraw): 
            # The floating point numbers, MAG, ERR, CHI, SHARP 
            arr = np.copy(magarr[i,:])
            arr = np.append(arr,raw['chi'][i])
            arr = np.append(arr,raw['sharp'][i])
            narr = len(arr) 
            # Loop over rows for this object 
            for j in range(nrows): 
                if narr>12:
                    thisarr = arr[0:12]
                    arr = arr[12:] 
                    narr = len(arr) 
                else: 
                    thisarr = arr 
                    arr = None
                    narr = 0 
                if j==0:
                    fmt = '%7d%9.3f%9.3f'+(len(thisarr)*'%9.4f')+'\n'
                    f.write(fmt % (raw['id'][i],raw['x'][i],raw['y'][i],*thisarr))
                else: 
                   # 27 leading spaces
                   fmt = '%25s'+(len(thisarr)*'%9.4f')+'\n'
                   f.write(fmt % ('',*thisarr))
     
     
    # Write out TFR file 
    #-------------------
    with open(mchbase+'.tfr','w') as f:
        for i in range(nfiles):
            fmt = '%1s%-30s%9.4f%9.4f\n'
            f.write(fmt % ('',bases[i]+'.als',99.9999,9.9999))
            #printf,unit,'',bases[i]+'.als',99.9999,9.9999,format='(A1,A-30,F9.4,F9.4)' 
        f.write(' ==============================\n')
        #format = '(I7,F9.2,F9.2,'+str(nfiles,2)+'I7)'
        fmt = '%7d%9.2f%9.2f'+(nfiles*'%7d')+'\n'
        for i in range(nraw): 
            f.write(fmt % (raw['id'][i],raw['x'][i],raw['y'][i],*list(tfr[i,:])))


def daomatch_tile(files,tile,group,mchbase=None,verbose=False,logfile=None,
                  maxshift=5000,fake=False):
    """
    This is very similar to the DAOMATCH.PRO program that 
    matches stars and finds the transformation but it does 
    it using the "tiling" coordinate system. 
 
    Parameters
    ----------
    files : list
       Array of ALS files,  The first file will be used 
         as the reference.  If /fake set then the first ALS file 
         in "files" should already have an associated MCH file. 
    tile : dict
       The tiling scheme information.
    groups : dict
       Grouping information.
    mchbase : str
       The base name of the main mch file.
    maxshift : int, optional
       Constraints on the initial X/Y shifts. 
    usewcs : bool, optional
       Use the WCS for initial matching.  This is the default.
    fake : bool, optional
       Run for artificial stars.  The MCH file should be input 
         and daomaster run to create raw/tfr files. 
    verbose : bool, optional
       Verbose output.  Default is False.
 
    Returns
    -------
    An MCH, TFR and RAW file with basename of the first file. 
 
    Example
    -------

    daomatch_tile(['obj1034_1.als','obj1035_1.als','obj1036_1.als'])
 
    Add options to freeze the scale and/or rotation. 
 
    By D. Nidever   December 2006  
    Translated to Python by D. Nidever, April 2022
    """

    t0 = time.time() 

    nfiles = dln.size(files)

    # Logfile 
    if logfile is not None:
        logf = logfile 
    else: 
        logf = -1 
     
    # Only one file, can't match 
    if nfiles == 1: 
        raise ValueError('ONLY ONE FILE INPUT.  No matching, simply creating .mch and .raw file')

    # Current directory
    curdir = os.getcwd()
    fdir = os.path.abspath(os.path.dirname(files[0]))
    os.chdir(fdir)

    bases = [os.path.splitext(os.path.basename(f))[0] for f in files]
    bases = np.char.array(bases)

    # FAKE, running for artificial star tests 
    if fake: 
        # Check that MCH file exists 
        if os.path.exists(mchbase+'.mch') == 0: 
            error = '/fake set but '+mchbase+'.mch NOT FOUND' 
            printlog,logf,error 
            return 
        # Skip the MCH creation process 
        goto,rundaomaster 
     
     
    # Gather information on all of the files 
    print('Gathering file information')
    fitsfiles = bases+'.fits' 
    bdfits, = np.where(dln.exists(fitsfiles)==False)
    if len(bdfits)>0: 
        for b in bdfits:
            fitsfiles[b] += '.fz'
    info = io.fileinfo(fitsfiles)
    ntrans = 6
    info['catfile'] = 100*' '
    info['resamptrans'] = np.zeros((nfiles,ntrans),float)
    info['resamptransrms'] = 0.0
    info['catfile'] = np.char.array(bases)+'.als' 
     
    # Creating MCH file in the tile coordinate system
    mchfinal = []
    for i in range(nfiles): 
        # Get the header
        if info['file'][i][-7:]=='fits.fz':
            fhead = io.readfile(info['file'][i],exten=1,header=True)
            # Fix the NAXIS1/2 in the header
            fhead['NAXIS1'] = fhead['ZNAXIS1']
            fhead['NAXIS2'] = fhead['ZNAXIS2']  
        else:
            fhead = io.readfile(info['file'][i],header=True)            
         
        # Convert X/Y of this system into the combined reference frame 
        #  The pixel values are 1-indexed like DAOPHOT uses. 
        #  Use a 2D grid of points in the image and the WCS to get 
        #  the transformation. 
        ngridbin = 50 
        nxgrid = info['nx'][i] // ngridbin 
        nygrid = info['ny'][i] // ngridbin 
        xgrid = (np.arange(nxgrid)*ngridbin+1).reshape(-1,1) + np.zeros(nygrid,int).reshape(1,-1)
        ygrid = np.zeros(nxgrid,int).reshape(-1,1) + (np.arange(nygrid)*ngridbin+1).reshape(1,-1)

        fwcs = WCS(fhead)
        coogrid = fwcs.pixel_to_world(xgrid-1,ygrid-1)
        ragrid = coogrid.ra.deg
        decgrid = coogrid.dec.deg
        wcs = WCS(tile['head'])
        refxgrid,refygrid = wcs.world_to_pixel(coogrid)
        refxgrid += 1  # convert 0-indexed to 1-indexed 
        refygrid += 1 
        # Convert to tile X/Y values 
        refxgrid -= group['x0']
        refygrid -= group['y0']
         
        # Now fit the transformation 
        xdiff = refxgrid-xgrid 
        ydiff = refygrid-ygrid 
        xmed = np.median(xdiff) 
        ymed = np.median(ydiff) 
        # Fit rotation with linear fits if enough points
        slp1 = dln.mediqrslope(ygrid,xdiff)  # fit rotation term
        slp1rms = dln.mad(xdiff-ygrid*slp1)
        slp2 = dln.mediqrslope(xgrid,ydiff)
        slp2rms = dln.mad(ydiff-xgrid*slp2)
        theta = dln.wtmean(np.array([-slp1,slp2]),np.array([slp1rms,slp2rms]))
        #coef1 = robust_poly_fitq(ygrid,xdiff,1)# fit rotation term 
        #coef1b = dln_poly_fit(ygrid,xdiff,1,measure_errors=xdiff*0+0.1,sigma=coef1err,/bootstrap) 
        #coef2 = robust_poly_fitq(xgrid,ydiff,1)# fit rotation term 
        #coef2b = dln_poly_fit(xgrid,ydiff,1,measure_errors=ydiff*0+0.1,sigma=coef2err,/bootstrap) 
        ##theta = mean([-coef1[1],coef2[1]]) 
        #WMEANERR,[-coef1[1],coef2[1]],[coef1err[1],coef2err[1]],theta,thetaerr 
         
        # [xoff, yoff, cos(th), sin(th), -sin(th), cos(th)] 
        trans = np.array([xmed, ymed, 1.0-theta**2, theta, -theta, 1.0-theta**2])
        # Adjust Xoff, Yoff with this transformation 
        xout,yout = utils.trans_coo([xgrid,ygrid],*trans) 
        trans[0] += np.median(refxgrid-xout)
        trans[1] += np.median(refygrid-yout)

        # trans rotates xgrid/ygrid to refxgrid/ygrid

        # Fit full six parameters if there are enough stars 
        #fa = {x1:(refxgrid)(*),y1:(refygrid)(*),x2:(xgrid)(*),y2:(ygrid)(*)} 
        #initpar = trans 
        #fpar = MPFIT('trans_coo_dev',initpar,functargs=fa, perror=perror, niter=iter, status=status,
        # bestnorm=chisq,:f=dof, autoderivative=1, /quiet) 
        # trans = fpar 
        xdata = [[refxgrid,refygrid],[xgrid,ygrid]]
        null = np.zeros(refxgrid.size,float)
        fpar,cov = curve_fit(utils.trans_coo_dev,xdata,null,p0=trans)
        trans = fpar
        transerr = np.sqrt(np.diag(cov))
        diff = utils.trans_coo_dev(xdata,*fpar)
        rms = np.sqrt(np.mean(diff**2.)) 
        info['resamptrans'][i] = trans 
        info['resamptransrms'][i] = rms 
         
        # The output is: 
        # filename, xshift, yshift, 4 trans, mag offset, magoff sigma 
        #format = '(A2,A-30,A1,2A10,4A12,F9.3,F8.4)' 
        # In daomaster.f the translations are 10 digits with at most 4 
        # decimal places (with a leading space), the transformation 
        # coefficients are 12 digits with at most 9 decimal places. 
        # Need a leading space to separate the numbers. 
        strans = ['%30.4f' % trans[0], '%30.4f' % trans[1], '%30.9f' % trans[2], '%30.9f' % trans[3],
                  '%30.9f' % trans[4], '%30.9f' % trans[5]]
        strans = [' '+s.strip() for s in strans]
        fmt = '%2s%-30s%1s%10.10s%10.10s%12.12s%12.12s%12.12s%12.12s%9.3f%8.4s'
        data = "'",info['catfile'][i],"'",*strans, 0.0, rms,
        newline = fmt % data
        mchfinal += [newline]
         
        # Printing the transformation 
        if verbose:
            print('%-22s%10.10s%10.10s%12.12s%12.12s%12.12s%12.12s%9.3f%8.4f' % (info['catfile'][i],*strans,0.0,rms))

    # Write to the new MCH file 
    mchbase = bases[0] 
    mchfile = mchbase+'.mch'
    dln.writelines(mchfile,mchfinal)

    # Create the raw and tfr files
    daomaster_tile(mchbase,info,tile,group,verbose=verbose)
     
    # Back to the original directory
    os.chdir(curdir)
     
