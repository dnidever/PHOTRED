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
from astropy.time import Time
from astropy.table import Table
from astropy.wcs import WCS
from dlnpyutils import utils as dln
import struct
from itertools import zip_longest
from itertools import accumulate
from io import StringIO

def isfloat(s):
    """ Returns True is string is a number. """
    try:
        float(s)
        return True
    except ValueError:
        return False

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

def file_wait(filename,wait=5,timeout=600,silent=False):
    """
    Wait until a file exists. 
 
    Parameters
    ----------
    filename : str
       Filename to check. 
    wait : int, optional
       Wait time between checks.  Default is 5 sec. 
    timeout : int, optional
       Stop trying and throw an error after this time.  Default 
         is 600 sec. 
    silent : boolean, optional
       Don't print anything to the screen 
 
    Returns
    -------
    None 
 
    Example
    -------

    file_wait('file.txt')
 
    By D. Nidever  April 2020 
    Translated to Python by D. Nidever,  April 2022
    """

    # While loop 
    t0 = time.time() 
    while (os.path.exists(filename) == False):
        if silent==False:
            print(filename+' NOT FOUND.  Waiting '+str(wait)+' sec.')
        time.sleep(wait)
        if time.time()-t0 > timeout: 
            raise ValueError('Timeout ('+str(timeout)+' sec) reached on '+filename)


def fitsext(files,isfpack=False,basename=False,full=False):
    """
    Gets the fits extension (.fits or .fits.fz) and optionaly 
    the base name. 
 
    Parameters
    ----------
    files : str or list
       Scalar or array of file names. 
    isfpack : boolean
       Returns a boolean True if the file is a fpack compressed 
         FITS file and False if not.  Default is False.
    basename : boolean, optional
       Returns the basename only.  Default is False.
    full : boolean, optional
       Returns an array of basename and extension.  If a 
         single file is input then the output will a 2-element 
         array, otherwise it is [Nfiles,2].  Default is False.
 
    Returns
    -------
    output : str
      The output of the program which depends on the inputs. 
        By default this is the FITS extension (.fits or .fits.fz) 
        and '' if none of those. 
          isfpack   Returns a boolean True or False if the file(s)
                       are fits.fz extensions. 
          basename  Returns the file basename (without the 
                       extension or directory). 
          full      Returns the file basename and extension in 
                       2-element or [Nfiles,2] array depending 
                      on how many files were input. 
 
    Example
    -------

    ext = fitsext('F2-04958679_01.fits.fz') 
 
    By D. Nidever  August 2017 
    Translated to Python by D. Nidever,  April 2022
    """
 
    # Cases 
    # 1 - extension only 
    # 2 - isfpack, boolean 1/0 
    # 3 - basename only 
    # 4 - full, basename and extension 
    usecase = 1 # extension by default
    if isfpack:
        usecase = 2
    if basename:
        usecase = 3
    if full:
        usecase = 4 

    nfiles = np.array(files).size
    if nfiles==1 and (type(files)==str or type(files)==np.str or type(files)==np.str_):
        files = [files]
    
    # Initializing the output array 
    out = None
    # Extension only     
    if usecase==1:
        out = np.zeros(nfiles,(np.str,100))
    # isfpack, boolean 1/0 
    elif usecase==2:
        out = np.zeros(nfiles,bool)
    # basename only 
    elif usecase==3:
        out = np.zeros(nfiles,(np.str,100))
    # full, basename and extension
    elif usecase==4:
        if nfiles == 1:
            out = np.zeros(2,(np.str,100))
        else:
            out = np.zeros((nfiles,2),(np.str,100))
    else: # Not supported 
        pass
         
    # Loop through the files 
    for i in range(nfiles): 
        file1 = os.path.basename(files[i].strip())
             
        # Get the extension 
        ext = ''
        if file1[-7:] == 'fits.fz':
            ext = '.fits.fz'
        elif file1[-4:]=='fits':
            ext = '.fits' 
             
        # Extension only 
        if usecase==1:
            out[i] = ext 
        # isfpack, boolean
        elif usecase==2:
            if ext == '.fits.fz':
                out[i] = True
        # basename only
        elif usecase==3:
            out[i] = file1[0:-len(ext)]
        # full, basename and extension 
        elif usecase==4:
            if nfiles == 1: 
                out[0] = file1[0:-len(ext)]
                out[1] = ext 
            else: 
                out[i,0] = file1[0:-len(ext)]
                out[i,1] = ext 
        else:  # Not supported 
            pass

    if nfiles==1 and out.ndim==1:
        out = out[0]
     
    return out 

def date2jd(dateobs,mjd=False):
    """ Converte DATE-OBS to JD or MJD."""
    t = Time(dateobs, format='fits')
    if mjd:
        return t.mjd
    else:
        return t.jd

def trans_coo(xdata,*par):
    """ Apply the transformation to X/Y"""

    xin = xdata[0]
    yin = xdata[1]

    A = par[0]
    B = par[1]
    C = par[2]
    D = par[3]
    E = par[4]
    F = par[5]
    
    # from ccdpck.txt
    #              x(1) = A + C*x(n) + E*y(n)
    #              y(1) = B + D*x(n) + F*y(n)

    # Apply transformation
    xout = A + C*xin + E*yin
    yout = B + D*xin + F*yin
    
    return xout,yout


def trans_coo_dev(xdata,*par):

    # Rotate coordinates(2) to coordinate system 1
    # and return deviates
    x1,y1 = xdata[0]
    x2,y2 = xdata[1]

    newx2,newy2 = trans_coo([x2,y2],*par)

    diff = np.sqrt( (x1-newx2)**2 + (y1-newy2)**2 )

    # Do robust outlier rejection
    std = dln.mad(diff)
    med = np.median(diff)
    bd = (diff > (med+3.0*std))
    if np.sum(bd)>0:
        diff[bd] = 0.0

    return diff.flatten()

def trans_coo_outlier(xdata,*par):

    # Rotate coordinates(2) to coordinate system 1
    # and return deviates
    x1,y1 = xdata[0]
    x2,y2 = xdata[1]

    newx2,newy2 = trans_coo([x2,y2],*par)

    diff = np.sqrt( (x1-newx2)**2 + (y1-newy2)**2 )

    # Do robust outlier rejection
    std = dln.mad(diff)
    med = np.median(diff)
    bd = (diff > (med+3.0*std))
    if np.sum(bd)>0:
        newx2[bd] = x1[bd]
        newy2[bd] = y1[bd]

    return np.append(newx2,newy2)
 
def validtile(tile):
    """
    This double-checks if the TILE structure is valid. 
 
    Parameters
    ----------
    tile : dict
      The tile structure 
 
    Returns
    -------
    check : boolean
      1-if the tile is good and 0-if there is 
        a problem with it or it doesn't exist. 
 
    Example
    -------

    check = validtile(tile) 
 
    By D. Nidever  Oct 2016 
    Translated to Python by D. Nidever, April 2022
    """

    # Must be a dictionary
    if type(tile) != dict:
        return 0
     
    # Must have TYPE column 
    if 'type' not in tile.keys():
        return 0
     
    # Do we have enough information for each type 
    if tile['type'].lower()=='orig':
        return 1 
    elif tile['type'].lower()=='wcs':
        # wcs exists 
        if 'wcs' in tile.keys():
            # wcs must be a WCS object
            if type(tile['wcs']) != WCS:
                return 0         
        # NO wcs, check other needed values 
        else: 
            # Need NAXIS, CRVAL, CRPIX, CTYPE, CDELT 
            for k in ['naxis','crval','crpix','ctype','cdelt']:
                if k not in tile.keys():
                    return 0
    # Need XRANGE, YRANGE 
    if 'xrange' not in tile.keys():
        return 0
    if 'yrange' not in tile.keys():
        return 0 
 
    return 1 

def photred_mkopt(inpfiles,hilimit=6.4e4,va=2,fitradius_fwhm=None,
                  inp_fwhm=None,verbose=True):
    """
    This makes opt files for FITS files to be used with 
    DAOPHOT and ALLSTAR in the PHOTRED pipeline 
 
    Parameters
    ----------
    inpfiles : str
      Input files. Three formats can be used (1) Name of file 
        with a list of input filenames.  Must start with an '@'; 
        (2) A name with wildcard characters, such as '*'; 
        (3) An array of filenames. 
    hilimit : float, optional
      The saturation upper limit, 64,000 by default. 
    va : int, optional
      The spatial variable PSF setting to use.  Default is 2.
    fitradius_fwhm : float, optional
      The value to use for the fitting radius (FI), in 
         units of the FWHM. 
    inp_fwhm : float, optional
      Use this FWHM. 
    verbose : boolean, optional
      Output information about what is happening. 
 
    Returns
    -------
    Makes .opt and .als.opt files for each FITS file, in the same 
    directory that the FITS file is in. 
 
    fwhm : float
      The image FWHM.

    Example
    -------

    fwhm = mkopt('mkopt.lst')
 
    Very similar to Tony Sohn's mkopt.f fortran program 
    but estimate FWHM automatically with IMFWHM.PRO 
    and gets RDNOISE and GAIN directly from the image 
    headers. 
    
    By D.Nidever  May 2008   basically a copy of MKOPT.PRO 
                             which was copied from Tony's mkopt 
                             fortran program 
    Translated to Python by D. Nidever, May 2022
    """ 
 
     
    # Loading input 
    files = dln.loadinput(inpfiles)
    count = len(files)
     
    # Not enough inputs 
    if nfiles == 0: 
        raise ValueError('No files') 
     
    # More than one name input 
    if nfiles > 1: 
        fwhm = fltarr(nfiles) 
        for i in range(nfiles): 
            fwhm1 = mkopt(files[i],hilimit=hilimit,va=va,fitradius_fwhm=fitradius_fwhm,verbose=verbose)
            fwhm[i] = fwhm1
            if verbose):
                print('')
        return fwhm
     
    filename = str(files[0]).strip()
     
    # Default settings
    va = np.maximum(va,2)
    if fitradius_fwhm is not None:
        fitradius_fwhm = 1.0
    fitradius_fwhm = np.minimum(fitradius_fwhm,0.0)
         
         
    # Processing ONE file 
    #-------------------- 
    if os.path.exits(filename)==False:
        raise ValueError(filename+' NOT FOUND')

    if verbose:
        print('Running MKOPT on ',filename)
         
    # Get the base
    base = fitsext(filename,basename=True)
    #if strmid(str(filename,2),6,7,/reverse_offset) == 'fits.fz': 
    #    fpack = 1 
    #    base = str(os.path.basename(filename,'.fits.fz'),2) 
    #else: 
    #    fpack = 0 
    #    base = str(os.path.basename(filename,'.fits'),2) 
    fdir = os.path.dirname(filename) 
         
             
    # Get the FITS header 
    if fpack == 1:
        head = io.readfile(filename,exten=1,header=True)
    else: 
        head = io.readfile(filename,header=True)
         
    # We need GAIN, READNOISE, FWHM, and Hi-limit 
    #-------------------------------------------- 
         
    # Getting GAIN
    gain = io.getgain(filename)
    # Getting READNOISE 
    rdnoise = io.getrdnoise(filename)
         
    # Run IMFWHM to get the FWHM
    if inp_fwhm is None:
        fwhm,gtab = imfwhm(filename,im=im)
        # somtimes imfwhm has problems if the saturation level in the 
        # header is too low, run without header 
        if (fwhm > 20 or len(gtab) < 10) and len(im) == 0: 
            fwhm1 = fwhm 
            print('Running IMFWHM again.  FWHM too high or number of sources too small')
            im = fits.getdata(filename)
            fwhm,tab = imfwhm('',im=im,verbose=False)
            # Still bad, using original one if possible 
            if fwhm > 90.0 and fwhm1 < 20: 
                print('New FWHM is bad but original FWHM was acceptable.  Using it')
                fwhm = fwhm1 
             
        if fwhm > 90.0: 
            print('Error with FWHM')
            return 
             
    else:
        fwhm = inp_fwhm 
             
    # Load the image 
    im,head = io.readfile(filename,header=True)
             
    # Getting saturation limit from the header 
    lolimit = 10000.0  # just in case 
    saturate = head.get('SATURATE')
    #if nsaturate eq 0 then saturate=(max(im) < hilimit)  ; if not found 
    if saturate is None:
        saturate = dln.limit(np.max(im)-1000,lolimit, hilimit)
             
    # Minimum of all saturation levels 
    #hi = lolimit > ( (saturate - 4000.0) < hilimit ) 
    #hi = lolimit > ( (saturate - 1000.0) < hilimit ) 
    # Don't constrain the saturation value that is input 
    hi = saturate 
             
    # Verbose output 
    if verbose:
        print('gain = ',str(gain)) 
        print('rdnoise = ',str(rdnoise)) 
        print('fwhm = ',str(fwhm)) 
        print('saturation = ',str(hi))
             
             
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    #% MAKING THE OPT FILES 
             
             
    # (1) DAOPHOT parameters 
    # 
    # LO    : Low good datum (7. works fine on most imags) 
    # TH    : Threshold (3.5 works fine) 
    # LS,HS : Low and high sharpness (default : 0.2 - 1.0) 
    # LR,HR : Low roundness and high roundness (default : -1.0 - 1.0) 
    # WA    : Watch progress 
    # VA    : Variable PSF 
    # AN    : Analytic model PSF 
    # EX    : Extra PSF cleaning passes 
    # PE    : Percent error 
    # PR    : Profile error 
             
    # (2) ALLSTAR parameters 
    # 
    # CR    : Clipping range (leave it) 
    # CE    : Clipping exponent (leave it) 
    # MA    : Maximum group size 
    # RED   : Redetermine centroid (0 = no, 1 = yes) 
             
    # Frame-specific parameters. 
    # 
    # GA    : gain (e/ADU) 
    # RD    : readout noise (e) 
    # RE    : readout noise (ADU) 
    # FW    : FWHM 
    # HI    : hi good datum in ADU - saturation level 
    # FI    : fitting radius 
    # PS    : PSF radius 
    # IS,OS : inner and outer sky annalus 
             
    LO =  7.0 
    TH =  3.5 
    LS =  0.2 
    HS =  1.0 
    LR = -1.0 
    HR =  1.0 
    WA = -2 
    # VA  defined above 
    AN = -6  # It will try all PSF models (#1-6) and use the one with the lowest chi value 
    EX =  5  # extra PSF passes 
    PE =  0.75 
    PR =  5.00 
    CR =  2.5 
    CE =  6.0 
    MA = 50. 
    RED = 1.0 
    WA2 = 0.0 
             
    # Frame specific parameters 
    GA = gain 
    RD = rdnoise 
    FW = fwhm 
    #HI = hi 
             
    # Calculating some things 
    FW = np.minimum(FW, 20)  # daophot won't accept anything higher than 20 
    RE = np.maximum(RD/GA, 0.01) 
    FI = np.minimum(fitradius_fwhm*FW, 51)  # daophot won't accept anything higher than 51 
    PS = np.minimum((4.0*FW), 51)      # daophot won't accept anything higher than 51 
    IS = np.minimum((FI - 1.0), 35)    # daophot won't accept anything higher than 35 
    OS = np.minimum((PS + 1.0), 100)   # daophot won't accept anything higher than 100 
             
    # Writing the DAOPHOT parameter 
    #------------------------------ 
    # 
    # RE    : readout noise (ADU) 
    # GA    : gain (e/ADU) 
    # LO    : Low good datum (7. works fine on most imags) 
    # HI    : hi good datum in ADU - saturation level 
    # FW    : FWHM 
    # TH    : Threshold (3.5 works fine) 
    # LS,HS : Low and high sharpness (default : 0.2 - 1.0) 
    # LR,HR : Low roundness and high roundness (default : -1.0 - 1.0) 
    # WA    : Watch progress 
    # FI    : fitting radius 
    # PS    : PSF radius 
    # VA    : Variable PSF 
    # AN    : Analytic model PSF 
    # EX    : Extra PSF cleaning passes 
    # PE    : Percent error 
    # PR    : Profile error 
             
    outarr = [RE,GA,LO,HI,FW,TH,LS,HS,LR,HR,WA,FI,PS,VA,AN,EX,PE,PR] 
    anotarr = ['RE','GA','LO','HI','FW','TH','LS','HS','LR','HR','WA','FI','PS','VA','AN','EX','PE','PR'] 
    anotarr = np.char.array(anotarr)+' = ' 
    nanot = len(anotarr) 

    with open(fdir+'/'+base+'.opt','w') as f:
        for j in range(nanot):
            form = '%5s%8.2f'
            if anotarr[j] == 'HI = ':
                form = '%5s,%8d'
            f.write(form % (anotarr[j],outarr[j]))
             
    # Writing the ALLSTAR parameter file 
    #----------------------------------- 
             
    # FI    : fitting radius 
    # IS    :  ?? 
    # OS    :  ?? 
    # RED   : Redetermine centroid (0 = no, 1 = yes) 
    # WA2   : Watch progress 
    # PE    : Percent error 
    # PR    : Profile error 
    # CR    : Clipping range (leave it) 
    # CE    : Clipping exponent (leave it) 
    # MA    : Maximum group size 
             
    outarr2 = [FI,IS,OS,RED,WA2,PE,PR,CR,CE,MA] 
    anotarr2 = ['FI','IS','OS','RE','WA','PE','PR','CR','CE','MA'] 
    anotarr2 = np.char.array(anotarr2)+' = ' 
    nanot2 = len(anotarr2) 

    with open(fdir+'/'+base+'.als.opt','w') as f:
        for j in range(nanot2):
            form = '%5s%8.2f'
            f.write(form % (anotarr2[j],outarr2[j]))
             
    # Verbose output 
    if verbose:
        print('Created ',fdir+'/'+base+'.opt')
        print('Created ',fdir+'/'+base+'.als.opt') 
             
    return fwhm

