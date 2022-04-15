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


def pixscale(filename,head=None):
    """
    Get the pixel scale for an image. 
 
    Parameters
    ----------
    file    FITS filename 
    =head   The image header for which to determine the pixel scale. 
    /stp    Stop at the end of the program. 
 
    Returns
    -------
    scale   The pixel scale of the image in arcsec/pix. 
    =error  The error if one occurred. 
 
    Example
    -------

    scale = pixscale('ccd1001.fits')
 
    BY D. Nidever   February 2008 
    Translated to Python by D. Nidever,  April 2022
    """

    scale = None  # bad until proven good 

     
    # No header input, read from fits file 
    fpack = 0 
    if head is None:
        # Check that the file exists
        if os.path.exists(filename)==False:
            raise ValueError(filename+' NOT FOUND')

        # Open the file
        hdu = fits.open(filename)
        
        # Fpack or regular fits
        if filename[-7:]=='fits.fz':
            fpack = 1 
            exten = 1 
        else: 
            fpack = 0 
            exten = 0 
         
        # Read the header
        if head is None:
            head = readfile(filename,exten=exten,header=True) 
         
        # Fix NAXIS1/2 in header 
        if fpack == 1:
            head['NAXIS1'] = head['ZNAXIS1']
            head['NAXIS2'] = head['ZNAXIS2']            
     
    # Does the image have a SCALE parameter
    if head.get('scale') is not None:
        scale = head['scale']
    # Does the image have a PIXSCALE parameter 
    if scale is None:
        pixscale = head.get('pixscale')
        if pixscale is not None:
            scale = pixscale
    # Does the image have a PIXSCALE1 parameter 
    if scale is None: 
        pixscale1 = head.get('pixscale1')
        if npixscale1 is not None:
            scale = pixscale1
     
    # Try the WCS 
    if scale is None:
        try:
            wcs = WCS(head)
             
            # Get the coordinates for two positions 
            #  separated by 1 pixel 
            #head_xyad,head,0.0,0.0,ra1,dec1,/degree 
            #head_xyad,head,1.0,0.0,ra2,dec2,/degree 
            #dist = sphdist(ra1,dec1,ra2,dec2,/deg)*3600. 
            #scale = dist
            ra1,dec1 = wcs.pixel_to_world(0,0,0)
            ra2,dec2 = wcs.pixel_to_world(1,0,0)            
            dist = dln.sphdist(ra1,dec1,ra1,dec1)*3600
            scale = dist
            
            if scale == 0.0: 
                scale = None 
        except:
            raise ValueError('No WCS in header')
                
    # Couldn't determine the pixel scale 
    if scale == None:
        error = 'WARNING! COULD NOT DETERMINE THE PIXEL SCALE' 
        print(error)

    return scale


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
