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
    if nfiles==1 and type(files)==str:
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
     
    return out 


def trans_coo(xin,yin,par):
    """ Apply the transformation to X/Y"""

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


def trans_coo_dev(par,x1=x1,y1=y1,x2=x2,y2=y2):

    # Rotate coordinates(2) to coordinate system 1
    # and return deviates

    newx2,newy2 = trans_coo(x2,y2,par)

    diff = np.sqrt( (x1-newx2)**2 + (y1-newy2)**2 )

    # Do robust outlier rejection
    std = dln.mad(diff)
    med = dln.median(diff)
    bd, = np.where(diff > (med+3.0*std))
    if len(bd)>0:
        diff[bd] = 0.0

    return diff




