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
