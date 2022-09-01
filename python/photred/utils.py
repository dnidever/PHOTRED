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


def meanclip(image,clipsig=3,maxiter=5,converge_num=0.02,verbose=False):
    """
    Computes an iteratively sigma-clipped mean on a data set,

    Parameters
    ----------
    image : numpy array
      Input data, any numeric array.
    CLIPSIG:  Number of sigma at which to clip.  Default=3
    MAXITER:  Ceiling on number of clipping iterations.  Default=5
    CONVERGE_NUM:  If the proportion of rejected pixels is less
        than this fraction, the iterations stop.  Default=0.02, i.e.,
        iteration stops if fewer than 2% of pixels excluded.
    VERBOSE:  Set this flag to get messages.

    Returns
    -------
    Mean:     N-sigma clipped mean.
    Sigma:    Standard deviation of remaining pixels.
    SUBS:     Subscript array for pixels finally used.

    """

    subs, = np.where(np.isfinite(image.ravel()))
    ct = len(subs)
    niter = 0
    endflag = False
    while (endflag==False):
        skpix = image.ravel()[subs]
        niter += 1
        lastct = ct
        medval = np.median(skpix)
        sig = np.sqrt(np.var(skpix))
        wsm = (np.abs(skpix-medval) < clipsig*sig)
        ct = np.sum(wsm)
        if ct > 0:
            subs = subs[wsm]         
        if (np.float(np.abs(ct-lastct))/lastct <= converge_num) or (niter > maxiter) or (ct == 0):
            endflag = True

    mean = np.mean(image.ravel()[subs])
    sigma = np.sqrt(np.var(image.ravel()[subs]))

    return mean, sigma, subs



def mmm(sky_vector, highbad=None,debug=False,readnoise=None,
        integer=False,maxiter=50,silent=False,minsky=20):
    """
    Estimate the sky background in a stellar contaminated field. 

    MMM assumes that contaminated sky pixel values overwhelmingly display 
    POSITIVE departures from the true value.  Adapted from DAOPHOT 
    routine of the same name. 
     
    CALLING SEQUENCE: 
           MMM, sky, [ skymod, sigma, skew, HIGHBAD = , READNOISE=, /DEBUG, 
                      MINSKY=, NSKY=, /INTEGER,/SILENT] 
     
    INPUTS: 
           SKY - Array or Vector containing sky values.  This version of 
                   MMM does not require SKY to be sorted beforehand.  SKY 
                   is unaltered by this program. 
     
    OPTIONAL OUTPUTS: 
           skymod - Scalar giving estimated mode of the sky values (float) 
           SIGMA -  Scalar giving standard deviation of the peak in the sky 
                   histogram.  If for some reason it is impossible to derive 
                   skymod, then SIGMA = -1.0 
           SKEW -   Scalar giving skewness of the peak in the sky histogram 
     
                   If no output variables are supplied or if /DEBUG is set 
                   then the values of skymod, SIGMA and SKEW will be printed. 
     
    OPTIONAL KEYWORD INPUTS: 
           HIGHBAD - scalar value of the (lowest) "bad" pixel level (e.g. cosmic 
                    rays or saturated pixels) If not supplied, then there is 
                    assumed to be no high bad pixels. 
           MINSKY - Integer giving mininum number of sky values to be used.   MMM 
                    will return an error if fewer sky elements are supplied. 
                    Default = 20. 
           MAXITER - integer giving maximum number of iterations allowed,default=50 
           READNOISE - Scalar giving the read noise (or minimum noise for any 
                    pixel).     Normally, MMM determines the (robust) median by 
                    averaging the central 20% of the sky values.     In some cases 
                    where the noise is low, and pixel values are quantized a 
                    larger fraction may be needed.    By supplying the optional 
                    read noise parameter, MMM is better able to adjust the 
                    fraction of pixels used to determine the median. 
           /INTEGER - Set this keyword if the  input SKY vector only contains 
                    discrete integer values.    This keyword is only needed if the 
                    SKY vector is of type float or double precision, but contains 
                    only discrete integer values.     (Prior to July 2004, the 
                    equivalent of /INTEGER was set for all data types) 
           /DEBUG - If this keyword is set and non-zero, then additional 
                   information is displayed at the terminal. 
           /SILENT - If set, then error messages will be suppressed when MMM 
                    cannot compute a background.    Sigma will still be set to -1 
     OPTIONAL OUTPUT KEYWORD: 
          NSKY - Integer scalar giving the number of pixels actually used for the 
                 sky computation (after outliers have been removed). 
     NOTES: 
           (1) Program assumes that low "bad" pixels (e.g. bad CCD columns) have 
           already been deleted from the SKY vector. 
           (2) MMM was updated in June 2004 to better match more recent versions 
           of DAOPHOT. 
           (3) Does not work well in the limit of low Poisson integer counts 
           (4) MMM may fail for strongly skewed distributions. 
     METHOD: 
           The algorithm used by MMM consists of roughly two parts: 
           (1) The average and sigma of the sky pixels is computed.   These values 
           are used to eliminate outliers, i.e. values with a low probability 
           given a Gaussian with specified average and sigma.   The average 
           and sigma are then recomputed and the process repeated up to 20 
           iterations: 
           (2) The amount of contamination by stars is estimated by comparing the 
           mean and median of the remaining sky pixels.   If the mean is larger 
           than the median then the true sky value is estimated by 
           3*median - 2*mean 
     
     REVISION HISTORY: 
           Adapted to IDL from 1986 version of DAOPHOT in STSDAS, 
           W. Landsman, STX Feb 1987 
           Added HIGHBAD keyword, W. Landsman January, 1991 
           Fixed occasional problem with integer inputs    W. Landsman  Feb, 1994 
           Avoid possible 16 bit integer overflow   W. Landsman  November 2001 
           Added READNOISE, NSKY keywords,  new median computation 
                              W. Landsman   June 2004 
           Added INTEGER keyword W. Landsman July 2004 
           Improve numerical precision  W. Landsman  October 2004 
           Fewer aborts on strange input sky histograms W. Landsman October 2005 
           Added /SILENT keyword  November 2005 
           Fix too many /CON keywords to MESSAGE  W.L. December 2005 
           Fix bug introduced June 2004 removing outliers when READNOISE not set 
             N. Cunningham/W. Landsman  January 2006 
           Make sure that MESSAGE never aborts  W. Landsman   January 2008 
           Add mxiter keyword and change default to 50  W. Landsman August 2011 
           Added MINSKY keyword W.L. December 2011 
           Always return floating point sky mode  W.L.  December 2015 
    """ 
     
    nsky = np.array(sky_vector).size  # Get number of sky elements 
     
    if nsky < minsky: 
        raise ValueError('ERROR -Input vector must contain at least '+str(minsky)+' elements')
     
    nlast = nsky-1  # Subscript of last pixel in SKY array
    if debug:
        print('Processing '+str(nsky) + ' element array')
    if ~integer: 
        integer = str(sky_vector.ravel()[0]).isnumeric()
        
    sky = np.sort(sky_vector.flatten()) # Sort SKY in ascending values 
     
    skymid = 0.5*sky[(nsky-1)//2] + 0.5*sky[nsky//2]  # Median value of all sky values 
     
    cut1 = np.min( np.array([skymid-sky[0],sky[nsky-1] - skymid]) ) 
    if highbad is not None:
        cut1 = np.minimum(cut1, (highbad - skymid))
    cut2 = skymid + cut1 
    cut1 = skymid - cut1 
     
    # Select the pixels between Cut1 and Cut2 
    good, = np.where( (sky <= cut2) & (sky >= cut1))
    ngood = len(good)
    if ( ngood == 0 ): 
        raise ValueError('ERROR - No sky values fall within ' + str(cut1)+' and ' + str(cut2))
     
    delta = sky[good] - skymid   # Subtract median to improve arithmetic accuracy 
    sm = np.sum(delta)
    sumsq = np.sum(delta**2) 
     
    maximm = np.max(good)  # Highest value accepted at upper end of vector 
    minimm = np.min(good)
    minimm = minimm -1   # Highest value reject at lower end of vector 
     
    # Compute mean and sigma (from the first pass). 
    skymed = 0.5*sky[(minimm+maximm+1)//2] + 0.5*sky[(minimm+maximm)//2 + 1]  # median 
    skymn = float(sm/(maximm-minimm))  # mean 
    sigma = np.sqrt(sumsq/(maximm-minimm)-skymn**2)  # sigma 
    skymn = skymn + skymid  # Add median which was subtracted off earlier 
     
    # If mean is less than the mode, then the contamination is slight, and the 
    # mean value is what we really want. 
    if skymed < skymn:
        skymod = 3.*skymed - 2.*skymn
    else:
        skymod = skymn
     
    # Rejection and recomputation loop: 
    niter = 0 
    clamp = 1 
    old = 0 
    redo = True
    while (redo):
        niter += 1 
        if ( niter > maxiter ): 
            sigma = -1.0
            skew = 0.0 
            raise ValueError('ERROR - Too many ('+str(maxiter)+') iterations, unable to compute sky')
        if ( maximm-minimm < minsky ): # Error?
            raise ValueError('RROR - Too few ('+str(maximm-minimm)+') valid sky elements, unable to compute sky')
     
        # Compute Chauvenet rejection criterion. 
        r = np.log10( float( maximm-minimm ) ) 
        r = np.max( [ 2., ( -0.1042*r + 1.1695)*r + 0.8895 ] ) 
     
        # Compute rejection limits (symmetric about the current mode). 
        cut = r*sigma + 0.5*np.abs(skymn-skymod) 
        if integer : 
            cut = cut > 1.5 
        cut1 = skymod - cut
        cut2 = skymod + cut 

        # 
        # Recompute mean and sigma by adding and/or subtracting sky values 
        # at both ends of the interval of acceptable values. 
        redo = False
        newmin = minimm 
        tst_min = (sky[newmin+1] >= cut1)  # Is minimm+1 above current CUT? 
        done = (newmin == -1) and tst_min  # Are we at first pixel of SKY? 
        if ~done : 
            done = (sky[np.maximum(newmin,0)] < cut1) and tst_min 
        if ~done: 
            istep = 1 - 2*int(tst_min) 
            while (done==False):
                newmin += istep 
                done = (newmin == -1) or (newmin == nlast) 
                if ~done : 
                    done = (sky[newmin] <= cut1) and (sky[newmin+1] >= cut1) 

            if tst_min: 
                delta = sky[newmin+1:minimm+1] - skymid 
            else: 
                delta = sky[minimm+1:newmin+1] - skymid 
            sm = sm - istep*np.sum(delta) 
            sumsq = sumsq - istep*np.sum(delta**2) 
            redo = True
            minimm = newmin 
        # 
        newmax = maximm 
        tst_max = (sky[maximm] <= cut2)       # Is current maximum below upper cut? 
        done = (maximm == nlast) and tst_max  # Are we at last pixel of SKY array? 
        if ~done: 
            done = ( tst_max ) and (sky[np.minimum((maximm+1),nlast)] > cut2) 
        if ~done:  # Keep incrementing NEWMAX 
            istep = -1 + 2*int(tst_max)  # Increment up or down? 
            while (done==False):
                newmax += istep 
                done = (newmax == nlast) or (newmax == -1) 
                if ~done : 
                    done = ( sky[newmax] <= cut2 ) and ( sky[newmax+1] >= cut2 ) 
            if tst_max: 
                delta = sky[maximm+1:newmax+1] - skymid 
            else: 
                delta = sky[newmax+1:maximm+1] - skymid 
            sm = sm + istep*np.sum(delta) 
            sumsq = sumsq + istep*np.sum(delta**2) 
            redo = True
            maximm = newmax 

        # 
        # Compute mean and sigma (from this pass). 
        # 
        nsky = maximm - minimm 
        if ( nsky < minsky ):   # Error? 
            raise ValueError('ERROR - Outlier rejection left too few sky elements')
 
        skymn = float(sm/nsky) 
        sigma = float( np.sqrt( np.maximum((sumsq/nsky - skymn**2),0) )) 
        skymn = skymn + skymid 
 
        #  Determine a more robust median by averaging the central 20% of pixels. 
        #  Estimate the median using the mean of the central 20 percent of sky 
        #  values.   Be careful to include a perfectly symmetric sample of pixels about 
        #  the median, whether the total number is even or odd within the acceptance 
        #  interval 
 
        center = (minimm + 1 + maximm)/2. 
        side = int(np.round(0.2*(maximm-minimm)))/2.  + 0.25 
        j = int(np.round(center-side))
        k = int(np.round(center+side))
 
        #  In case  the data has a large number of of the same (quantized) 
        #  intensity, expand the range until both limiting values differ from the 
        #  central value by at least 0.25 times the read noise. 
 
        if readnoise is not None:
            l = int(np.round(center-0.25))
            m = int(np.round(center+0.25)) 
            r = 0.25*readnoise 
            while ((j > 0) and (k < nsky-1) and ( ((sky[l] - sky[j]) < r) or ((sky[k] - sky[m]) < r))): 
                j -= 1 
                k += 1 
        skymed = np.sum(sky[j:k+1])/(k-j+1) 
            
        #  If the mean is less than the median, then the problem of contamination 
        #  is slight, and the mean is what we really want. 
 
        if skymed < skymn:
            dmod = 3.*skymed-2.*skymn-skymod 
        else:
            dmod = skymn - skymod
 
        # prevent oscillations by clamping down if sky adjustments are changing sign 
        if dmod*old < 0: 
            clamp = 0.5*clamp 
        skymod += clamp*dmod 
        old = dmod 
 
    # 
    skew = float( (skymn-skymod)/np.max([1.,sigma]) ) 
    nsky = maximm - minimm 
 

    if debug:
        print( '% MMM: Number of unrejected sky elements: '+str(nsky)+'    Number of iterations: '+str(niter)) 
        print( '% MMM: Mode, Sigma, Skew of sky vector:'+str(skymod)+','+str(sigma)+','+str(skew))
 
    return skymod, sigma, skew, nsky
 
