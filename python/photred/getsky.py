#!/usr/bin/env python

import os
import time
import numpy as np
from . import utils


def getsky(image,silent=False,circlerad=False,meanback=False,highbad=None,
           histback=False,nan=False):
    """
           Determine the sky level in an image 
     EXPLANATION: 
           Approximately 10000 uniformly spaced pixels are selected for the 
           computation.  Adapted from the DAOPHOT routine of the same name. 
     
           The sky is computed either by using the procedure mmm.pro (default) 
           or by sigma clipping (if /MEANBACK is set) 
     
     CALLING SEQUENCE: 
           SKY, image, [ skymode, skysig ,/SILENT, /MEANBACK, /NAN, CIRCLERAD= ] 
     
             Keywords available  when MEANBACK is not set (passed to mmm.pro): 
                       /DEBUG, HIGHBAD=, /INTEGER, MAXITER=. READNOISE= 
             Keywords available when /MEANBACK is set: 
                       CLIPSIG=, /DOUBLE, CONVERGE_NUM=, MAXITER=, /VERBOSE 
     INPUTS: 
           IMAGE - One or two dimensional array 
     
     OPTIONAL OUTPUT ARRAYS: 
           SKYMODE - Scalar, giving the mode of the sky pixel values of the 
                   array IMAGE, as determined by the procedures MMM or MEANCLIP 
           SKYSIG -  Scalar, giving standard deviation of sky brightness.   If it 
                   was not possible to derive a mode then SKYSIG is set to -1 
     
     INPUT KEYWORD PARAMETERS: 
    	CIRCLERAD - Use this keyword to have SKY only select pixels within 
    		specified pixel radius of the center of the image.  If 
    		CIRCLERAD =1, then the radius is set equal to half the image 
    		width.   Can only be used with square images. 
           /MEANBACK - if set, then the background is computed using the 3 sigma 
                 clipped mean (using meanclip.pro) rather than using the mode 
                 computed with mmm.pro.    This keyword is useful for the Poisson 
                 count regime or where contamination is known  to be minimal. 
           /HISTBACK - use mode from histogram to find background value. 
           /NAN - This keyword must be set to  ignore NaN values when computing 
                  the sky. 
           /SILENT - If this keyword is supplied and non-zero, then SKY will not 
                   display the sky value and sigma at the terminal 
     
          The _EXTRA facility can is used to pass optional keywords to the programs 
                 that actually perform the sky computation: either mmm.pro 
                 (default) or meanclip.pro (if /MEANBACK) is set.    The following 
                 keywords are available with the mmm.pro (default) setting 
     
           HIGHBAD - scalar value of the (lowest) "bad" pixel level (e.g. cosmic 
                    rays or saturated pixels) If not supplied, then there is 
                    assumed to be no high bad pixels. 
           READNOISE - Scalar giving the read noise (or minimum noise for any 
                    pixel).     Normally, MMM determines the (robust) median by 
                    averaging the central 20% of the sky values.     In some cases 
                    where the noise is low, and pixel values are quantized a 
                    larger fraction may be needed.    By supplying the optional 
                    read noise parameter, MMM is better able to adjust the 
                    fraction of pixels used to determine the median. 
           /INTEGER - Set this keyword if the  input SKY image only contains 
                    discrete integer values.    This keyword is only needed if the 
                    SKY image is of type float or double precision, but contains 
                    only discrete integer values. 
     
         If the /MEANBACK keyword is set then the following keywords are available 
     
           CLIPSIG:  Number of sigma at which to clip.  Default=3 
    	MAXITER:  Ceiling on number of clipping iterations.  Default=5 
           CONVERGE_NUM:  If the proportion of rejected pixels is less 
               than this fraction, the iterations stop.  Default=0.02, i.e., 
               iteration stops if fewer than 2% of pixels excluded. 
           /DOUBLE - if set then perform all computations in double precision. 
                     Otherwise double precision is used only if the input 
                     data is double 
     
     PROCEDURE: 
           A grid of points, not exceeding 10000 in number, is extracted 
           from the srray.  The mode of these pixel values is determined 
           by the procedure mmm.pro or meanclip.pro.   In a 2-d array the grid is 
           staggered in each row to avoid emphasizing possible bad columns 
     
     PROCEDURE CALLS: 
           MEANCLIP, MMM, DIST_CIRCLE 
     REVISION HISTORY: 
           Written, W. Landsman   STX Co.            September, 1987 
           Changed INDGEN to LINDGEN                 January, 1994 
           Fixed display of # of points used         March, 1994 
           Stagger beginning pixel in each row, added NSKY, READNOISE, HIGHBAD 
              W. Landsman        June 2004 
          Adjustments for unbiased sampling  W. Landsman June 2004 
          Added /NAN keyword, put back CIRCLERAD keyword W. Landsman July 2004 
          Added MEANBACK keyword, _EXTRA kewyord ,preserve data type in 
                 calculations       W. Landsman November 2005 
          Fix problem for very large images by requiring at least 2 pixels to 
           be sampled per row.    March 2007    W. Landsman 
          Avoid possible out of bounds if /NAN set   W. Landsman   Jan 2008 
          Use  TOTAL(/INTEGER)      June 2009 
          Fix occasional out of bounds problem when /NAN set W. Landsman Jul 2013 
          Use HIGHBAD in selecting data points  W. Landsman  Nov 2016 
          Fixed a bug and made photred_sky.pro version  D.Nidever  Jan 2019 
    """ 

    if image.ndim not in [1,2]:
        raise ValueError('ERROR - Input array (first parameter) must be 1 or 2 dimensional')
        
    checkbad = ( (highbad is not None) or circlerad or nan)
    sh = image.shape
    if image.ndim==1:
        ncol = 1
        nrow = image.size
    else:
        ncol,nrow = image.shape

    if circlerad:
        if ncol != nrow: 
            raise ValueError('ERROR - The CIRCLERAD keyword only applies to a 2-d square array')
 
    if checkbad: 
        mask = np.ones(image.shape)
        if highbad is not None:
            mask = mask and (image < highbad) 
        if nan: 
            mask = mask and finite(image) 
        if circlerad: 
            if circlerad == 1: 
                rad = nrow/2 
            else: 
                rad = int(circlerad) 
            # Make image where each value is its distance to a given center
            xv,yv = np.meshgrid(np.arange(nrow),np.arange(nrow))
            cen = (nrow-1)/2.
            drad = np.sqrt((xv-cen)**2+(yv-cen)**2)            
            #dist_circle,drad, nrow 
            mask = mask and (drad < rad) 
        npts = np.sum(mask)
    else:
        npts = image.size
     
    #  Use ~10000 data points or  at least 2 points per row 
    maxsky = np.maximum(2*npts//(nrow-1), 10000)  # Maximum # of pixels to be used in sky calculation 
    # Maintain the same data type as the input image Nov 2005 
    istep = npts//maxsky +1
    skyvec = np.zeros(maxsky+500,dtype=image.dtype)
    #skyvec = make_array(maxsky+200,type=size(image,/type)) 
    nstep = (nrow//istep) 
     
    jj = 0 
    index0 = istep*np.arange(nstep) 
    if nstep > 1: 
        i0 = np.maximum((nrow-1 - max(index0)  - istep)//2, 0)  # Adjust margin for symmetry 
        index0  = index0 + i0 
     
    # The beginning index in each row is staggered to avoid emphasizing possible 
    # bad columns 
     
    for i in range(ncol): 
        index = index0 + (i % istep) 
        row = image[i,:] 
        if checkbad: 
            g, = np.where(mask[i,:]) 
            ng = len(g)
            if ng==0:
                break
            row = row[g] 
        else:
            ng = nrow 
        imax = np.maximum(np.searchsorted(index, ng-1), 0)
        #imax = value_locate( index, ng-1) > 0 
        ix = np.minimum( index[0:imax], ng-1)
        skyvec[jj:jj+len(ix)] = row[ix] 
        jj += imax
        if jj > maxsky: 
            break 

    skyvec = skyvec[0:jj] 

    meanback = True
         
    if meanback: 
        skymode, skysig, subs = utils.meanclip(skyvec)
        nsky = len(subs) 
    else:
        MMM, skyvec, skymode, skysig, _EXTRA = _extra, nsky = nsky 
             
    # Use histogram around median to get mode 
    if histback:
        gd = (np.abs(image-skymode) < 4*skysig) 
        xhist = np.arange(np.min(image[gd]),np.max(image[gd]),skysig/40)
        hist,bin_edges = np.histogram(image[gd],bins=xhist)
        xhist2 = np.linspace(np.min(xhist),np.max(xhist),1000)
        hist2 = np.interp(xhist2,xhist[:-1],hist)
        bestind = np.argmax(hist2)
        skymode1 = np.copy(skymode)  # save original one 
        skymode = xhist2[bestind] 
             
             
    skymode = float(skymode)
    skysig = float(skysig) 
    if silent==False:
        print('Number of points used to find sky = ',nsky)
        print('Approximate sky value for this frame = ',skymode)
        print('Standard deviation of sky brightness = ',skysig)
        
    return skymode,skysig
 
        
