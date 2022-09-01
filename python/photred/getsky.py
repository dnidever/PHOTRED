#!/usr/bin/env python

import os
import time
import numpy as np

def photred_sky,image,skymode,skysig, SILENT=silent, CIRCLERAD = circlerad,       _EXTRA = _EXTRA, NAN = nan, MEANBACK = meanback, HIGHBAD = highbad,      histback=histback 
    #+ 
    # NAME: 
    #       PHOTRED_SKY 
    # PURPOSE: 
    #       Determine the sky level in an image 
    # EXPLANATION: 
    #       Approximately 10000 uniformly spaced pixels are selected for the 
    #       computation.  Adapted from the DAOPHOT routine of the same name. 
    # 
    #       The sky is computed either by using the procedure mmm.pro (default) 
    #       or by sigma clipping (if /MEANBACK is set) 
    # 
    # CALLING SEQUENCE: 
    #       SKY, image, [ skymode, skysig ,/SILENT, /MEANBACK, /NAN, CIRCLERAD= ] 
    # 
    #         Keywords available  when MEANBACK is not set (passed to mmm.pro): 
    #                   /DEBUG, HIGHBAD=, /INTEGER, MAXITER=. READNOISE= 
    #         Keywords available when /MEANBACK is set: 
    #                   CLIPSIG=, /DOUBLE, CONVERGE_NUM=, MAXITER=, /VERBOSE 
    # INPUTS: 
    #       IMAGE - One or two dimensional array 
    # 
    # OPTIONAL OUTPUT ARRAYS: 
    #       SKYMODE - Scalar, giving the mode of the sky pixel values of the 
    #               array IMAGE, as determined by the procedures MMM or MEANCLIP 
    #       SKYSIG -  Scalar, giving standard deviation of sky brightness.   If it 
    #               was not possible to derive a mode then SKYSIG is set to -1 
    # 
    # INPUT KEYWORD PARAMETERS: 
    #	CIRCLERAD - Use this keyword to have SKY only select pixels within 
    #		specified pixel radius of the center of the image.  If 
    #		CIRCLERAD =1, then the radius is set equal to half the image 
    #		width.   Can only be used with square images. 
    #       /MEANBACK - if set, then the background is computed using the 3 sigma 
    #             clipped mean (using meanclip.pro) rather than using the mode 
    #             computed with mmm.pro.    This keyword is useful for the Poisson 
    #             count regime or where contamination is known  to be minimal. 
    #       /HISTBACK - use mode from histogram to find background value. 
    #       /NAN - This keyword must be set to  ignore NaN values when computing 
    #              the sky. 
    #       /SILENT - If this keyword is supplied and non-zero, then SKY will not 
    #               display the sky value and sigma at the terminal 
    # 
    #      The _EXTRA facility can is used to pass optional keywords to the programs 
    #             that actually perform the sky computation: either mmm.pro 
    #             (default) or meanclip.pro (if /MEANBACK) is set.    The following 
    #             keywords are available with the mmm.pro (default) setting 
     
    #       HIGHBAD - scalar value of the (lowest) "bad" pixel level (e.g. cosmic 
    #                rays or saturated pixels) If not supplied, then there is 
    #                assumed to be no high bad pixels. 
    #       READNOISE - Scalar giving the read noise (or minimum noise for any 
    #                pixel).     Normally, MMM determines the (robust) median by 
    #                averaging the central 20% of the sky values.     In some cases 
    #                where the noise is low, and pixel values are quantized a 
    #                larger fraction may be needed.    By supplying the optional 
    #                read noise parameter, MMM is better able to adjust the 
    #                fraction of pixels used to determine the median. 
    #       /INTEGER - Set this keyword if the  input SKY image only contains 
    #                discrete integer values.    This keyword is only needed if the 
    #                SKY image is of type float or double precision, but contains 
    #                only discrete integer values. 
    # 
    #     If the /MEANBACK keyword is set then the following keywords are available 
    # 
    #       CLIPSIG:  Number of sigma at which to clip.  Default=3 
    #	MAXITER:  Ceiling on number of clipping iterations.  Default=5 
    #       CONVERGE_NUM:  If the proportion of rejected pixels is less 
    #           than this fraction, the iterations stop.  Default=0.02, i.e., 
    #           iteration stops if fewer than 2% of pixels excluded. 
    #       /DOUBLE - if set then perform all computations in double precision. 
    #                 Otherwise double precision is used only if the input 
    #                 data is double 
    # 
    # PROCEDURE: 
    #       A grid of points, not exceeding 10000 in number, is extracted 
    #       from the srray.  The mode of these pixel values is determined 
    #       by the procedure mmm.pro or meanclip.pro.   In a 2-d array the grid is 
    #       staggered in each row to avoid emphasizing possible bad columns 
    # 
    # PROCEDURE CALLS: 
    #       MEANCLIP, MMM, DIST_CIRCLE 
    # REVISION HISTORY: 
    #       Written, W. Landsman   STX Co.            September, 1987 
    #       Changed INDGEN to LINDGEN                 January, 1994 
    #       Fixed display of # of points used         March, 1994 
    #       Stagger beginning pixel in each row, added NSKY, READNOISE, HIGHBAD 
    #          W. Landsman        June 2004 
    #      Adjustments for unbiased sampling  W. Landsman June 2004 
    #      Added /NAN keyword, put back CIRCLERAD keyword W. Landsman July 2004 
    #      Added MEANBACK keyword, _EXTRA kewyord ,preserve data type in 
    #             calculations       W. Landsman November 2005 
    #      Fix problem for very large images by requiring at least 2 pixels to 
    #       be sampled per row.    March 2007    W. Landsman 
    #      Avoid possible out of bounds if /NAN set   W. Landsman   Jan 2008 
    #      Use  TOTAL(/INTEGER)      June 2009 
    #      Fix occasional out of bounds problem when /NAN set W. Landsman Jul 2013 
    #      Use HIGHBAD in selecting data points  W. Landsman  Nov 2016 
    #      Fixed a bug and made photred_sky.pro version  D.Nidever  Jan 2019 
    #- 
    On_error,2#Return to caller 
    compile_opt idl2 
     
    if N_params() == 0: 
        print('Syntax - photred_sky, image, [ skymode, skysig , HIGHBAD= ' 
        print( '                    READNOISE = , /NAN, CIRCLERAD = , /SILENT ]' 
        return 
     
    checkbad = (len(highbad) > 0) || keyword_set(circlerad) ||               keyword_set(nan) 
    s = size(image) 
    nrow = s[1] 
    if s[0] == 1: 
        ncol = 1 
    else: 
        begin 
    if s[0] != 2 : 
        message,           'ERROR - Input array (first parameter) must be 1 or 2 dimensional' 
    ncol = s[2] 
if keyword_set(circlerad) : 
    if ncol != nrow : 
        message,        'ERROR - The CIRCLERAD keyword only applies to a 2-d square array' 
 
if checkbad: 
    mask = replicate(1b, nrow, ncol) 
    if len(highbad) > 0 : 
        mask = mask and (image < highbad) 
    if keyword_set(nan) : 
        mask = mask and finite(image) 
    if keyword_set(circlerad): 
        if circlerad == 1: 
            rad = nrow/2 
        else: 
            rad = int(circlerad) 
        dist_circle,drad, nrow 
        mask = mask and (temporary(drad) < rad) 
    npts = np.sum(mask,/integer) 
 else  npts = len(image) 
     
    #  Use ~10000 data points or  at least 2 points per row 
    maxsky = 2*npts/(nrow-1) > 10000#Maximum # of pixels to be used in sky calculation 
    # Maintain the same data type as the input image Nov 2005 
    istep = npts/maxsky +1 
    skyvec = make_array(maxsky+500,type=size(image,/type)) 
    #skyvec = make_array(maxsky+200,type=size(image,/type)) 
    nstep = (nrow/istep) 
     
    jj = 0 
    index0 = istep*lindgen(nstep) 
    if nstep > 1: 
        i0 = (nrow-1 - max(index0)  - istep)/2 > 0#Adjust margin for symmetry 
        index0  = index0 + i0 
     
    # The beginning index in each row is staggered to avoid emphasizing possible 
    # bad columns 
     
    for i in range( Ncol): 
        index  = index0 + (i mod istep) 
        row = image[:,i] 
        if checkbad: 
            g , = np.where(mask[:,i],ng) 
            case ng of 
                0: goto,:ne 
                Nrow: 
                else: row = row[g] 
 
         else ng = nrow 
            imax = value_locate( index, ng-1) > 0 
            ix = index[0:imax] < (ng-1) 
            skyvec[jj] = row[ix] 
            jj = jj + imax + 1 
            if jj > maxsky : 
                break 
           :NE: 
             
        skyvec = skyvec[0:jj-1] 
         
         
        if keyword_set(meanback): 
            meanclip, skyvec, skymode, skysig,sub=sub, _EXTRA = _extra 
            nsky = len(sub) 
         else $ MMM, skyvec, skymode, skysig, _EXTRA = _extra, nsky = nsky 
             
            # Use histogram around median to get mode 
            if keyword_set(histback): 
                gd , = np.where(abs(image-skymode) < 4*skysig,ngd) 
                hist = histogram(image[gd],bin=skysig/40,locations=xhist) 
                xhist2 = scale_vector(findgen(1000),min(xhist),max(xhist)) 
                interp,xhist,hist,xhist2,hist2 
                bestind = first_el(maxloc(hist2)) 
                skymode1 = skymode# save original one 
                skymode = xhist2[bestind] 
             
             
            skymode = float(skymode)  &  skysig = float(skysig) 
            if ~keyword_set(SILENT): 
                print('Number of points used to find sky = ',nsky 
                print('Approximate sky value for this frame = ',skymode 
                print('Standard deviation of sky brightness = ',skysig 
             
            return 
 
        