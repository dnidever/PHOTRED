#!/usr/bin/env python

import os
import time
import numpy as np
from astropy.io import fits
from dlnpyutils import utils as dln
from . import utils,io,sky
 
def get_subim(im,xind,yind,hwidth,sky=0.0):
     
    # This function returns a subimage 
    #  hwidth is the half-width.  Total width = 2*hwidth+1 
     
    # Getting the image size
    ny,nx = im.shape
     
    # Initializing the subimage
    subim = np.zeros(2*hwidth+1,2*hwidth+1,float) + sky
     
    # Checking that we're getting any part of the image 
    # The center can be off the edge 
    if ( (xind+hwidth) >= 0 ) and ( (xind-hwidth) <= (nx-1) )    and ( (yind+hwidth) >= 0 ) and ( (yind-hwidth) <= (ny-1) ): 
         
        # Indices for the big image 
        xbg0 = np.maximum((xind-hwidth), 0)
        xbg1 = np.minimum((xind+hwidth), nx)
        ybg0 = np.maximum((yind-hwidth), 0)
        ybg1 = np.minimum((yind+hwidth), ny)
         
        # Indices for the subimage 
        xsm0 = hwidth+xbg0-xind 
        xsm1 = hwidth+xbg1-xind 
        ysm0 = hwidth+ybg0-yind 
        ysm1 = hwidth+ybg1-yind 
         
        # Getting part of the image 
        subim[ysm0:ysm1,xsm0:xsm1] = im[ybg0:ybg1,xbg0:xbg1] 
     
    return subim 
     
 
def get_fluxcenter(subim): 
     
    # Calculate the flux-weighted center 
     
    # Center-of-Mass like centroid.  Weight by flux 
    # Add up all position along columns since they will have the same 
    # X value. Similar for Y
    ny,nx = subim.shape
    mask = (subim >= 0.0).astype(float)  # mask out negative pixels
    xind = np.sum( np.sum(subim*mask,axis=2)*np.arange(nx) )/np.sum(subim*mask) 
    yind = np.sum( np.sum(subim*mask,axis=1)*np.arange(ny) )/np.sum(subim*mask) 
     
    return xind,yinb
 
 
def imfwhm(inpfiles=None,outfile=None,exten=None,im=None,head=None,
           skymode=None,skysig=None,backgim=None,nsigdetect=8,
           verbose=True):
    """ 
    The program estimates the FWHM of stars in images. 

    Parameters
    ----------
    input      The name of the FITS image file. Can be a glob, 
               i.e. '*.fits', or a list of files if it starts 
               with an '@' 
    outfile   A file to print the output to 
    nsigdetect  The source detection limit in background sigma. 
                 The default is 8. 
    verbose : boolean, option
      Print to the screen.  By default the 
        filename and FWHM are printed to the screen. 

    Returns
    -------
    fwhm : float
       The median fwhm of stars in the image.  If multiple 
         images are processed then this will be an array. 
    If "outfile" is set then the FWHM values are written to this file. 
    ellipticity  The ellipticity (1-a/b). 
    gtab : table
       The table of Gaussian fits to the sources used 
         to measure FWHM and ELLIPTICITY.  This only works 
         if there's only ONE input file. 
    peaktab : table
       The structure of detected peaks and quick shape 
               measurements. 

    Example
    -------

    fwhm,gtab = imfwhm('test.fits')
    test.fits     4.973 
 
    PROCEDURES USED: 
    sky.pro          To measure the background level (in IDL Astro User's Library) 
    meanclip.pro     Iteratively sigma-clipped mean (in IDL Astro User's Library) 
    readcol.pro      Read an ASCII text file (in IDL Astro User's Library) 
    mmm.pro          Estimate the sky backround (in IDL Astro User's Library) 
    resistant_mean.pro  Compute the mean of an array with outlier rejection (in IDL Astro User's Library) 
    loadinput.pro    Load input files (by D.Nidever) 
    robust_mean.pro  Compute the mean of an array robustly with iterative outlier rejection 
                        and input uncertainties (by D.Nidever) 
    undefine.pro     Make a variable undefined  (by D.Nidever) 
    mad.pro          A resistant Standard Deviation (by D.Nidever) 
    push.pro         Add an element to an array (by D.Nidever) 
    readline.pro     Read an ASCII file into a string array (by D.Nidever) 
    wmeanerr.pro     Compute a weighted mean and Standard Deviation (originally by M.Buie) 
 
    MODIFICATION HISTORY: 
    Created on July 22, 2006 by D.Nidever 
    August 31, 2006  D. Nidever, won't crash on FITS_READ errors anymore 
                      and checks heights of neighboring pixels, which 
                      gets rid of most cosmic rays.  Also doesn't use 
                      regions close to saturated stars. 
 
    Translated to Python by D. Nidever, May 2022
    """ 
    
    # Load the input 
    files = '' 
    nfiles = 0
    inpim = None
    if im is not None:
        inpim = im
    else:
        inpim = None
    if head is not None:
        head0 = head
    else:
        head0 = None
    if inpfiles is not None:
        files = dln.loadinput(inpfiles)
        nfiles = len(files)
     
    # Using the input image
    if im is not None and nfiles <= 1:
        if nfiles == 1: 
            print('One file input and one image (=im) input.  Using the input image')
        nfiles = 1 
        files = 'Input Image' 
        inpim = im
        
    # Multiple files AND image input with =im 
    if (nfiles > 1 and im is not None):
        print('Multiple files AND one image (=im) input.  Using the multiple files')
        im = 0
     
    # Starting ALLFWHM
    allfwhm = np.zeros(nfiles,float)+np.inf
    allellip = np.zeros(nfiles,float)+np.inf
     
    # Detection threshoold
    nsigdetect = np.maximum(nsigdetect,0)
      
    # Opening output file
    if outfile is not None:
        fout = open(outfile,'w')

    # Looping through the files 
    for f in range(nfiles):         
        fwhm = 99.99 # bad until proven good 
        ellipticity = 99.99 

        # Loading image from file
        if inpim is None:
            # Test that file exists
            if os.path.exists(files[f])==False:
                if verbose:
                    print(files[f]+' NOT FOUND')
                    continue
            img,head = io.readfile(files[f],exten=exten)
             
        # Using input image 
        else:    
            files = 'Input Image' 
            img = inpim.copy()
            if head0 is not None: # input header, otherwise make fake header 
                head = head0 
            else:
                head = fits.Header()

        # We have an image 
        if (img is not None):
            # Need float 
            img = img.astype(float)             
            # Image size
            ny,nx = img.shape
            # Saturation limit 
            satlim = np.max(img)  # starting point 
            nsaturate = 0
            if head is not None:
                saturate = head.get('SATURATE')
                if saturate is not None:
                    satlim = saturate 
            if nsaturate == 0: 
                satlim = dln.limit(satlim,40000.,65000.)  # this is a realistic limit for now
            if verbose:
                print('satlim = %.1f' % satlim)
             
            # Set NAN/Inf pixels to above the saturation limit 
            bdnan = (~np.isfinite(img))
            if np.sum(bdnan) > 0: 
                img[bdnan] = satlim+5000. 
             
            # Not enough "good" pixels 
            # sometimes "bad" pixels are exactly 0.0 
            gdpix = ((img < satlim*0.90) & (img != 0.0))
            if np.sum(gdpix) < 2:
                if verbose:
                    print(files[f]+' NOT ENOUGH GOOD PIXELS')
                fwhm = 99.99 
                ellipticity = 99.99
                continue
             
            # Computing sky level and sigma
            skymode,skysig1 = sky.getsky(img,highbad=satlim*0.95,silent=True)
            if skysig1 < 0.0: 
                skysig1 = dln.mad(img[gdpix])
            if skysig1 < 0.0: 
                skysig1 = dln.mad(img) 
            maxim = np.max(img) 

            import pdb; pdb.set_trace()

             
            #-- Compute background image -- 
             
            # First pass, no clipping (except for saturated pixels) 
            backgim_large = img.copy()
            # Set saturated pixels to NaN so they won't be used in the smoothing 
            bdpix = (backgim_large > satlim*0.95) 
            if np.sum(bdpix) > 0: 
                backgim_large[bdpix] = np.nan 
            sm = np.minimum(np.minimum(400,(nx//2)),(ny//2))
            #backgim_large = dln.smooth(backgim_large,[sm,sm],/edge_truncate,/nan,missing=skymode) 
            backgim_large = dln.smooth(backgim_large,sm,fillvalue=skymode) 
            backgim_large[~np.isfinite(backgim_large)] = skymode

            # Second pass, use clipping, and first estimate of background 
            backgim1 = img.copy()
            # Setting hi/low pixels to NaN, they won't be used for the smoothing 
            #bd = where(abs(im-skymode) gt 2.0*skysig1,nbd) 
            bd1 = (np.abs(backgim1-backgim_large) > 3.0*skysig1)
            if np.sum(bd1) > 0: 
                backgim1[bd1] = np.nan
            sm = np.minimum(np.minimum(400, nx/2.0),ny/2.0)
            #backgim1 = dln.smooth(backgim1,[sm,sm],/edge_truncate,/nan,missing=skymode) 
            backgim1 = dln.smooth(backgim1,[sm,sm],missing=skymode) 
             
            # Third pass, use better estimate of background 
            backgim2 = img.copy()
            # Setting hi/low pixels to NaN, they won't be used for the smoothing 
            #bd = where(abs(im-skymode) gt 2.0*skysig1,nbd) 
            bd2 = (np.abs(backgim2-backgim1) > 3.0*skysig1) 
            if len(bd2) > 0: 
                backgim2[bd2] = np.inf
            sm = np.minimum(np.minimum(400,nx/2.0),ny/2.0)
            #backgim2 = dln.smooth(backgim2,[sm,sm],/edge_truncate,/nan,missing=skymode) 
            backgim2 = dln.smooth(backgim2,[sm,sm],missing=skymode) 
             
            backgim = backgim2.copy()
             
            # Setting left-over NaNs to the skymode 
            b = n(~np.isfinite(backgim)) 
            if np.sum(b)>0 : 
                backgim[b] = skymode 
             
            # Creating background subtracted image 
            img2 = img - backgim 
             
            # Mask out bad pixels, set to background 
            bpmask = (img2 >= (0.9*satlim)).astype(float)  # making bad pixel mask 
            bpmask_orig = bpmask 
            satpix = (img2 >= 0.9*satlim) 
            if np.sum(satpix) > 0 : 
                img2[satpix] = backgim[satpix]
             
            # Computing sky level and sigma AGAIN with 
            #  background subtracted image
            skymode2,skysig = getsky(img2,highbad=satlim*0.95,verbose=False)
            if verbose:
                print('skymode = %.2f skysig = %.2f' % (skymode,skysig))
             
            # Gaussian smooth the image to allow detection of fainter sources
            gx = np.arange(5).reshape(-1,1) + np.zeros(5).reshape(1,-1)
            gy = np.zeros(5).reshape(-1,1) + np.arange(5).reshape(1,-1)
            gkernel = np.exp(-0.5*( (gx-2.0)**2.0 + (gy-2.0)**2.0 )/1.0**2.0 ) 
            gkernel = gkernel/np.sum(gkernel) 
            smim = CONVOL(img2,gkernel) #,/center,/edge_truncate) 
             
            # Getting maxima points 
            diffx1 = smim-shift(smim,1,0) 
            diffx2 = smim-shift(smim,-1,0) 
            diffy1 = smim-shift(smim,0,1) 
            diffy2 = smim-shift(smim,0,-1) 
             
            # Get the gain (electrons/ADU)
            gain = head.get('GAIN')
            if gain is None:
                # use scatter in background and Poisson statistic 
                # to calculate gain empirically 
                # Nadu = Ne/gain 
                # Poisson scatter in electrons = sqrt(Ne) = sqrt(Nadu*gain) 
                # scatter (ADU) = scatter(e)/gain = sqrt(Nadu*gain)/gain = sqrt(Nadu/gain) 
                # gain = Nadu / scatter(ADU)^2 
                 
                # of course we should remove the RDNOISE in quadrature from the empirical scatter 
                rdnoise = head.get('RDNOISE')  # normally in electrons
                if rdnoise is not None:
                    rdnoise_adu = 0.0 
                if nrdnoise > 0 and ngain > 0 : 
                    rdnoise_adu = rdnoise/gain 
                skyscatter = np.sqrt( skysig**2 - rdnoise_adu**2) 
                gain = np.median(backgim)/skyscatter**2 
            if gain < 0: 
                gain = 1.0
            if verbose:
                print('gain = %.2f' % gain)
             
            # Make the SIGMA map 
            #sigmap = sqrt(backgim>1) > skysig 
            sigmap = np.maximum(np.sqrt(np.maximum(img/gain,1)), skysig)
            sigmap = sigmap*(1.0-bpmask) + bpmask*65000. 
             
             
            # Getting the "stars" 
            # Must be a maximum, 8*sig above the background, but 1/2 the maximum (saturation) 
            niter = 0 
            detendflag = False
            while (detendflag==False):
                if niter > 0 :# restore bpmask since we modify it below to specify these pixels as "done" 
                    bpmask = bpmask_orig.copy()
                diffth = 0.0#skysig  ; sigmap 
                ind, = np.where((diffx1 > diffth) & (diffx2 > diffth) & (diffy1 > diffth) & (diffy2 > diffth) &
                                (im2 > (nsig*sigmap)) & (im < 0.5*satlim))
                # No "stars" found, try lower threshold 
                if nind < 2 and niter < 5: 
                    nsig *= 0.7  # smaller drop 
                    if verbose:
                        print('No good sources detected.  Lowering detection limts to '+str(nsig)+' sigma')
                    detendflag = False
                else:
                    detendflag = True
                niter += 1 

            # No "stars" found, giving up 
            if nind < 2:
                if verbose:
                    print(files[f]+' NO STARS FOUND. GIVING UP!')
                fwhm = 99.99 
                ellipticity = 99.99 
                continue
            if verbose:
                print(str(nind)+' initial peaks found')
            ind2 = array_indices(img,ind) 
            xind = reform(ind2[0,:]) 
            yind = reform(ind2[1,:]) 
             
            offL = 50  # 1/2 width of large subimage 
            offs = 10  # 1/2 width of small subimage 
            minfrac = 0.5  # minimum fraction that neighbors need to be of the central pixel 
             
            # Creating peak structure
            dt = [('backg',float),('flux',float),('fwhm',float),('xcen',float),('ycen',float),
                  ('round',float),('max',float),('elip',float),('nbelow',float)]
            peaktab = np.zeros(nind,dtype=np.dtype(dt))
            peaktab = Table(peaktab)
             
            # Loop through the peaks 
            for i in range(nind): 
                # Checking neighboring pixels 
                # Must >50% of central pixel 
                cen = img2[yind[i],xind[i]] 
                xlo = np.maximum(xind[i]-1,0)
                xhi = np.minimum(xind[i]+1,nx)
                ylo = np.maximum(yind[i]-1,0)
                yhi = np.minimum(yind[i]+1,ny)
                nbrim = img2[ylo:yhi,xlo:xhi] 
                lofrac = np.min(nbrim/cen) 
                nbelow = np.sum((nbrim/np.max(nbrim) <= minfrac).astype(float))
                #nbelowarr[i] = nbelow 
                peaktab['nbelow'][i] = nbelow 
                #cenvalarr[i] = cen 
                 
                # Checking the bad pixel mask 
                bpmask2 = get_subim(bpmask,xind[i],yind[i],offs) 
                nbadpix = np.sum(bpmask2) 
                 
                # Smaller image 
                subims = get_subim(img2,xind[i],yind[i],5) 
                maxsubims = np.max(subims) 
                 
                maxnbrim = np.max(nbrim)   # maximum within +/-1 pixels 
                 
                # Good so far 
                #  -No saturated pixels 
                #  -Not a cosmic ray 
                #  -Must be the maximum within +/-5 pix 
                #  -Not close to the edge 
                background = backgim[yind[i],xind[i]] 
                 
                # Getting Large image 
                subimL = get_subim(img2,xind[i],yind[i],offL,background) 
                 
                # Local background in image 
                #backgarr[i] = median(subimL) 
                peaktab['backg'][i] = np.median(subimL) 
                 
                # Getting small image 
                subimS = get_subim(subimL,offL,offL,offS) 
                 
                # Getting flux center 
                xcen,ycen = get_fluxcenter(subimS)
                 
                xcen2 = int(np.round(xcen-offS))  # wrt center 
                ycen2 = int(np.round(ycen-offS))  # wrt center 
                 
                # Getting flux-centered small image 
                subim = get_subim(subimL,offL+xcen2,offL+ycen2,offS) 
                maxsubim = np.max(subim) 
                 
                # What is the flux and magnitude 
                #fluxarr[i] = total(subimS-median(subimL)) 
                peaktab['flux'][i] = np.sum(subimS-np.median(subimL)) 
                #magarr[i] = 25.0-2.5*alog10(fluxarr[i] > 1) 
                 
                # Getting the contours 
                CONTOUR(subim,levels=[maxsubim*0.5],path_xy=pathxy) #,/path_data_coords 
                 
                # Getting the path 
                xpath = reform(pathxy[0,:]) 
                ypath = reform(pathxy[1,:]) 
                xmnpath = np.mean(xpath) 
                ymnpath = np.mean(ypath) 
                 
                #peakstr[i].xcen = xcen+xlo  ; THIS WAS WRONG!! 
                #peakstr[i].ycen = ycen+ylo 
                peaktab['xcen'][i] = xcen + (xind[i]-offS) 
                peaktab['ycen'][i] = ycen + (yind[i]-offS) 
                 
                # Calculating the FWHM 
                dist = np.sqrt((xpath-xmnpath)**2.0 + (ypath-ymnpath)**2.0) 
                fwhm1 = 2.0 * np.mean(dist) 
                 
                # Measuring "ellipticity" 
                # DLN 5/9/16, added 2x factor to make it close 
                #  to "real" ellipticity (1-a/b) 
                elip = 2*np.std(dist-fwhm1)/fwhm1 
                peaktab['elip'][i] = elip 
                 
                # Putting it in the structure 
                peaktab['fwhm'][i] = fwhm1 
                peaktab['max'][i] = maxsubim 
                 
                # Computing the "round" factor 
                # round = difference of the heights of the two 1D Gaussians 
                #         ------------------------------------------------- 
                #                 average of the two 1D Gaussians 
                # 
                # Where the 1D Gaussians are of the marginal sums, i.e. sum 
                # along either the x or y dimensions 
                # round~0 is good 
                # round<0 object elongated in x-direction 
                # round>0 object elongated in y-direction 
                htx = np.max(np.sum(subim,axis=0)) 
                hty = np.max(np.sum(subim,axis=1)) 
                rnd = (hty-htx)/np.mean([htx,hty]) 
                 
                peaktab['round'][i] = rnd
                 
                # Making these pixels "bad pixels" so they won't be used again 
                xlo = np.maximum(xind[i]+xcen2-offS,0)
                xhi = np.minimum(xind[i]+xcen2+offS,nx)
                ylo = np.maximum(yind[i]+ycen2-offS,0)
                yhi = np.minimum(yind[i]+ycen2+offS,ny)
                bpmask[ylo:yhi,xlo:xhi] = 1.0 
         
        # Getting the good ones: 
        #  -If they were bad then FWHM=0 
        #  -Making sure they are "round" stars 
        #  -Ellipticity is low 
        if len(gd) == 0: 
            gd, = np.where((peaktab['fwhm'] > 0.0) & (np.abs(peaktab['round']) < 1.0))
         
        # Retry with lower detection limit 
        #if ngd lt 10 and niter lt 5 and nsig gt 2 then begin 
        if len(gd) < 50 and niter < 5 and nsig > 2: 
            nsig *= 0.7  # smaller drop
            if verbose:
                print('No good sources detected.  Lowering detection limts to '+str(nsig)+' sigma')
            niter += 1
            goto,detection 
         
        if len(gd) == 0:
            if verbose:
                print('No good sources detected')
            fwhm = 99.99 
            ellipticity = 99.99 
            continue 
        if verbose:
            print(str(ngd)+' final good peaks')
         
        # Fit Gaussians to the good sources 
        #------------------------------------
        ny,nx = img.shape
        x = np.arange(nx) 
        y = np.arange(ny)
        gdt = [('x',float),('y',float),('pars',(float,7)),('chisq',float),('dof',int),('status',int),
               ('fwhm',float),('elip',float),('round',float),('sharp',float)]
        gtab = np.zeros(ngd,dtype=np.dtype(gdt))
        gtab = Table(gtab)
        gtab['x'] = peaktab['xcen'][gd]
        gtab['y'] = peaktab['ycen'][gd]
        for i in range(ngd): 
            peaktab1 = peaktab[gd[i]] 
            ix = gtab['x'][i]
            iy = gtab['y'][i]
            xlo = np.maximum(int(np.round(ix)-10),0)
            xhi = np.minimum(int(np.round(ix)+10),nx)
            ylo = np.maximum(int(np.round(iy)-10),0)
            yhi = np.minimum(int(np.round(iy)+10),ny)
            subim = img2[ylo:yhi,xlo:xhi]  # use background subtracted image 
            xarr = x[xlo:xhi] 
            yarr = y[ylo:yhi] 
             
            parinfo = replicate({limited:[0,0],limits:[0,0],fixed:0},7) 
            parinfo[1].limited=[1,1] 
            parinfo[1].limits=[0,np.max(subim)*2]# only positive peaks 
            parinfo[4].limited=[1,1]# constrain X and Y 
            #parinfo[4].limits=[-1,1]+ix 
            parinfo[4].limits = [np.min(xarr),np.max(xarr)] 
            parinfo[5].limited=[1,1] 
            #parinfo[5].limits=[-1,1]+iy 
            parinfo[5].limits = [np.min(yarr),np.max(yarr)] 
            estimates = [peaktab1['backg'], np.maximum((peaktab1['max']-peak1['backg']), 0.5),
                         np.maximum(peaktab1['fwhm']/2.35, 0.5), np.maximum(peaktab1['fwhm/']/2.35, 0.5), ix, iy, 0.0] 
            #fit = MPFIT2DPEAK(subim,pars,xarr,yarr,estimates=estimates,chisq=chisq,dof=dof,perror=perror,/gaussian,
            #                  parinfo=parinfo,status=status) 
            gtab['status'][i] = status 
            pars,cov = curve_fit()
            if status > 0: 
                gtab['x'][i] = ix 
                gtab['y'][i] = iy 
                gtab['pars'][i] = pars 
                gtab['perror'][i] = perror 
                gtab['chisq'][i] = chisq 
                gtab['dof'][i] = dof
                gtab['fwhm'][i] = 0.5*(pars[2]+pars[3]) * 2# FWHM=average of half-widths (times 2) 
                gtab['elip'][i] = (1-np.min(pars[2:4])/np.max(pars[2:4]))# 1-a/b 
                htx = np.max(np.sum(fit-pars[0],axis=0)) 
                hty = np.max(np.sum(fit-pars[0],axis=1)) 
                rnd = (hty-htx)/np.mean([htx,hty])# see definition above 
                gtab['round'][i] = rnd 
                # This is a slightly different "sharp" since the normal 
                #  Stetson one is height of delta function / height 
                #  of symmetric Gaussian 
                gtab['sharp'][i] = peaktab1['max'] / pars[1] 
                 
                # The 2D Gaussian parameters are: 
                #   A(0)   Constant baseline level 
                #   A(1)   Peak value 
                #   A(2)   Peak half-width (x) -- gaussian sigma or half-width at half-max 
                #   A(3)   Peak half-width (y) -- gaussian sigma or half-width at half-max 
                #   A(4)   Peak centroid (x) 
                #   A(5)   Peak centroid (y) 
                #   A(6)   Rotation angle (radians) if TILT keyword set 
                 
                #display,subim,position=[0,0,0.5,1.0] 
                #display,fit,position=[0.5,0,1.0,1.0],/noerase 
                #wait,0.5 
             
         
        # Some Gaussian fits converged 
        gdgtab, = np.where(gstr['status'] > 0) 
        if len(gdgtab) > 0: 
            gtab0 = gtab 
            gtab = gtab[gdgtab]# only keep the good ones 
        else:    # none converged
            if verbose:
                print('No Gaussian fits converged')
            fwhm = 99.99 
            ellipticity = 99.99 
            continue
         
        # First make a bad pixel or CR cut on width, FWHM>1.1 pixels 
        gdgtabwid, = np.where((gtab['pars'][:,2]*2.0 > 1.1) & (gtab['pars'][:,3]*2.0 > 1.1))
        if len(gdgtabwid) > 0: 
            gtab = gtab[gdgtabwid]   # only keep ones with FWHM>1.1 pixels 
        else:
            if verbose:
                print('No Gaussian fits with FWHM>1.1 pixels')
            fwhm = 99.99 
            ellipticity = 99.99 
            continue
         
        # Now pick out the "good" ones 
        medpar2 = np.median(gtab['pars'][:,2]) 
        sigpar2 = dln.mad(gtab['pars'][:,2]) 
        medpar3 = np.median(gtab['pars'][:,3]) 
        sigpar3 = dln.mad(gtab['pars'][:,3]) 
        medchisq = np.median(gtab['chisq']) 
        sigchisq = dln.mad(gtab['chisq']) 
        #  Positive amplitude,  Low Chisq, and Half-widths within "normal" values 
        okay, = np.where((gtab['pars'][:,1] > 0.0) & (gtab['chisq'] < medchisq+3*sigchisq) &
                           (np.abs(gtab['pars'][:,2]-medpar2) < 3*sigpar2) &
                           (np.abs(gtab['pars'][:,3]-medpar3) < 3*sigpar3))
        # No stars passed, increase threshold, 4 sigma 
        if len(okay) == 0 : 
            okay, = np.where((gtab['pars'][:,1] > 0.0) & (gtab['chisq'] < medchisq+4*sigchisq) &
                             (np.abs(gtab['pars'][:,2]-medpar2) < 4*sigpar2) &
                             (np.abs(gtab['pars'][:,3]-medpar3) < 4*sigpar3))
        # No stars passed, increase threshold, 5 sigma 
        if len(okay) == 0 : 
            okay, = np.where((gtab['pars'][:,1] > 0.0) & (gtab['chisq'] < medchisq+5*sigchisq) &
                             (np.abs(gtab['pars'][:,2]-medpar2) < 5*sigpar2) &
                             (np.abs(gtab['pars'][:,3]-medpar3) < 5*sigpar3))
         
        if nokay > 0: 
            gd_orig = gd 
            gd = gd[okay] 
            ngd = nokay 
            gtab2 = gtab[okay] 
        else:
            ngd = 0 
        if verbpse:
            print(str(ngd)+' sources passed the Gaussian fitting cuts')
             
        # There are some good stars 
        if len(gd) >= 2: 
                 
            # Maybe randomly sample ~50 peaks and fit them 
            # with a Gaussian 
                 
            # Use the 2D fit values instead! 
            #fwhmarr2 = fwhmarr[gd] 
                 
            # ; Use sigma-clipped mean to get good indices 
            # meanclip,fwhmarr2,mean,sigma,subs=subs,clipsig=3 
            # 
            # ; Use the median 
            # fwhm = median(fwhmarr2(subs)) 
                 
            # Use MEDIAN 
            medfwhm = np.median(gtab2['fwhm']) 
            medellip = np.median(gtab2['elip']) 
                 
            # Use RESISTANT_MEAN 
            RESISTANT_MEAN(gtab2['fwhm'],3.0,resfwhm,fwhm_sigma)
            RESISTANT_MEAN(gtab2['elip'],3.0,resellip,ellip_sigma)
                 
            # Weighted mean (by flux) 
            #wt = fluxarr[gd]>0.0001 
            wt = np.maximum(gtab2['pars'][:,0],0.0001)
            wtfwhm = np.sum(wt*gtab2['fwhm'])/np.sum(wt) 
            wtellip = np.sum(wt*gtab2['elip'])/np.sum(wt) 
                 
            # Weighted mean with outlier rejection 
            sig = 1.0/np.sqrt(wt) 
            ROBUST_MEAN(gtab2['fwhm'],robfwhm,robfwhmsigma,sig=sig)
            ROBUST_MEAN(gtab2['elip'],robellip,robellipsigma,sig=sig)
                 
            # The four methods don't agree 
            # Use flux-weighted mean with outlier rejection. 
            #  Should be the most reliable if anything weird is going on. 
            fw = [medfwhm,resfwhm,wtfwhm,robfwhm] 
            if (np.max(fw)-np.min(fw) > 2.0): 
                if verbose:
                    print('Using Flux-weighted Mean with outlier rejection')
                fwhm = robfwhm 
                ellipticity = robellip          
            # Use resistant mean 
            else: 
                fwhm = resfwhm 
                ellipticity = resellip 
        # NO good stars 
        else: 
            fwhm = 99.99 
            ellipticity = 99.99 
             
        # Print results
        if verbose:
            form = '%15s%10.3f%10.3f\n'
            maxlen = np.max(np.array([len(f) for f in files]))
            if maxlen > 15:
                form = '%'+str(maxlen)+'s%10.3f%10.3f\n'
            f.write(form % (files(f),fwhm,ellipticity))
             
            # Input into IMFWHMARR 
            allfwhm[f] = fwhm 
            allellip[f] = ellipticity 
             
            if len(outfile) > 0 : 
                printf,unit,format=form,files[f],fwhm,ellipticity 
     
    # Closing output file
    if outfile is not None:
        f.close()
     
    # Copy ALLFWHM to FWHM 
    fwhm = allfwhm 
    if len(fwhm) == 1: 
        fwhm = fwhm[0] 
    ellipticity = allellip 
    if len(ellipticity) == 1: 
        ellipticity = ellipticity[0] 
    # Output GTAB 
    if len(gtab2) > 0 : 
        gtab = gtab2 
 
    return fwhm, ellipticity, gtab, peaktab
