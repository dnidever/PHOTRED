#!/usr/bin/env python

import os
import time
import numpy as np
import glob as glob
import shutil
import warnings
import subprocess
import tempfile
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.convolution import convolve,Box2DKernel
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import SkyCoord
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import curve_fit
from dlnpyutils import utils as dln,coords
from . import utils,io,iraf

# Ignore these warnings
#warnings.filterwarnings(action="ignore", message=r'FITSFixedWarning:*')
warnings.simplefilter('ignore', category=AstropyWarning)

def ia_trim(xshift,yshift,xsize,ysize,trimsection,vignette):
    """ Helper function to deal with trim sections."""
     
   # Found this in the IMCENTROID source code: 
    # /iraf/iraf/pkg/images/immatch/src/listmatch/t_imctroid.x 
     
    ## IA_TRIM -- Compute the trim section. 
    # 
    #procedure ia_trim (cp) 
    # 
    #pointer cp                      #I center structure pointer 
    # 
    #real    xlo, xhi, ylo, yhi, xmin, ymin 
    #int     ixlo, ixhi, iylo, iyhi, ixlonew, ixhinew, iylonew, iyhinew, i 
    #int     vxlo, vxhi, vylo, vyhi          # vignetted versions 
    #bool    firsttime 
    # 
    #begin 
     
    nimages = len(xshift) 
     
    firsttime = 1 
    for i in range(nimages): 
         
        #firsttime = true 
        #do i = 1, NIMAGES(cp) { 
         
        #if (IS_INDEFR(XSHIFT(cp,i)) || IS_INDEFR(YSHIFT(cp,i))) 
        #    next 
         
        ## Compute limits. 
        #xlo = 1. + XSHIFT(cp,i) 
        #ylo = 1. + YSHIFT(cp,i) 
        #xhi = XSIZE(cp,i) + XSHIFT(cp,i) 
        #yhi = YSIZE(cp,i) + YSHIFT(cp,i) 
        xlo = 1.0 + xshift[i] 
        ylo = 1.0 + yshift[i] 
        xhi = xsize[i] + xshift[i] 
        yhi = ysize[i] + yshift[i] 
         
        #ixlonew = int (xlo) 
        #if (xlo > ixlonew)                  # round up 
        #    ixlonew = ixlonew + 1 
        ixlonew = np.ceil(xlo) 
         
        #ixhinew = int (xhi) 
        #if (xhi < ixhinew)                  # round down 
        #    ixhinew = ixhinew - 1 
        ixhinew = np.floor(xhi) 
         
        #iylonew = int (ylo)                 # round up 
        #if (ylo > iylonew) 
        #    iylonew = iylonew + 1 
        iylonew = np.ceil(ylo) 
         
        #iyhinew = int (yhi)                 # round down 
        #if (yhi < iyhinew) 
        #    iyhinew = iyhinew - 1 
        iyhinew = np.floor(yhi) 
         
        if (firsttime): 
            ixlo = ixlonew 
            ixhi = ixhinew 
            iylo = iylonew 
            iyhi = iyhinew 
             
            #xmin = XSIZE(cp,i) 
            #ymin = YSIZE(cp,i) 
            xmin = xsize[i] 
            ymin = ysize[i] 
             
            firsttime = 0 
        else: 
            #ixlo = max (ixlo, ixlonew) 
            #ixhi = min (ixhi, ixhinew) 
            #iylo = max (iylo, iylonew) 
            #iyhi = min (iyhi, iyhinew) 
            ixlo = np.max([ixlo,ixlonew]) 
            ixhi = np.min([ixhi,ixhinew]) 
            iylo = np.max([iylo,iylonew]) 
            iyhi = np.min([iyhi,iyhinew]) 
             
            #xmin = min (XSIZE(cp,i), xmin) 
            #ymin = min (YSIZE(cp,i), ymin) 
            xmin = np.min([xsize[i], xmin]) 
            ymin = np.min([ysize[i], ymin]) 
     
    ## Don't bother to complain. 
    #if (firsttime) 
    #    return 
     
    #call printf ("\n") 
     
    ## Vignetting is possible downstream since imshift and other tasks 
    ## preserve the size of the input image. 
     
    #vxlo = max (1, min (ixlo, int(xmin))) 
    #vxhi = max (1, min (ixhi, int(xmin))) 
    #vylo = max (1, min (iylo, int(ymin))) 
    #vyhi = max (1, min (iyhi, int(ymin))) 
    vxlo = np.max([1, np.min([ixlo, int(xmin)]) ]) 
    vxhi = np.max([1, np.min([ixhi, int(xmin)]) ]) 
    vylo = np.max([1, np.min([iylo, int(ymin)]) ]) 
    vyhi = np.max([1, np.min([iyhi, int(ymin)]) ]) 
     
    #if (vxlo != ixlo || vxhi != ixhi || vylo != iylo || vyhi != iyhi) { 
    #    call eprintf ("#Vignette_Section = [%d:%d,%d:%d]\n") 
    #        call pargi (vxlo) 
    #        call pargi (vxhi) 
    #        call pargi (vylo) 
    #        call pargi (vyhi) 
    #} 
    if (vxlo != ixlo or vxhi != ixhi or vylo != iylo or vyhi != iyhi): 
        vignette = [vxlo,vxhi,vylo,vyhi] 
        print('Vignette_Section = [',str(vxlo,2),':',str(vxhi,2),',',str(vylo,2),':',str(vyhi,2),']')
     
    ## Output the trim section. 
    #call printf ("#Trim_Section = [%d:%d,%d:%d]\n") 
    #    call pargi (ixlo) 
    #    call pargi (ixhi) 
    #    call pargi (iylo) 
    #    call pargi (iyhi) 
    trimsection = [ixlo,ixhi,iylo,iyhi] 
    print('Trim_Section = [',str(ixlo,2),':',str(ixhi,2),',',str(iylo,2),':',str(iyhi,2),']')
     
    # call flush (STDOUT) 

 
def calcweights(mag,err,fwhm,rdnoise,medsky):
    """
    Helper function to calculate weights
    """

    shape = mag.shape
    if mag.ndim==1:
        nfiles = len(mag)
        ngd = 1
    else:
        nfiles = shape[0]
        ngd = shape[1]
     
    # Make 2D arrays for fwhm, rdnoise and medsky
    fwhm2 = fwhm.reshape(-1,1) + np.ones(ngd,float).reshape(1,-1)
    rdnoise2 = rdnoise.reshape(-1,1) + np.ones(ngd,float).reshape(1,-1)
    medsky2 = medsky.reshape(-1,1) + np.ones(ngd,float).reshape(1,-1)
     
    # Computing intensity and weight for each star 
    #C Compute the intensity from the magnitude, easy enough 
    #            intensity=10**( (Mag(i,num(i))-25.0)/(-2.5)) 
    #C Now compute the S/N (called the weight here) 
    #            weight(i,n)=(intensity/(Pi*FWHM(i)**2)) / 
    #     &      (sqrt(RDNOISE(i)**2 
    #     &      + MEDSKY(i) + intensity/(Pi*FWHM(i)**2) ) ) 
    #C            print*, id(i,num(i)) 
    # magnitude zero-point: 1 star ADU == 25.0 mag 
    intensity = 10.0**( (mag-25.0)/(-2.5) ) 
    weight = (intensity/(np.pi*fwhm2**2.0))/ (np.sqrt(rdnoise2**2 + medsky2 + intensity/(np.pi*fwhm2**2) ) ) 
     
    # You get similar results if you use: weight = 1.0/err 
    # since magnitude errors are basically N/S 
     
    #C Now lets normalize the S/N (weight) for each star, take the maximimum S/N 
    #C and normalize so that is 1 S/N(i)/maximum S/N 
    #       do j=1,n-1 
    #          maxweight=-99999999.98 
    #          do i=1,FMAX 
    #           if(weight(i,j).gt.maxweight) maxweight=weight(i,j) 
    #          end do 
    #           weight(1:FMAX,j)=weight(1:FMAX,j)/maxweight 
    #C           print*,weight(1:FMAX,j),maxweight 
    #       end do 
     
    maxw = np.max(weight,axis=0)
    maxw2 = np.ones(nfiles,float).reshape(-1,1) + maxw.reshape(1,-1)
    nweight = weight/maxw2 
     
     
    #C Finally computed the actual weight, by summing up the normalized weights, 
    #C and dividing thru by the sum 
    #       do j=1,FMAX 
    #        avgweight(j)=sum(weight(j,1:n)) 
    #C/dble(n-1) 
    #C        print*, avgweight(j),n 
    #       end do 
    #        actweight(1:FMAX)=avgweight(1:FMAX)/sum(avgweight(1:FMAX)) 
    #C Print them out, so we can put them somewhere useful 
    #       do j=1,FMAX 
    #          print*,actweight(j) 
    #       end do 
    #       stop 
    #       end 
     
    avgweight = np.sum(nweight,axis=1) 
    actweight = avgweight/np.sum(avgweight) 
    # They sum to 1.0 
     
     
    # Compute the scaling for each frame.  scale=1 for the first frame 
    scales = np.zeros(nfiles,float)
    for i in range(nfiles): 
        ratio = intensity[i,:]/intensity[0,:]
        ratio_error = np.sqrt(err[i,:]**2 + err[0,:]**2)  # mag errors are ~fractional 
        med_ratio = np.median(ratio) 
        #wmeanerr,ratio,ratio_error,xmean,xsigma 
        scales[i] = med_ratio 
     
    return actweight,scales


def getweights_raw(tab):
    """ 
    This calculates weights, scales and sky values from the arrays from 
    the DAOPHOT .raw file.  Called by allframe_getweights.pro 
     
    Parameters
    ----------
    tab         Input structure (one element per file) giving 
                   MAG[Nstars], ERR[Nstars], FWHM, RDNOISE, MEDSKY 
     
    Returns
    -------
    out      Output structure similar to input structure but 
                    with WEIGHTS and SCALES added. 
     
    Example
    -------

    out = getweights_raw(tab)
     
    By D.Nidever   Jan 2019, broke out functionality from allframe_getweights.pro 
    """

     
    mag = tab['mag'] 
    err = tab['err']
    fwhm = tab['fwhm']
    rdnoise = tab['rdnoise']
    medsky = tab['medsky']

    if mag.ndim==2:
        nstars = mag.shape[1]
        nfiles = mag.shape[0]
    else: 
        nstars = len(mag)
        nfiles = 1

    # Getting the reference sources 
    totstars = np.sum(mag < 50,axis=0) 
    si = np.flip(np.argsort(totstars))# get the stars with the most detections 
    gdrefstars = si[0:np.minimum(nstars,50)] 
    nrefstars = len(gdrefstars) 
    # Getting the "good" frames 
    totframe = np.sum(mag[:,gdrefstars] < 50,axis=1) # of ref stars detected per frame 
    gdframe, = np.where(totframe == nrefstars)
    ngdframe = len(gdframe)
    bdframe, = np.where(totframe != nrefstars)
    nbdframe = len(bdframe)
    # No good frames, lower number of reference stars 
    if ngdframe == 0: 
        gdrefstars = si[0:np.minimum(nstars,30)] 
        nrefstars = len(gdrefstars) 
        totframe = np.sum(mag[:,gdrefstars] < 50,axis=1) 
        gdframe, = np.where(totframe == nrefstars)
        ngdframe = len(gdframe)
        bdframe, = np.where(totframe != nrefstars)
        nbdframe = len(bdframe)
                        
    # No good frames again, say the frame with the most detections is "good" 
    #   get weights relative to that one for the others 
    if ngdframe == 0: 
        # say the frame with the most detections is "good" and 
        #  the rest are bad, get weights relative to this one frame 
        totstarsframe = np.sum(mag < 50,axis=1) 
        gdframe = np.argmax(totstarsframe)
        gdframe = np.array([gdframe])
        ngdframe = 1
        bdframe = np.arange(nfiles)
        bdframe = np.delete(bdframe,gdframe)
        nbdframe = len(bdframe) 
        # Get stars that are good in this frames and in ALL of the others 
        gdrefstars, = np.where((mag[gdframe[0],:] < 50) & (totstars == nfiles))
        nrefstars = len(gdrefstars)
        #  lower threshold, at least half 
        if nrefstars == 0: 
            gdrefstars, = np.where((mag[gdframe[0],:] < 50) & (totstars > 0.5*nfiles))
            nrefstars = len(gdrefstars)
        #  just the good ones 
        if nrefstars == 0: 
            gdrefstars, = np.where(mag[gdframe[0],:] < 50)
            nrefstars = len(gdrefstars)
            si = np.flip(np.argsort(totstars[gdrefstars]))      # order by how many other frames they are detected in 
            gdrefstars = gdrefstars[si[0:np.minimum(50,nrefstars-1)]]   # only want 50 
            nrefstars = len(gdrefstars) 
     
    # Calculate the weights
    weights = np.zeros(nfiles,float)
    scales = np.zeros(nfiles,float)
    mag2 = mag[gdframe,:]
    err2 = err[gdframe,:]
    mag2 = mag2[:,gdrefstars]
    err2 = err2[:,gdrefstars]
    weights1,scales1 = calcweights(mag2,err2,fwhm[gdframe],rdnoise[gdframe],medsky[gdframe])
    weights[gdframe] = weights1 
    scales[gdframe] = scales1
     
    # If there are some "bad" frames calculate their weights 
    #  relative to some of the "good" ones 
    for i in range(nbdframe):      
        iframe = bdframe[i] 
         
        # First round of "good" stars 
        totgdstars = np.sum(mag[gdframe,:] < 50,axis=0)# stars in good frames 
        igdrefstars1, = np.where((mag[iframe,:] < 50) & (totgdstars >= 5))
        nigdrefstars1 = len(igdrefstars1)
        if nigdrefstars1 < 3: 
            igdrefstars1 , = np.where((mag[iframe,:] < 50) & (totgdstars >= 1))
            nigdrefstars1 = len(igdrefstars1)
        if nigdrefstars1 < 2: 
            continue
         
        totgdstars1 = np.sum(mag[gdframe,:][:,igdrefstars1] < 50,axis=0) 
        si1 = np.flip(np.argsort(totgdstars1))# get the best ones 
        igdrefstars = igdrefstars1[si1[0:(49<(nigdrefstars1-1))]] 
        nirefstars = len(igdrefstars)
         
        itotframe = np.sum(mag[gdframe,:][:,igdrefstars] < 50,axis=1) 
        igdframe1, = np.where(itotframe == nirefstars)
        nigdframe1 = len(igdframe1)
         
        # Whittle down to the best stars/frames 
        if nigdframe1 == 0: 
            igdframe = gdframe 
            igdrefstars = igdrefstars 
             
            # whittle down to best stars and frames
            count = 0
            while endflag==0:  
                # sum over frames 
                itot1 = np.sum(mag[igdframe,:][:,igdrefstars] < 50,axis=0) 
                si1 = np.argsort(itot1) 
                # remove worst 5% stars 
                if len(igdrefstars) > 3: 
                    bd1 = si1[0:int(np.round(len(igdrefstars)*0.05))]
                    #remove,bd1,igdrefstars
                    igdrefstars = np.delete(igdrefstars,bd1)
                              
                # sum over stars 
                itot2 = np.sum( mag[igdframe,:][:,igdrefstars] < 50,axis=1) 
                si2 = np.argsort(itot2) 
                # remove worst 5% frames 
                if len(igdframe) > 1: 
                    bd2 = si2[0:int(np.round(len(igdframe)*0.05))] 
                    #remove,bd2,igdframe
                    igdframe = np.delete(igdframe,bd2)
             
                # Testing if we should end 
                itotframe = np.sum( (mag[igdframe,:])[:,igdrefstars] < 50,axis=1) 
                igdframe1, = np.where(itotframe == len(igdrefstars))
                nigdframe1 = len(igdframe1)
                if nigdframe1 > 0:
                    endflag = 1
                if endflag == 1: 
                    igdframe = igdframe[igdframe1] 
                if count > 50:  # not converging, end it 
                    break
         
                count += 1
     
        else:
            igdframe = gdframe[igdframe1] 
     
        # Get weights relative to some "good" frames 
        indframes = np.append(igdframe,iframe)
        mag3 = mag[indframes,:][:,igdrefstars] 
        err3 = err[indframes,:][:,igdrefstars] 
        weights3,scales3 = calcweights(mag3,err3,fwhm[indframes],rdnoise[indframes],medsky[indframes])
     
        # Scale these relative to the original ones 
        weights3a = weights[igdframe]           # original 
        weights3b = weights3[0:nigdframe1]      # new 
        wtfrac = np.median(weights3a/weights3b) 
        scales3a = scales[igdframe]             # original
        scales3b = scales3[0:nigdframe1  ]      # new 
        sclfrac = np.median(scales3a/scales3b) 
        new_weights = weights3[nigdframe1] * wtfrac 
        new_scale = scales3[nigdframe1] * sclfrac 
        weights[iframe] = new_weights 
        scales[iframe] = new_scale 
     
        #print,iframe,new_weights,new_scale 
     
 
    # Fix images with bad weights likely due to no overlap 
    #  use the FLUX values to get a weight 
    bdweights, = np.where(weights <= 0.0)
    nbdweights = len(bdweights)
    gdweights, = np.where(weights > 0.0)
    ngdweights = len(gdweights)                          
    if nbdweights > 0: 
        # Use fluxrate10 to get weights and flux10 to get scales 
        # relative to "good" frame 
        if ngdweights > 0: 
            weights[bdweights] = tab['fluxrate10'][bdweights] * np.median([weights[gdweights]/tab['fluxrate10'][gdweights]]) 
            scales[bdweights] = tab['flux10'][bdweights] * np.median([scales[gdweights]/tab['flux10'][gdweights]]) 
        # all bad 
        else: 
            weights = tab['fluxrate10']
            scales = tab['flux10']
 
    # Normalize the weights 
    weights /= np.sum(weights) 
 
    # Rescale SCALES so the reference frames has SCALE=1.0 
    if scales[0] > 0.0: 
        scales /= scales[0] 
    else: 
        scales /= np.max(scales) 
 
    # Create the output structure
    out = tab.copy()
    if 'weight' not in out.colnames:
        out['weight'] = 0.0
    if 'scale' not in out.colnames:
        out['scale'] = 0.0
    out['weight'] = weights
    out['scale'] = scales

    return out

                              
def getweights(mchfile,imager=None,setup=None,logger=None,silent=False):
    """
    This calculates weights, scales and sky values from DAOPHOT 
    photometry files for image combination.  Used by combine and allframe.
 
    You need to have run daophot, allstar, daomatch and daomaster 
    already.  There need als, mch and raw files. 
 
    Parameters
    ----------
    mchfile : str
       The MCH filename 
    imager : dict
       Imager structure with basic information 
    setup : dict
       The setup information contained in the photred.setup file.
    logger : logging object, optional
       A logging object to use for printing.
    silent : bool, optional
       Don't print anything to the screen.   Default is False.
 
    Returns
    -------
    actweight   The array of weights. 
    scales      The array of scales. 
    medsky      The array of skys. 
    raw2       The RAW photometry structure. 
 
    Example
    -------

    weights,scales,sky = getweights('ccd1001.mch')
 
    By D.Nidever   February 2008, broken out into its own program 4/8/15 
    Translated to Python by D. Nidever, April 2022
    """

    # OUTPUTS: 
    #  actweight  The weight for each frame 
    #  scales     The scale for each frame 
    #  medsky     The sky value for each frame 
     
    tilesep = '+' 
    #tilesep = '.' 
    #btilesep = int(byte(tilesep)) 

    # MCH file not found 
    if os.path.exists(mchfile) == False:
        raise ValueError(mchfile+' NOT FOUND')

    if logger is None:
        logger = dln.basiclogger()
     
    # Load the MCH file
    files,trans,magoff = io.readmch(mchfile)
    nfiles = len(files) 
     
    #----------------------------------- 
    # Get the information that we need 
     
    # Load the opt files
    dt = [('name',(np.str,200)),('filter',(np.str,10)),('exptime',float),('fwhm',float),('rdnoise',float),
          ('mnsky',float),('medsky',float),('mag10',float),('flux10',float),('fluxrate10',float),
          ('weight',float),('scale',float)]
    info = np.zeros(nfiles,dtype=np.dtype(dt))
    info['name'] = files 
    for i in range(nfiles): 
        fdir = os.path.abspath(os.path.dirname(mchfile))
        base = os.path.splitext(os.path.basename(files[i]))[0]
        optfile = fdir+'/'+base+'.opt' 
        logfile1 = fdir+'/'+base+'.log' 

        opt = io.readopt(optfile)
        info['rdnoise'][i] = opt['RE']
        info['fwhm'][i] = opt['FW']
         
        loglines1 = dln.readlines(logfile1)
        out = dln.grep(loglines1,'Clipped')
        #out = subprocess.check_output(['grep','Clipped',logfile1],shell=False)
        #              Clipped mean and median =  187.442  187.215 
         
        # daophot.sh log files are clipped on Tortoise for some reason 
        #  Get mean/median sky level 
        if len(out) == 1 and out[0] == '': 
            print('Getting mean/median sky levels for ',base)
            if os.path.exists('daophot.opt')==False:
                shutil.copyfile(base+'.opt','daophot.opt')
            cmdlines = []
            cmdlines += ['#!/bin/sh']
            cmdlines += ['export image=${1}']
            cmdlines += ['daophot << END_DAOPHOT >> ${image}.find.log']
            cmdlines += ['OPTIONS']
            cmdlines += ['${image}.opt']
            cmdlines += [' ']
            cmdlines += ['ATTACH ${image}.fits']
            cmdlines += ['FIND']
            cmdlines += ['1,1']
            cmdlines += ['${image}.find.temp']
            cmdlines += ['y']
            cmdlines += ['EXIT']
            cmdlines += ['END_DAOPHOT']
            tid,tempscript = tempfile.mkstemp(prefix="dfind")  # absolute filename 
            dln.writelines(tempscript,cmdlines)
            os.chmod(tempscripts,0o755)
            for f in [base+'.find.temp',base+'.find.log']:
                if os.path.exists(f): os.remove(f)
            # Run DAOPHOT/FIND 
            out1 = subprocess.check_output(tempscript+' '+base)
             
            logfile2 = base+'.find.log'
            loglines2 = dln.readlines(logfile2)
            out = dln.grep(loglines2,'Clipped')
            #out = subprocess.check_output(['grep','Clipped',logfile2],shell=False)
            #              Clipped mean and median =  187.442  187.215 
             
            # Delete temporary files
            for f in [base+'.find.temp',tempscript]:
                if os.path.exists(f): os.remove(f)
         
        arr = out[0].split()
        info['mnsky'][i] = float(arr[5]) 
        info['medsky'][i] = float(arr[6]) 
         
        # Get exptime and filter 
        fitsfile = base+'.fits' 
        if os.path.exists(fitsfile)==False: 
            fitsfile+='.fz' 
        info['exptime'][i] = io.getexptime(fitsfile) 
        info['filter'][i] = io.getfilter(fitsfile,setup=setup) 

    info = Table(info)  # convert to astropy table
     
    # Only ONE file, return 1s 
    if nfiles == 1: 
        actweight = 1.0 
        scales = 1.0 
        medsky = info['medsky'][0]
        return actweight,scales,medsky
     
     
    #      program getweight 
    #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
    #C     Program computes weights according to signal to noise 
    #C     Input in aperture photometry at FWHM, the readnoise, 
    #C     and sky, and the program computes the S/N 
    #C        S/N = I/(Pi*FWHM**2) / sqrt(RN**2 + sky + I/(Pi*FWHM**2) ) 
    #C     This S/N is scaled for each set of stars, such that the maximum is 1 
    #C      then the "scaled" S/N are added up, to get a "total" S/N for the frame 
    #C      (dividing by n would give you the average, "scaled" S/N), lastly 
    #C      all of the "total" S/N for each frame are summed, and that sum is 
    #C      then divided into each "total" S/N to give a weight, good for use 
    #C      in imcombine. 
    #C      JCO -- 2000-ish 
    #C 
    #C      Use companion shell scripts to generate input files: 
    #C       getweights.sh 
    #C       rmweights.sh 
    #C      Altered fortran (I hope) to read inpfile in format created by 
    #C       these scripts. 
    #C      RLB -- 06/08/2007 
    #C 
    #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
     
    # Load the RAW file 
    fdir = os.path.abspath(os.path.dirname(mchfile))
    base = os.path.splitext(os.path.basename(mchfile))[0]
    rawfile = fdir+'/'+base+'.raw'
    raw,rawhead = io.readraw(rawfile)
     
    # Making mag and err arrays 
    nstars = len(raw) 
    tags = raw.colnames
    mag = np.zeros((nfiles,nstars),float)
    err = np.zeros((nfiles,nstars),float)
    for i in range(nfiles): 
        mag[i,:] = raw['mag'+str(i+1)]
        err[i,:] = raw['err'+str(i+1)] 
     
    # Calculate the magnitude and flux at the 10sigma magnitude 
    for i in range(nfiles): 
        gd, = np.where(mag[i,:] < 50)
        if len(gd) > 0:
            mag1 = mag[i,gd] 
            snr1 = 1.087/err[i,gd] 
            gdsnr , = np.where(np.abs(snr1-10.0) < 1)
            if len(gdsnr) < 5: 
                gdsnr, = np.where(np.abs(snr1-10.0) < 2) 
            if len(gdsnr) < 5: 
                si = np.argsort(abs(snr1-10.0)) 
                gdsnr = si[0:np.minimum(ngd,99)]
                ngdsnr = len(gdsnr) 
            mag10 = np.median([mag1[gdsnr]]) 
            info['mag10'][i] = mag10 
            info['flux10'][i] = 10.0**( (mag10-25.0)/(-2.5) )  # total counts 
            info['fluxrate10'][i] = info['flux10'][i] / info['exptime'][i]  # flux rate = counts / sec 
     
     
    # Using TILES 
    #-------------- 
    # We are using TILES and have multiple chips/amps 
    #   F1-00507800_39+T2.als, '+T' and two dots 
    if imager is not None:
        namps = imager['namps']
    else: 
        namps = 1
    if len(dln.grep(files,'\\'+tilesep+'T'))==nfiles and (np.sum(np.char.array(list(files[0]))==tilesep) >= 2) and (namps > 1): 
        usetiles = 1 
         
        # Number of unique exposures 
        expname = np.zeros(nfiles,(np.str,200))
        chip = np.zeros(nfiles,(np.str,10))
        for i in range(nfiles): 
            base1 = os.path.splitext(os.path.basename(files[i]))[0]  # F1-00507800_39+T2 
            field1,expchptile = base1.split('-')[0]                  # F1 and 00507800_39+T2 
            expchp = expchptile.split(tilesep)[0]                    # 00507800_39 
            expname[i] = expchp.split(imager.separator)[0] 
            chip[i] = expchp.split(imager.separator)[1] 
        # Unique exposures 
        uexpname,uiexp = np.unique(expname,return_index=True)
        nexp = len(uexpname) 
         
        # Combine all the catalogs for a given exposure 
        # Calculate the weights
        dt = [('mag',(float,nstars)),('err',(float,nstars)),('nfiles',int),('index',(int,nfiles)),
              ('exptime',float),('filter',(np.str,10)),('fwhm',float),('rdnoise',float),
              ('medsky',float),('mag10',float),('flux10',float),('fluxrate10',float)]
        expstr = np.zeros(nexp,dtype=np.dtype(dt))
        expstr = Table(expstr)
        expstr['mag10'] = 99.99
        expmag = np.zeros((nexp,nstars),float)
        experr = np.zeros((nexp,nstars),float)        
        for i in range(nexp): 
            ind, = np.where(expname == uexpname[i])
            nind = len(ind)
            expstr['nfiles'][i] = nind 
            expstr['index'][i][0:nind] = ind
            #expstr[i]['index'][0:nind-1] = ind             
            expstr['filter'][i] = info[ind[0]].filter 
            expstr['exptime'][i] = info[ind[0]].exptime 
            expstr['fwhm'][i] = np.median([info[ind].fwhm]) 
            expstr['rdnoise'][i] = np.median([info[ind].rdnoise]) 
            expstr['medsky'][i] = np.median([info[ind].medsky]) 
            expstr['mag10'][i] = np.median([info[ind].mag10]) 
            expstr['flux10'][i] = np.median([info[ind].flux10]) 
            expstr['fluxrate10'][i] = np.median([info[ind].fluxrate10]) 
            # Combine the photometry 
            mag1 = mag[ind,:] 
            err1 = err[ind,:] 
            # Multiple chips 
            #   they shouldn't overlap, so just use the mean/median 
            #   and that should pick up the detections 
            if nind>1: 
                bd, = np.where(mag1 > 50) 
                if len(bd)>0: 
                    mag1[bd] = np.inf
                    err1[bd] = np.inf
                expmag[i,:] = np.median(mag1,axis=0) 
                experr[i,:] = np.median(err1,axis=0) 
            else: 
                expmag[i,:] = mag1 
                experr[i,:] = err1 
        # Replace NANs with 99.9999 
        bdmag, = np.where(~np.isfinite(expmag))
        if len(bdmag)>0:
            expmag[bdmag] = 99.99 
            experr[bdmag] = 9.99 
        expstr['mag'] = expmag.T
        expstr['err'] = experr.T
        # Perform the weights and scales calculations 
        outexpstr = getweights_raw(expstr)
        # Give each chip the weight of its respective exposure 
        for i in range(nexp):
            info['weight'][expstr['index'][i][0:expstr['nfiles']]] = outexpstr['weight'][i]
            info['scake'][expstr['index'][i][0:expstr['nfiles']]] = outexpstr['scale'][i]
         
    # REGULAR Method 
    #--------------- 
    else:
        dt = [('mag',(float,nstars)),('err',(float,nstars)),
              ('exptime',float),('filter',(np.str,10)),('fwhm',float),('rdnoise',float),
              ('medsky',float),('flux10',float),('fluxrate10',float)]
        tab = np.zeros(nfiles,dtype=np.dtype(dt))
        tab = Table(tab)
        for n in tab.colnames:
            if n in info.colnames:
                tab[n] = info[n]
        tab['mag'] = mag
        tab['err'] = err
        # Perform the weights and scales calculations 
        outstr = getweights_raw(tab)
        info['weight'] = outstr['weight'] 
        info['scale'] = outstr['scale']
     
     
    # Print out the information 
    if silent==False:
        logger.info('        FILE         FILTER EXPTIME FWHM RDNOISE MEDSKY WEIGHT SCALE')
        for i in range(nfiles):
            #fmt = '(A-23,A4,F6.1,F6.2,F6.2,F8.1,F7.3,F7.3)'
            fmt = '%-23s%4s%6.1f%6.2f%6.2f%8.1f%7.3f%7.3f'
            data = (info['name'][i],info['filter'][i],info['exptime'][i],info['fwhm'][i],info['rdnoise'][i],
                    info['medsky'][i],info['weight'][i],info['scale'][i])
            logger.info(fmt % data)
     
    # Final output 
    actweight = info['weight']
    scales = info['scale']
    medsky = info['medsky']

    return actweight,scales,medsky

 
def adxyinterp(head,rr,dd,nstep=10,xyad=False):
    """
    Instead of transforming the entire large 2D RA/DEC 
    arrays to X/Y do a sparse grid and perform linear 
    interpolation to the full size. 
 
    Parameters
    ----------
    head : header
       The FITS header with the WCS. 
    rr : numpy array
       The 2D RA array. 
    dd  : numpy array
       The 2D DEC array. 
    nstep : int, optional
       The undersampling to use, default nstep=10. 
    xyad : bool, optional
       Perform X/Y->RA/DEC conversion instead. 
          In this case the meaning of the coordinate 
          arrays are flipped, i.e. rr/dd<->xx/yy 
 
    Returns
    -------
    xx : numpy array
       The 2D X array. 
    yy : numpy array
       The 2D Y array. 
 
    Example
    -------

    xx,yy = adxyinterp(head,ra,dec,nstep=10)
            
    By D. Nidever  Oct. 2016 
    Translated to Python by D. Nidever, April 2022
    """

    nx,ny = rr.shape
    nxs = (nx-1)//nstep + 1 
    nys = (ny-1)//nstep + 1 

    wcs = WCS(head)
     
    # Small dimensions, just transform the full grid 
    if nx <= nstep or ny <= nstep: 
        if xyad==False:
            xx,yy = wcs.world_to_pixel(SkyCoord(ra=rr,dec=dd,unit='deg'))
        else:
            coo = wcs.pixel_to_world(rr,dd)
            xx = coo.ra.deg
            yy = coo.dec.deg

    # Subsample the RA/DEC arrays 
    rrs = rr[0:nx:nstep,0:ny:nstep] 
    dds = dd[0:nx:nstep,0:ny:nstep] 
    if xyad==False:
        xxs,yys = wcs.world_to_pixel(SkyCoord(ra=rrs,dec=dds,unit='deg'))
    else:
        coos = wcs.pixel_to_world(rrs,dds)
        xxs = coos.ra.deg
        yys = coos.dec.deg
     
    # Start final arrays 
    xx = np.zeros((nx,ny),float)
    yy = np.zeros((nx,ny),float)
     
    # Use CONGRID to perform the linear interpolation 
    #   congrid normally does nx_out/nx_in scaling 
    #   if /minus_one set it does (nx_out-1)/(nx_in-1) 

    xx0,yy0 = np.arange(xxs.shape[0]),np.arange(xxs.shape[1])
    xx1 = np.arange((nxs-1)*nstep+1)/((nxs-1)*nstep)*(xxs.shape[0]-1)
    yy1 = np.arange((nys-1)*nstep+1)/((nys-1)*nstep)*(xxs.shape[1]-1)
    ixx = RectBivariateSpline(xx0,yy0,xxs,kx=1,ky=1)(xx1,yy1)
    iyy = RectBivariateSpline(xx0,yy0,yys,kx=1,ky=1)(xx1,yy1)
    xx[0:(nxs-1)*nstep+1,0:(nys-1)*nstep+1] = ixx 
    yy[0:(nxs-1)*nstep+1,0:(nys-1)*nstep+1] = iyy 
    #ixx = CONGRID(xxs,(nxs-1)*nstep+1,(nys-1)*nstep+1,/interp,/minus_one) 
    #iyy = CONGRID(yys,(nxs-1)*nstep+1,(nys-1)*nstep+1,/interp,/minus_one) 
     
    # Deal with right edge 
    if (nxs-1)*nstep+1 < nx: 
        # Make a short grid in X at the right edge 
        rrs_rt = rr[nx-nstep-1:nx:nstep,0:ny:nstep] 
        dds_rt = dd[nx-nstep-1:nx:nstep,0:ny:nstep] 
        if xyad==False:
            xxs_rt,yys_rt = wcs.world_to_pixel(SkyCoord(ra=rrs_rt,dec=dds_rt,unit='deg'))
        else: 
            coo_rt = wcs.pixel_to_world(rrs_rt,dds_rt)
            xxs_rt = coo_rt.ra.deg
            yys_rt = coo_rt.dec.deg

        xx0_rt,yy0_rt = np.arange(xxs_rt.shape[0]),np.arange(xxs_rt.shape[1])
        xx1_rt = np.arange(nstep+1)/nstep*(xxs_rt.shape[0]-1)
        yy1_rt = np.arange((nys-1)*nstep+1)/((nys-1)*nstep)*(xxs_rt.shape[1]-1)
        ixx_rt = RectBivariateSpline(xx0_rt,yy0_rt,xxs_rt,kx=1,ky=1)(xx1_rt,yy1_rt)
        iyy_rt = RectBivariateSpline(xx0_rt,yy0_rt,yys_rt,kx=1,ky=1)(xx1_rt,yy1_rt)
        xx[nx-nstep-1:nx,0:(nys-1)*nstep+1] = ixx_rt
        yy[nx-nstep-1:nx,0:(nys-1)*nstep+1] = iyy_rt 
        #ixx_rt = CONGRID(xxs_rt,nstep+1,nys*nstep+1,/interp,/minus_one) 
        #iyy_rt = CONGRID(yys_rt,nstep+1,nys*nstep+1,/interp,/minus_one) 
        #xx[nx-nstep-1:nx-1,0:(nys-1)*nstep] = ixx_rt
        #yy[nx-nstep-1:nx-1,0:(nys-1)*nstep] = iyy_rt

    # Deal with top edge 
    if (nys-1)*nstep+1 < ny: 
        # Make a short grid in Y at the top edge 
        rrs_tp = rr[0:nx:nstep,ny-nstep-1:ny:nstep] 
        dds_tp = dd[0:nx:nstep,ny-nstep-1:ny:nstep] 
        if xyad==False:
            xxs_tp, yys_tp = wcs.world_to_pixel(SkyCoord(ra=rrs_tp,dec=dds_tp,unit='deg'))
        else: 
            coo_tp = wcs.pixel_to_world(rrs_tp,dds_tp)
            xxs_tp = coo_tp.ra.deg
            yys_tp = coo_tp.dec.deg

        xx0_tp,yy0_tp = np.arange(xxs_tp.shape[0]),np.arange(xxs_tp.shape[1])
        xx1_tp = np.arange((nxs-1)*nstep+1)/((nxs-1)*nstep)*(xxs_tp.shape[0]-1)
        yy1_tp = np.arange(nstep+1)/nstep*(xxs_tp.shape[1]-1)
        ixx_tp = RectBivariateSpline(xx0_tp,yy0_tp,xxs_tp,kx=1,ky=1)(xx1_tp,yy1_tp)
        iyy_tp = RectBivariateSpline(xx0_tp,yy0_tp,yys_tp,kx=1,ky=1)(xx1_tp,yy1_tp)
        xx[0:(nxs-1)*nstep+1,ny-nstep-1:ny] = ixx_tp 
        yy[0:(nxs-1)*nstep+1,ny-nstep-1:ny] = iyy_tp 
        #ixx_tp = CONGRID(xxs_tp,nxs*nstep+1,nstep+1,/interp,/minus_one) 
        #iyy_tp = CONGRID(yys_tp,nxs*nstep+1,nstep+1,/interp,/minus_one) 
        #xx[0:(nxs-1)*nstep,ny-nstep-1:ny-1] = ixx_tp
        #yy[0:(nxs-1)*nstep,ny-nstep-1:ny-1] = iyy_tp

    # Deal with top/right corner 
    if (nxs-1)*nstep+1 < nx and (nys-1)*nstep+1 < ny: 
        # Make a short grid in X and Y at the top-right corner 
        rrs_tr = rr[nx-nstep-1:nx:nstep,ny-nstep-1:ny:nstep] 
        dds_tr = dd[nx-nstep-1:nx:nstep,ny-nstep-1:ny:nstep] 
        if xyad==False:
            xxs_tr, yys_tr = wcs.world_to_pixel(SkyCoord(ra=rrs_tr,dec=dds_tr,unit='deg'))
        else: 
            coo_tr = wcs.pixel_to_world(rrs_tr,dds_tr)
            xxs_tr = coo_tr.ra.deg
            yys_tr = coo_tr.dec.deg

        xx0_tr,yy0_tr = np.arange(xxs_tr.shape[0]),np.arange(xxs_tr.shape[1])
        xx1_tr = np.arange(nstep+1)/nstep*(xxs_tr.shape[0]-1)
        yy1_tr = np.arange(nstep+1)/nstep*(xxs_tr.shape[1]-1)
        ixx_tr = RectBivariateSpline(xx0_tr,yy0_tr,xxs_tr,kx=1,ky=1)(yy1_tr,xx1_tr)
        iyy_tr = RectBivariateSpline(xx0_tr,yy0_tr,yys_tr,kx=1,ky=1)(yy1_tr,xx1_tr)
        xx[nx-nstep-1:nx,ny-nstep-1:ny] = ixx_tr
        yy[nx-nstep-1:nx,ny-nstep-1:ny] = iyy_tr
        #ixx_tr = CONGRID(xxs_tr,nstep+1,nstep+1,/interp,/minus_one) 
        #iyy_tr = CONGRID(yys_tr,nstep+1,nstep+1,/interp,/minus_one) 
        #xx[nx-nstep-1:nx-1,ny-nstep-1:ny-1] = ixx_tr
        #yy[nx-nstep-1:nx-1,ny-nstep-1:ny-1] = iyy_tr

    return xx,yy

def fiximages(input,satlevel=50000.0,nofix=False):
    """
    This programs fixes bad pixels in images. 
     
    Parameters
    ----------
    input : list or st
       The filenames 
    satlevel float, optional
       The saturation level to use if NOT 
         in the FITS header.  Default is 50000.
    nofix : boolean, optional
       Don't fix the images. 
     
    Returns
    -------
    The files are overwritten with the fixed images. 
     
    Example
    -------

    fiximages('*.fits',satlevel=5e4)
     
    By D.Nidever   August 2008 
    Translated to Python by D. Nidever, April 2022
    """


    # Load the input 
    files = dln.loadinput(input)
    nfiles = len(files)

    # No files to process 
    if nfiles == 0:
        raise ValueError('No files to process')
     
    # Loop through the files 
    for i in range(nfiles): 
        ifile = files[i] 
        ibase = os.path.splitext(os.path.basename(ifile))[0]
        idir = os.path.dirname(ifile) 
         
        # Does the file exist 
        if os.path.exists(ifile)==False:
            print(ifile,' NOT FOUND')
            continue
         
        # Read the FITS file 
        im,head = fits.getdata(ifile)
         
        # Add mask name to header 
        maskname = idir+'/'+ibase+'.mask.fits' 
        head['BPM'] = maskname,' BPM mask file' 
         
        # Do we have bad pixels? 
        saturate = head.get('SATURATE')
        if saturate is None:
            saturate = satlevel 
        bd, = np.where(im > saturate)
         
        # Fix bad pixels 
        if len(bd) > 0:
            gd, = np.where(im < saturate)             
            # Background value 
            backg = np.median(im[gd]) 
             
            # Set all bad pixels to the background value 
            first = im 
            first[bd] = backg 
             
            # Do one smoothing 
            smlen = 3
            #sm = smooth(first,smlen,/nan,/edge_truncate)
            kernel = Box2DKernel(smlen,mode='center')
            sm = convolve(first,kernel)
            
            # Use the convolved image for the bad pixels 
            newim = im 
            newim[bd] = sm[bd] 
             
            # Print 
            print(ifile,'  '+str(nbd)+' bad pixels fixed')
             
            # Write the output 
            outname = ifile
            if nofix==False:
                fits.PrimaryHDU(newim,head).writeto(outname,overwrite=True)
             
            # Make the bad pixel mask image, 0-bad, 1-good 
            mask = np.ones(im.shape,int)
            mask[bd] = 0
            fits.PrimaryHDU(mask,head).writeto(maskname,overwrite=True)
            
        # NO bad pixels 
        else: 
            print(ifile,' has no bad pixels')
            # Make the bad pixel mask image, 0-bad, 1-good 
            mask = np.ones(im.shape,int)
            fits.PrimaryHDU(mask,head).writeto(maskname,overwrite=True)
         
                      
def combine(filename,tile=None,setup=None,scriptsdir=None,logger=None,irafdir=None,
            satlevel=6e4,nocmbimscale=False,fake=False,usecmn=False,imager=None,
            verbose=True):
    """
    This combines/stacks images for ALLFRAME. 
 
    You need to have run daophot, allstar, daomatch and daomaster 
    already.  There needs to be a fits, opt, als.opt, ap and als 
    file for each file in the MCH file.  There also needs to be an 
    associated RAW file for the MCH file. 
 
    Parameters
    ----------
    filename : str
       The MCH filename
    tile : dict
       Tile information on the sky.
    setup : dict
       The information contained in the photred.setup file.
    nocmbimscale : boolean, optional
       Don't scale the images when combining them.  Not 
        recommended, but the old way of doing it.  Bright 
        stars can be missed this way.  Default is False.
    scriptsdir : str
       The directory that contains all of the necessary scripts. 
    irafdir : str
       The IRAF home directory. 
    satlevel : float, optional
       Saturation level.  Default is 6e4.
    logger : logging object
       A logging object to use for printing to the screen.
    fake : boolean, optiona
       Run for artificial star tests.  Default is False.
    usecmn : boolean, optional
       Use the individual cmn.lst files to construct a 
         cmn.lst file for the combined image. 
    imager : dict
       Imager structure with basic information. 
 
    Returns
    -------
    The combined image and mask: 
    FILEBASE_comb.fits       combined image 
    FILEBASE_comb.bpm.fits   bad pixel mask 
    FILEBASE_comb.mask.fits  SExtractor weight map (very similar to bpm) 
    FILEBASE_comb.mch        The transformations of the individual frames 
                               to the combined frame. 
    maskdatalevel : float
       The "bad" data level above the highest "good" value 
         in the combined image. 
    fileinfo : table
       Information on all of the files used for the 
         stacked iamge. 
 
    Example
    -------

    combine('ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/')
 
 
    By D.Nidever   February 2008 
    Automation of steps and scripts by J.Ostheimer and Rachael Beaton 
    Major upgrade of the resampling and combination code, Oct 2016  
    Translated to Python by D. Nidever,  April 2022
    """ 

    # Set up logging to screen and logfile
    if logger is None:
        logger = dln.basiclogger()
     
    # Getting scripts directory and iraf directory 
    if setup is not None:
        scriptsdir = setup['scriptsdir']
        irafdir = setup['irafdir']
     
    # No irafdir 
    if irafdir is None:
        raise ValueError('IRAFDIR NOT INPUT')
    # No irafdir 
    if scriptsdir is None:
        raise ValueError('SCRIPTSDIR NOT INPUT')    
     
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if iraf.check(irafdir=irafdir)==False:
        raise ValueError('IRAF TEST FAILED.  EXITING')

    # Check if the scripts exist in the current directory 
    scripts = ['getpsf.sh','photo.opt','apcor.opt','lstfilter.py','goodpsf.pro','allframe.opt',
               'default.sex','default.param','default.nnw','default.conv'] 
    nscripts = len(scripts) 
    # Loop through the scripts 
    for i in range(nscripts): 
        exists = os.path.exists(scriptsdir+'/'+scripts[i])
        if exists:
            size = os.path.getsize(scriptsdir+'/'+scripts[i])
        else:
            size = 0
        curexists = os.path.exists(scripts[i])
        if curexists:
            cursize = os.path.getsize(scripts[i])
        else:
            cursize = 0
         
        # No file 
        if exists == False or size == 0: 
            raise ValueError(scriptsdir+'/'+scripts[i],' NOT FOUND or EMPTY')
         
        # Check if the two files are the same size, if not copy it 
        if size != cursize:
            if os.path.exists(scripts[i]): os.remove(scripts[i])
            shutil.copyfile(scriptsdir+'/'+scripts[i],scripts[i])
     
    # FILENAME 
    mchfile = os.path.basename(filename) 
    mchdir = os.path.dirname(filename) 
    mchbase = os.path.splitext(filename)[0]
     
    # CD to the directory
    curdir = os.path.abspath(os.getcwd())
    os.chdir(mchdir)



    
    # Check that the mch, als, and opt files exist 
    if os.path.exists(mchfile)==False:
        raise ValueError(mchfile+' NOT FOUND')
    # Checking RAW file 
    if os.path.exists(mchbase+'.raw')==False:
        raise ValueError(mchbase+'.raw NOT FOUND')
     
     
    ############################################ 
    # CHECK NECESSARY FILES 
     
    # Load the MCH file
    files,trans,magoff = io.readmch(mchfile)
    nfiles = len(files) 
     
    # FAKE, check that we have all the files that we need 
    if fake: 
        # weights, scale, zero, comb_psf, _shift.mch 
        chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_comb.mch']
        exists = [os.path.exists(f) for f in chkfiles]
        nbdfiles = np.sum(np.array(exists)==False)
        if nbdfiles>0: 
            raise ValueError('FAKE.  Some necessary files not found. '+' '.join(chkfiles[bdfiles]))
     
    # Final files
    base = [os.path.splitext(os.path.basename(f))[0] for f in files]
    fitsfiles = [b+'.fits' for b in base]
    outfiles = [b+'.shft.fits' for b in base]
    outmaskfiles = [b+'.mask.shft.fits' for b in base]
          
    # Gather information on all of the files 
    # photred_gatherfileinfo.pro can do most of this 
    logger.info('Gathering file information')
    ntrans = len(trans[0,:])
    fileinfo1 = {'fitsfile':'','catfile':'','nx':0,'ny':0,'trans':np.zeros(ntrans,float),'magoff':np.zeros(2),
                 'head':None,'vertices_ra':np.zeros(4,float),'vertices_dec':np.zeros(4,float),
                 'pixscale':0.0,'saturate':0.0,'background':0.0,'comb_zero':0.0,'comb_scale':0.0,
                 'comb_weights':0.0,'resampfile':'','resampmask':'','resamptrans':np.zeros(ntrans,float),
                 'resamptransrms':0.0}
    fileinfo = nfiles*[None]
    for i in range(nfiles):
        fileinfo[i] = fileinfo1.copy()
        fileinfo[i]['fitsfile'] = fitsfiles[i]
        fileinfo[i]['catfile'] = files[i]
        fileinfo[i]['trans'] = trans[i,:]
        fileinfo[i]['magoff'] = magoff[i,:]
        fileinfo[i]['resampfile'] = outfiles[i] 
        fileinfo[i]['resampmask'] = outmaskfiles[i]
        im1,head1 = io.readfile(fileinfo[i]['fitsfile'])
        fileinfo[i]['head'] = head1
        fileinfo[i]['nx'] = head1['NAXIS1'] 
        fileinfo[i]['ny'] = head1['NAXIS2']
        wcs1 = WCS(head1)
        vcoo = wcs1.pixel_to_world([0,fileinfo[i]['nx']-1,fileinfo[i]['nx']-1,0],[0,0,fileinfo[i]['ny']-1,fileinfo[i]['ny']-1])
        vra = vcoo.ra.deg
        vdec = vcoo.dec.deg
        fileinfo[i]['vertices_ra'] = vra 
        fileinfo[i]['vertices_dec'] = vdec
        pixscale = io.getpixscale('',head=head1)
        fileinfo[i]['pixscale'] = pixscale 
        saturate = head1.get('SATURATE')
        if saturate is None:
            saturate = 50000 
        fileinfo[i]['saturate'] = saturate 
        gdpix = (im1 < saturate)
        background = np.median(im1[gdpix]) 
        fileinfo[i]['background'] = background 
     
     
    ################################################## 
    # Create default reference frame if TILE not input 
    if tile is None:
        tile = {'type':'WCS'} 
    if tile['type'] == 'WCS' and len(tile)==1: 
        logger.info('Creating TILE projection')
        # The default projection is a tangent plane centered 
        # at halfway between the ra/dec min/max of all of 
        # the images.  The mean pixel scale is used. 
        #  near RA=0 line 
        if dln.valrange(fileinfo['vertices_ra']) > 180: 
            vertices_ra = fileinfo['vertices_ra']
            over, = np.where(vertices_ra > 180)
            if len(over)>0: 
                vertices_ra[over]-=360 
            rar = dln.minmax(vertices_ra) 
            cenra = np.mean(rar) 
        else: 
            rar = dln.minmax(fileinfo['vertices_ra']) 
            cenra = np.mean(rar) 
        decr = dln.minmax(fileinfo['vertices_dec']) 
        cendec = np.mean(decr) 
        pixscale = np.mean(fileinfo['pixscale']) 
        # Set up the tangent plane projection 
        step = pixscale/3600.0
        delta_dec = dln.valrange(decr) 
        delta_ra = dln.valrange(rar)*np.cos(np.deg2rad(cendec))
        nx = np.ceil(delta_ra*1.01/step) 
        ny = np.ceil(delta_dec*1.01/step) 
        xref = nx/2 
        yref = ny/2
        tilehead = fits.Header()
        tilehead['NAXIS1'] = nx 
        tilehead['CDELT1'] = step 
        tilehead['CRPIX1'] = xref+1
        tilehead['CRVAL1'] = cenra 
        tilehead['CTYPE1'] = 'RA---TAN' 
        tilehead['NAXIS2'] = ny 
        tilehead['CDELT2'] = step 
        tilehead['CRPIX2'] = yref+1 
        tilehead['CRVAL2'] = cendec 
        tilehead['CTYPE2'] = 'DEC--TAN' 
        tilewcs = WCS(tilehead)
        tileast['equinox'] = 2000 
         
        logger.info('RA range = ['+str(rar[0])+','+str(rar[1])+'] deg')
        logger.info('DEC range = ['+str(decr[0])+','+str(decr[1])+'] deg')
        logger.info('Central RA = '+str(cenra))
        logger.info('Central DEC = '+str(cendec)) 
        logger.info('NX = '+str(nx))
        logger.info('NY = '+str(ny))
         
        # Create the TILE structure 
        tile = {'type':'WCS','naxis':np.array([nx,ny]),'cdelt':np.array([step,step]),'crpix':np.array([xref+1,yref+1]),
                'crval':np.array([cenra,cendec]),'ctype':['RA--TAN','DEC--TAN'],
                'head':tilehead,'wcs':tilewcs,'xrange':[0,nx-1],'yrange':[0,ny-1],'nx':nx,'ny':ny} 
     
    # Check that the TILE is valid 
    if utils.validtile(tile)==False:
        raise ValueError('Tile ERROR')

    # Add HEAD/AST to TILE if needed 
    if 'head' not in tile.keys():
        tilehead = fits.Header()
        tilehead['NAXIS1'] = tile['naxis'][0] 
        tilehead['CDELT1'] = tile['cdelt'][0] 
        tilehead['CRPIX1'] = tile['crpix'][0] 
        tilehead['CRVAL1'] = tile['crval'][0] 
        tilehead['CTYPE1'] = tile['ctype'][0] 
        tilehead['NAXIS2'] = tile['naxis'][1] 
        tilehead['CDELT2'] = tile['cdelt'][1] 
        tilehead['CRPIX2'] = tile['crpix'][1] 
        tilehead['CRVAL2'] = tile['crval'][1] 
        tilehead['CTYPE2'] = tile['ctype'][1] 
        if 'cd' in tile: 
            tilehead['CD1_1'] = tile['cd'][0,0] 
            tilehead['CD1_2'] = tile['cd'][0,1] 
            tilehead['CD2_1'] = tile['cd'][1,0] 
            tilehead['CD2_2'] = tile['cd'][1,1] 
        tile['head'] = tilehead
        tile['wcs'] = WCS(tilehead)
    if 'wcs' not in tile: 
        tile['wcs'] = WCS(tilehead)
    # Add XRANGE/YRANGE 
    if 'xrange' not in tile: 
        tile['xrange'] = [0,tile['naxis'][0]-1]
        tile['nx'] = tile['naxis'][0]
    if 'yrange' not in tile: 
        tile['yrange'] = [0,tile['naxis'][1]-1]
        tile['ny'] = tile['naxis'][1]
     
     
    ############################################ 
    # STEP 1: IMALIGN PREP 
     
    logger.info('------------------------')
    logger.info('STEP 1: GETTING WEIGHTS')
    logger.info('------------------------')
     
    #----------------------------------- 
    # Computs Weights 
    if fake==False:
        weights,scales,sky = getweights(mchfile,imager=imager,logger=logger,setup=setup)
        invscales = 1.0/scales 
        bdscale, = np.where((scales < 1e-5) | (invscales > 900))
        if len(bdscale) > 0: 
            scales[bdscale] = 1.0 
            invscales[bdscale] = 1.0 
            weights[bdscale] = 0.0 
        weightfile = mchbase+'.weights'
        dln.writelines(weightfile,weights)
        #WRITECOL,weightfile,weights,fmt='(F10.6)' 
        scalefile = mchbase+'.scale'
        dln.writelines(scalefile,invscales)
        #WRITECOL,scalefile,invscales,fmt='(F10.6)'# want to scale it UP 
        zerofile = mchbase+'.zero'
        dln.writelines(zerofile,-sky)
        #WRITECOL,zerofile,-sky,fmt='(F12.4)'# want to remove the background, set to 1st frame 
         
    # FAKE, use existing ones 
    else: 
        weightfile = mchbase+'.weights' 
        scalefile = mchbase+'.scale' 
        zerofile = mchbase+'.zero'
        weights = np.array(dln.readlines(weightfile)).astype(float)
        invscales = np.array(dln.readlines(scalefile)).astype(float)        
        #READCOL,weightfile,weights,format='F',/silent 
        #READCOL,scalefile,invscales,format='F',/silent 
        scales = 1.0/invscales
        sky = np.array(dln.readlines(zerofile)).astype(float)        
        #READCOL,zerofile,sky,format='F',/silent 
        sky = -sky 
    # Stuff the information into the FILEINFO structure 
    for i in range(nfiles):
        fileinfo[i]['comb_weights'] = weights[i]
        fileinfo[i]['comb_scale'] = invscales[i]
        fileinfo[i]['comb_zero'] = -sky[i] 
     
     
    ############################################ 
    # STEP 2: Resample the images 
     
    logger.info('---------------------------')
    logger.info('STEP 2: Resample the images')
    logger.info('---------------------------')

    # --- WCS projection on the sky ---     
    if tile['type']=='WCS':
        # Convert x/y and ra/dec grid for entire tile image 
        #  that we'll be using, ~15sec 
        #xb = (lindgen(tile.nx)+tile.xrange[0])#replicate(1,tile.ny) 
        #yb = replicate(1,tile.nx)#(lindgen(tile.ny)+tile.yrange[0]) 
        #HEAD_XYAD,tile.head,xb,yb,rab,decb,/deg
        xb = (np.arange(tile['nx'])+tile['xrange'][0]).reshape(-1,1) + np.zeros(tile['ny'],float).reshape(1,-1)
        yb = np.zeros(tile['nx'],float).reshape(-1,1) + (np.arange(tile['ny'])+tile['yrange'][0]).reshape(1,-1)
        wcs = WCS(tile['head'])
        bcoo = wcs.pixel_to_world(xb,yb)
        rab = bcoo.ra.deg
        decb = bcoo.dec.deg

        # Loop through the files 
        for i in range(nfiles):
            im1,head1 = io.readfile(fileinfo[i]['fitsfile'])
            # Make the mask 
            mask = np.ones((fileinfo[i]['ny'],fileinfo[i]['nx']),bool)
            gdpix = (im1 < fileinfo[i]['saturate'])
            bdpix = (im1 >= fileinfo[i]['saturate'])
            # False-bad, True-good 
            if np.sum(bdpix) > 0: 
                mask[bdpix] = False 
                im1[bdpix] = fileinfo[i]['background']
             
            # Get X/Y range for this image in the final coordinate system
            vx,vy = wcs.world_to_pixel(SkyCoord(ra=fileinfo[i]['vertices_ra'],dec=fileinfo[i]['vertices_dec'],unit='deg'))
            #HEAD_ADXY,tile.head,fileinfo[i].vertices_ra,fileinfo[i].vertices_dec,vx,vy,/deg 
            xout = [np.maximum(np.floor(np.min(vx))-2, tile['xrange'][0]),
                    np.minimum(np.ceil(np.max(vx))+2, tile['xrange'][1]-1)+1] 
            xout = np.array(xout)
            xoutrel = xout-tile['xrange'][0]  # relative to xrange[0] 
            xoutrel = xoutrel.astype(int)
            nxout = int(xout[1]-xout[0])
            yout = [np.maximum(np.floor(np.min(vy))-2, tile['yrange'][0]),
                    np.minimum(np.ceil(np.max(vy))+2, tile['yrange'][1]-1)+1] 
            yout = np.array(yout)
            youtrel = yout-tile['yrange'][0]  # relative to yrange[0] 
            youtrel = youtrel.astype(int)
            nyout = int(yout[1]-yout[0])
            rr = rab[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] 
            dd = decb[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] 
            #ALLFRAME_ADXYINTERP,head1,rr,dd,xx,yy,nstep=10
            xx,yy = adxyinterp(head1,rr,dd,nstep=10)

            # The x/y position to bilinear need to be in the original system, ~1sec 
            rim = np.zeros(xx.shape,float)+fileinfo[i]['background']
            rmask = np.zeros(xx.shape,bool)
            good = ((xx>=0) & (xx<=im1.shape[0]-1) & (yy>=0) & (yy<=im1.shape[1]-1))
            if np.sum(good)>0:
                rim[good] = RectBivariateSpline(np.arange(im1.shape[0]),np.arange(im1.shape[1]),im1,kx=1,ky=1).ev(xx[good],yy[good])
                rmask[good] = RectBivariateSpline(np.arange(im1.shape[0]),np.arange(im1.shape[1]),mask,kx=1,ky=1).ev(xx[good],yy[good])
            #rim = BILINEAR(im1,xx,yy,missing=fileinfo[i].background) 
            #rmask = BILINEAR(mask,xx,yy,missing=0) 
             
            # Contruct final image
            fim = np.zeros((tile['nx'],tile['ny']),float)+fileinfo[i]['saturate']
            fim[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] = rim 
            #fmask = bytarr(tile.nx,tile.ny)
            fmask = np.zeros((tile['nx'],tile['ny']),bool)
            fmask[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] = rmask 
             
            # Contruct the final header 
            fhead = head1.copy()
            # Delete any previous WCS keywords
            for n in ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2','CTYPE1','CTYPE2','CD1_1','CD1_2','CD2_1','CD2_2']:
                if n in fhead:
                    del fhead[n]
            cards = [f[0] for f in fhead.cards]
            pvnames = dln.grep(cards,'PV[0-9]_+[0-9]')
            #pvind, = np.where(stregex(strmid(fhead,0,5),'PV[0-9]_+[0-9]',/boolean) == 1,npvind)
            for p in pvnames:
                del fhead[p]
            # Add the new WCS 
            whead = wcs.to_header()
            fhead.extend(whead)
            fhead['NAXIS1'] = tile['nx']
            fhead['NAXIS2'] = tile['ny']
            fhead['BPM'] = fileinfo[i]['resampmask']
            logger.info('%d  %s [%d:%d,%d:%d] ' % (i+1,fileinfo[i]['resampfile'],xout[0],xout[1],yout[0],yout[1]))
            fits.PrimaryHDU(fim,fhead).writeto(fileinfo[i]['resampfile'],overwrite=True)
            #MWRFITS,fim,fileinfo[i].resampfile,fhead,/create 
            mhead = fhead.copy()
            mhead['BITPIX'] = 8 
            del mhead['BPM']
            fits.PrimaryHDU(fmask.astype(np.uint8),mhead).writeto(fileinfo[i]['resampmask'],overwrite=True)
            #MWRFITS,fmask,fileinfo[i].resampmask,mhead,/create 
            # this takes about ~37-50 sec for a 2kx4k image. 
            #  now it takes ~5 sec for a 2kx4k image 
 
     
    # --- Pixel based --- 
    elif tile['type']=='PIXEL':
    
        # Expand the images to sizes that will allow all of the shifts
        hd1 = io.readfile(fitsfiles[0],header=True)
        nx = hd1['NAXIS1'] 
        ny = hd1['NAXIS2']
        #  left, down, right, up 
        pix_expand = [ np.abs(np.minimum(np.floor(np.min(xshift)),0)), np.abs(np.minimum(np.floor(np.min(yshift)),0)),
                       np.maximum(np.ceil(np.max(xshift)),0), np.maximum(np.ceil(np.max(yshift)),0) ]
        outmaskfiles = [os.path.dirname(f)+'/'+os.path.splitext(os.path.basename(f))[0]+'.mask.shft.fits' for f in tempfits]
        logger.info('Expanding images by [',strjoin(str(pix_expand,2),','),'] pixels')
        files2,trans2,magoff2 = io.readmch(shiftmch+'.mch')
        nxf = nx+pix_expand[0]+pix_expand[2] 
        nyf = ny+pix_expand[1]+pix_expand[3]
        xx1 = np.arange(nx).reshape(-1,1) + np.zeros(ny,float).rshape(-1,1) 
        yy1 = np.zeros(nx,float).reshape(-1,1) + np.arange(ny).reshape(-1,1)
        #xx1 = lindgen(nx)#replicate(1,ny) 
        #yy1 = replicate(1,nx)#lindgen(ny) 
        for i in range(nfiles): 
            # Image
            tim,thead = io.readfile(tempfits[i])
            background = np.median(tim) 
            out = trans_coo(xx1[:],yy1[:],trans[i,:]) 
            xx2 = np.copy(xx1)*0. 
            yy2 = np.copy(yy1)*0. 
            xx2[:] = out[0,:] + pix_expand[0]  # shift to expanded grid 
            yy2[:] = out[1,:] + pix_expand[1]  # shift to expanded grid 
            triangulate,xx2,yy2,tr,b          # triangulate 
            xout = np.arange(nxf) 
            yout = np.arange(nyf) 
            tim2 = TRIGRID(xx2,yy2,tim, tr, XOUT = xout, YOUT = yout, missing=background) 
            thead2 = thead 
            thead2['NAXIS1'] = nxf 
            thead2['NAXIS2'] = nyf 
            thead2['CRPIX1'] = thead2['CRPIX1']+pix_expand[0] 
            thead2['CRPIX2'] = thead2['CRPIX2']+pix_expand[1]
            fits.PrimaryHDU(tim2,thead2).writeto(outfiles[i],overwrite=True)
            # Mask 
            mfile = os.path.basename(tempfits[i],'.fits')+'.mask.fits'
            mim,mhead = io.readfile(mfile)
            mim2 = TRIGRID(xx2,yy2,mim, tr, XOUT = xout, YOUT = yout, missing=background) 
            mhead2 = mhead 
            mhead2['NAXIS1'] = nxf 
            mhead2['NAXIS2'] = nyf 
            mhead2['CRPIX1'] = mhead2['CRPIX1']+pix_expand[0] 
            mhead2['CRPIX2'] = mhead2['CRPIX2']+pix_expand[1]
            fits.PrimaryHDU(mim2,mhead2).writeto(outmaskfiles[i],overwrite=True)
 
    else:
        raise ValueError(str(tile['type'])+' not implemented yet')
 
 
 
    # Creating new MCH file for the combined file 
    if fake==False:
        print('Deriving new transformation equations for the resampled coordinate system')
        mchfinal = []
        for i in range(nfiles): 
            # Convert X/Y of this system into the combined reference frame 
            #  The pixel values are 1-indexed like DAOPHOT uses. 
            ngridbin = 50 
            nxgrid = fileinfo[i]['nx'] // ngridbin 
            nygrid = fileinfo[i]['ny'] // ngridbin
            xgrid = (np.arange(nxgrid)*ngridbin+1).reshape(-1,1) + np.zeros(nygrid,float).reshape(1,-1)
            ygrid = np.zeros(nxgrid,float).reshape(-1,1) + (np.arange(nygrid)*ngridbin+1).reshape(1,-1)
            #xgrid = (lindgen(nxgrid)*ngridbin+1)#replicate(1,nygrid) 
            #ygrid = replicate(1,nxgrid)#(lindgen(nygrid)*ngridbin+1) 
            #HEAD_XYAD,(*fileinfo[i].head),xgrid-1,ygrid-1,ragrid,decgrid,/deg
            #HEAD_ADXY,tile.head,ragrid,decgrid,refxgrid,refygrid,/deg
            wcs1 = WCS(fileinfo[i]['head'])
            coogrid = wcs1.pixel_to_world(xgrid-1,ygrid-1)
            refxgrid,refygrid = wcs1.world_to_pixel(coogrid)
            refxgrid += 1    # convert 0-indexed to 1-indexed 
            refygrid += 1 
              
            # Now fit the transformation 
            xdiff = refxgrid-xgrid 
            ydiff = refygrid-ygrid 
            xmed = np.median([xdiff]) 
            ymed = np.median([ydiff]) 
            # Fit rotation with linear fits if enough points 
            slp1 = dln.mediqrslope(ygrid,xdiff)  # fit rotation term
            slp1rms = dln.mad(xdiff-ygrid*slp1)
            slp2 = dln.mediqrslope(xgrid,ydiff)
            slp2rms = dln.mad(ydiff-xgrid*slp2) 
            #coef1 = robust_poly_fitq(ygrid,xdiff,1)   # fit rotation term 
            #coef1b = dln.poly_fit(ygrid,xdiff,1,measure_errors=xdiff*0+0.1,sigma=coef1err) #,/bootstrap) 
            #coef2 = robust_poly_fitq(xgrid,ydiff,1)   # fit rotation term 
            #coef2b = dln.poly_fit(xgrid,ydiff,1,measure_errors=ydiff*0+0.1,sigma=coef2err) #,/bootstrap) 
            ##theta = mean([-coef1[1],coef2[1]]) 
            #theta,thetaerr = dln.wtmean([-coef1[1],coef2[1]],[coef1err[1],coef2err[1]])
            theta = dln.wtmean(np.array([-slp1,slp2]),np.array([slp1rms,slp2rms]))

            # [xoff, yoff, cos(th), sin(th), -sin(th), cos(th)] 
            #trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0] 
            trans = np.array([xmed, ymed, 1.0-theta**2, theta, -theta, 1.0-theta**2])
            # Adjust Xoff, Yoff with this transformation 
            #xyout = trans_coo(xgrid,ygrid,trans) 
            xout,yout = utils.trans_coo([xgrid,ygrid],*trans)
            trans[0] += np.median(refxgrid-xout) 
            trans[1] += np.median(refygrid-yout) 
     
            # Fit full six parameters if there are enough stars 
            xdata = [[refxgrid,refygrid],[xgrid,ygrid]]
            null = np.zeros(refxgrid.size,float)
            fpar,cov = curve_fit(utils.trans_coo_dev,xdata,null,p0=trans)
            trans = fpar
            transerr = np.sqrt(np.diag(cov))
            diff = utils.trans_coo_dev(xdata,*fpar)
            rms = np.sqrt(np.mean(diff**2.))

            #fa = {'x1':refxgrid.flatten(),'y1':refygrid.flatten(),'x2':xgrid.flatten(),'y2':ygrid.flatten()} 
            #fpar = MPFIT('trans_coo_dev',initpar,functargs=fa, perror=perror, niter=iter, status=status,
            #             bestnorm=chisq,:f=dof, autoderivative=1, /quiet) 
            #trans = fpar 
     
            #diff = trans_coo_dev(fpar,x1=refxgrid,y1=refygrid,x2=xgrid,y2=ygrid) 
            #rms = np.sqrt(np.mean(diff**2)) 
            fileinfo[i]['resamptrans'] = trans 
            fileinfo[i]['resamptransrms'] = rms 
     
            # The output is: 
            # filename, xshift, yshift, 4 trans, mag offset, magoff sigma 
            #format = '(A2,A-30,A1,2A10,4A12,F9.3,F8.4)' 
            # In daomaster.f the translations are 10 digits with at most 4 
            # decimal places (with a leading space), the transformation 
            # coefficients are 12 digits with at most 9 decimal places. 
            # Need a leading space to separate the numbers. 
            #strans = ' '+[str(string(trans[0:1],format='(F30.4)'),2),                 str(string(trans[2:5],format='(F30.9)'),2)] 
            #newline = STRING("'",fileinfo[i]['catfile'],"'", strans, fileinfo[i]['magoff'][0], rms, format=format) 
            #mchfinal += [newline]
            
            strans = ['%30.4f' % trans[0], '%30.4f' % trans[1], '%30.9f' % trans[2], '%30.9f' % trans[3],
                      '%30.9f' % trans[4], '%30.9f' % trans[5]]
            strans = [' '+s.strip() for s in strans]
            fmt = '%2s%-30s%1s%10.10s%10.10s%12.12s%12.12s%12.12s%12.12s%9.3f%8.4s'
            data = "'",fileinfo[i]['catfile'],"'",*strans, 0.0, rms,
            newline = fmt % data
            mchfinal += [newline]

        # Printing the transformation 
        #logger.info('(A-20,2A10,4A12,F9.3,F8.4)' % (fileinfo[i]['catfile'],strans,fileinfo[i]['magoff'][0],rms))
        logger.info('%-22s%10.10s%10.10s%12.12s%12.12s%12.12s%12.12s%9.3f%8.4f' % (fileinfo[i]['catfile'],*strans,fileinfo[i]['magoff'][0],rms))
        # Write to the new MCH file 
        combmch = mchbase+'_comb.mch'
        dln.writelines(combmch,mchfinal)
 
    # FAKE, use existing one 
    else: 
        combmch = mchbase+'_comb.mch' 
        # don't need to load the information 
 
 
    ############################################ 
    # STEP 5: COMBINE IMAGES 
    logger.info('-------------------')
    logger.info('STEP 5: IMCOMBINE')
    logger.info('-------------------') 
 
    # The imcombine input file 
    resampfile = mchbase+'.resamp' 
    dln.writelines(resampfile,[f['resampfile'] for f in fileinfo])
    #WRITELINE,resampfile,fileinfo.resampfile 
 
    # SCALE the images for combining 
    #------------------------------- 
    if nocmbimscale==False:
 
        # Put BPM mask names in file headers 
        #  these will be used by IMCOMBINE 
        #for i=0,nfiles-1 do begin 
        #  head = headfits(outfiles[i]) 
        #  sxaddpar,head,'BPM',outmaskfiles[i] 
        #  modfits,outfiles[i],0,head 
        #endfor 
        
        # Combine the frames WITH scaling/offset/masking, for the bright stars 
        #logger.info('Creating SCALED image' 
        combfile = mchbase+'_comb.fits'
        if os.path.exists(combfile): os.remove(combfile)
        if os.path.exists(mchbase+'_comb.bpm.pl'): os.remove(mchbase+'_comb.bpm.pl')
        iraf.imcombine('@'+resampfile,combfile,combine='average',reject='avsigclip',
                       weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',
                       irafdir=irafdir,scale='@'+scalefile,zero='@'+zerofile,
                       masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm')
  
        # Convert BPM mask from PL to FITS
        if os.path.exists(mchbase+'_comb.bpm.fits'): os.remove(mchbase+'_comb.bpm.fits')
        lines = []
        curdir = os.getcwd()
        lines += ['print("")']  # first line will be ignored 
        lines += ['cd '+curdir]
        lines += ['imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits']
        lines += ['logout']
        tid,tmpfile = tempfile.mkstemp(prefix="tiraf")  # absolute filename
        dln.writelines(tmpfile,lines)
        out = iraf.run(tmpfile,irafdir,verbose=verbose)
 
        # Delete temporary scripts and PL file
        for f in [tmpfile,mchbase+'_comb.bpm.pl']:
            if os.path.exists(f): os.remove(f)
 
 
        # Fix the rdnoise and background/sky level and saturate 
        #  the bad pixels for DAOPHOT 
        #------------------------------------------------------ 
 
        # 10/02/12 
        # THE IMAGES ARE (1) ZERO-SUBTRACTED, (2) SCALED, AND (3) WEIGHT AVERAGED 
        # The algorithm is: 
        # 1.) add zero-level correction.  im = im+zero 
        # 2.) scale the images.  im = im*scale 
        # 3.) take weighted average.  combim=total(weight*im) 
        #      there is also clipping that takes place during the averaging 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2)) 
        # A gain that changes from frame to frame could be problematic, 
        # but this shouldn't happen since it's the same chip from the same night. 
        
        # IMCOMBINE wants rdnoise in electrons and gain in electrons/DN. 
        # DAOPHOT expects rdnoise in DN.  That's why mkopt converts 
        #  it with the gain.  So we are fine.  The header should have 
        #  rdnoise in ELECTRONS. 
        
        # page 63-65 of daophot2.pdf shows how you need to modify rdnoise/gain 
        # when averaging/summing frames. in observing/mosaic/. 
        
        # Load the IMCOMBINE output combined file and BPM
        combim,combhead = io.readfile(combfile)
        badmask,maskhead = io.readfile(mchbase+'_comb.bpm.fits')   # 0-good, 1-bad 
 
        # Fix the gain 
        # For N averaged frames gain(N)=N*gain(1) 
        # Leave the gain as is!  We are scaling everything to the reference 
        # and using its gain.  It's nearly impossible to figure out the real 
        # gain since we are scaling the images and then taking a weighted 
        # average with outlier rejection.  Find a gain that properly 
        # describes/follows Poisson errors for the final combined image is 
        # difficult/impossible.  But that's okay.  This is just for source 
        # detection and DAOPHOT FIND just cares about the noise in the 
        # background.  We just need to ensure that the sky and rdnoise 
        # are correct. 
        
        # Fix the rdnoise 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2)) 
        rdnoisearr = np.zeros(nfiles,float)
        for i in range(nfiles): 
            rdnoisearr1,key = io.getrdnoise(base[i]+'.fits')
            rdnoisearr[i] = rdnoisearr1
        #  the "scales" array here is actually 1/scales used by IMCOMBINE. 
        rdnoise = np.sqrt(np.sum((weights*rdnoisearr/scales)**2)) 
        rdnoise = np.maximum(rdnoise, 0.01)  # must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
        dummy,rdnoisekey = io.getrdnoise(combfile)  # get keyword 
        combhead[rdnoisekey] = rdnoise 
 
        # Fix the sky 
        # DAOPHOT FIND computes the random error per pixel in ADU as 
        # noise = sqrt( sky level/gain + rdnoise^2) 
        # So it assumes that the noise in the background is sqrt(sky/gain) 
        # in ADU.  We need to set the sky level so this is correct. 
        # The final noise should be 
        # final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
        # So the final sky level should be 
        # final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
        gain,gainkey = io.getgain(combfile)
        comb_sky = np.sum((weights*np.sqrt(np.maximum(sky,0)/gain)/scales)**2)*gain 
        # the "scales" array here is actually 1/scale 
        combim += comb_sky.astype(float)   # keep it float 
        
 
        # set the maximum to a "reasonable" level 
        # Rescale the image and increase the gain 
        if np.max(combim) > 50000: 
            rescale = 50000./np.max(combim) 
            combim = combim*rescale
            combhead[gainkey] = gain/rescale
            # rdnoise does NOT get modified since it's in electrons 
            # we just need to modify the gain which takes you from ADUs to electrons 
 
        maskdatalevel = np.max(combim) + 10000    # set "bad" data level above the highest "good" value 
        combim2 = combim*(1-badmask) + maskdatalevel*badmask   # set bad pixels to maskdatalevel
        combhead['SATURATE'] = maskdatalevel
        fits.PrimaryHDU(combim2,combhead).writeto(combfile,overwrite=True) # fits_write can create an empty PDU 
        
        # Create the weight map for Sextractor using the BPM output by IMCOMBINE 
        #  bad only if bad in ALL images 
        weightmap = -2.0*(badmask == 1).astype(float) + 1.0 
        combweightfile = mchbase+'_comb.mask.fits' 
        fits.PrimaryHDU(weightmap,whead).writeto(combweightfile,overwrite=True)
        
    # NO SCALING of the images for combining 
    #--------------------------------------- 
    else:         
        combfile = mchbase+'_comb.fits'
        if os.path.exists(combfile): os.remove(combfile)
        iraf.imcombine('@'+resampfile,combfile,combine='average',reject='avsigclip',
                       weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',irafdir=irafdir)
 
        # Fix the rdnoise and background/sky level and saturate 
        #  the bad pixels for DAOPHOT 
        #------------------------------------------------------ 
 
        # See the explanations for all these steps above!! 
 
        # Load the IMCOMBINE output combined file and BPM
        combim,combhead = io.readfile(combfile)
        badmask,maskhead = io.readfile(mchbase+'_comb.bpm.fits')  # 0-good, 1-bad 
 
        # Fix the rdnoise 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2)) 
        rdnoisearr = np.zeros(nfiles,float) 
        for i in range(nfiles): 
            rdnoisearr[i] = utils.getrdnoise(base[i]+'.fits') 
        rdnoise = np.sqrt(np.sum((weights*rdnoisearr)**2)) 
        rdnoise = np.maximum(rdnoise, 0.01)  # must be >=0.01 or it will be 0.00 in the opt file and daophot will crash
        dummy,rdnoisekey = utils.getrdnoise(combfile) # get keyword
        combhead[rdnoisekey] = rdnoise
 
        # Fix the sky 
        # So the final sky level should be 
        # final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
        gain,gainkey = utils.getgain(combfile)
        comb_sky = np.sum((weights*np.sqrt(sky/gain))**2)*gain 
        # the "scales" array here is actually 1/scale 
        combim += comb_sky 
 
 
        # set the maximum to a "reasonable" level 
        # Rescale the image and increase the gain 
        if np.max(combim) > 50000: 
            rescale = 50000./np.max(combim) 
            combim = combim*rescale
            combhead[gainkey] = gain/rescale
            # rdnoise does NOT get modified since it's in electrons 
            # we just need to modify the gain which takes you from ADUs to electrons 
 
 
        # Making Sextractor "weight" map file 
        #------------------------------------ 
        # masks have 0-bad, 1-good. 
        # anything with less than 1.0 is considered bad 
        # weight map, -1 is bad, +1 is good 
        # "bpm" is the SUM of the bad pixel masks 
        # consider a pixel bad that is bad in ANY image 
        weightmap = -2.0*(bpm < nfiles).astype(float) + 1. 
        combweightfile = mchbase+'_comb.mask.fits'
        fits.PrimaryHDU(weightmap,whead).writeto(combweightfile,overwrite=True)
 
        #--------------------------------------------- 
        # SATURATE BAD pixels in the COMBINED IMAGE 
        # DAOPHOT needs to have the bad pixels "saturated", 
        # SExtractor will know which pixels are bad from the "weight" map. 
        # 
        # We could skip the fiximage.pro step but we still need the 
        # individual bpm masks and setting the bad pixels to the background 
        # probably helps in the IMALIGN/IMCOMBINE steps. 
        logger.info('')
        logger.info('"Saturating" bad pixels in the COMBINED image')
        logger.info('')
 
 
        badmask = float(weightmap < 0.5) 
        maskdatalevel = np.max(combim) + 10000  # set "bad" data level above the highest "good" value 
        combim2 = combim*(1.0-badmask) + maskdatalevel*badmask  # set bad pixels to 100,000 
        combhead['SATURATE'] = maskdatalevel 
        fits.PrimaryHDU(combim2,combhead).writeto(combfile,overwrite=True)
 
    # Add TILETYPE to the combined image
    combhead = io.readfile(combfile,header=True)
    combhead['AFTILTYP'] = tile['type']
    tempim = fits.getdata(combfile)
    fits.PrimaryHDU(tempim,combhead).writeto(combfile,overwrite=True)
    #MODFITS,combfile,0,combhead 
 
    # Delete the resampled images
    for f in fileinfo:
        if os.path.exists(f['resampfile']): os.remove(f['resampfile'])
        if os.path.exists(f['resampmask']): os.remove(f['resampmask'])        
 
    # Make the common source file for the combined image 
    #--------------------------------------------------- 
    if usecmn: 
        print('Combining COMMON SOURCE files for the combined image.')
        # Loop through the files and convert to coordinates to the comined file 
        allcmn = []
        for i in range(nfiles): 
            cmnfile1 = os.path.basename(fileinfo[i].fitsfile,'.fits')+'.cmn.lst' 
            if os.path.exists(cmnfile1) == 1:
                cmn1 = Table.read(cmnfile1,format='ascii')
                #cmn1 = IMPORTASCII(cmnfile1,fieldnames=['id','x','y','mag','err','sky','skysig','sharp','round','round2'],
                #                   skipline=3,/silent) 
                cmnlines1 = dln.readlines(cmnfile1)
                coohead1 = cmnlines1[0:1] 
                ncmn1 = len(cmn1) 
                # Get coordinates on the resampled/combined image grid 
                newx,newy = utils.trans_coo([cmn1['x'],cmn1['y']],*fileinfo[i]['resamptrans']) 
                cmn1['x'] = newx 
                cmn1['y'] = newy 
                if len(allcmn) == 0: 
                    allcmn = cmn1 
                else: 
                    # Remove any duplicates 
                    ind1,ind2,dist = coords.xmatch(allcmn['x'],allcmn['y'],cmn1['x'],cmn1['y'],2.0)
                    if len(ind1)>0:
                        if nmatch < ncmn1: 
                            cmn1 = np.delete(cmn1,ind2)
                            #remove,ind2,cmn1 
                        else: 
                            undefine,cmn1 
                    if len(cmn1) > 0 : 
                        allcmn += [cmn1]
        if len(allcmn) > 0: 
            dln.writecol(mchbase+'_comb.cmn.lst',allcmn['id'],allcmn['x'],allcmn['y'],allcmn['mag'],
                         allcmn['err'],allcmn['sky'],allcmn['skysig'],allcmn['sharp'],allcmn['round'],
                         allcmn['round2'],fmt='(I7,2F9.2,3F9.3,F9.2,3F9.3)')
            dln.writelines(mchbase+'_comb.cmn.lst',[coohead1,''],prepend=True)  #  prepend the COO header  
            #WRITELINE,mchbase+'_comb.cmn.lst',[coohead1,''],/prepend  # prepend the COO header 
        else:
            logger.info('No common file to combine')
     
    # Don't use common file 
    else: 
        # Make sure it doesn't exist otherwise it will be used
        if os.path.exist(mchbase+'_comb.cmn.lst'): os.remove(mchbase+'_comb.cmn.lst')

    return maskdatalevel,fileinfo

 
def combine_orig(filename,scriptsdir=None,logfile=None,
                 irafdir=None,satlevel=6e4,nocmbimscale=False,
                 trimcomb=False,fake=False,usecmn=False):
    """
    This is the old version of the code that combines images ALLFRAME. 
 
    You need to have run daophot, allstar, daomatch and daomaster 
    already.  There needs to be a fits, opt, als.opt, ap and als 
    file for each file in the MCH file.  There also needs to be an 
    associated RAW file for the MCH file. 
 
    Parameters
    ----------
    filename : str
       The MCH filename 
    nocmbimscale : boolean, optional
       Don't scale the images when combining them.  Not 
         recommended, but the old way of doing it.  Bright 
         stars can be missed this way. 
    trimcomb : boolean, optional
       Trim the combined images to the overlapping region. 
         This used to be the default, but now the default 
         is to keep the entire original region. 
    scriptsdir : str
       The directory that contains all of the necessary scripts. 
    irafdir : str
       The IRAF home directory. 
    logfile : str
       A logfile to print to output to. 
    fake : boolean, optional
       Run for artificial star tests. 
    usecmn : boolean, optional
       Use the cmn.lst file of the reference image for the 
         combined image. 
 
    Returns
    -------
    The combined image and mask: 
    FILEBASE.comb.fits       combined image 
    FILEBASE.comb.bpm.fits   bad pixel mask 
    FILEBASE.comb.mask.fits  SExtractor weight map (very similar to bpm) 
    FILEBASE.comb.mch        The transformations of the individual frames 
                               to the combined frame. 
    maskdatalevel : float
       The "bad" data level above the highest "good" value 
         in the combined image. 
    fileinfo : catalog
       Information about all of the files.
 
    Example
    -------

    maskdatalevel,fileinfo = combine_orig('ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/',finditer=2)
 
 
    By D.Nidever   February 2008 
    Automation of steps and scripts by J.Ostheimer and Rachael Beaton 
    Translated to Python by D. Nidever, April 2022
    """
 
              
    global setup 

     
    # Logfile 
    if logfile is not None:
        logf = logfile 
    else: 
        logf = -1 
     
    # Getting scripts directory and iraf directory
    if setup is not None:
        scriptsdir = setup['SCRIPTSDIR']
        irafdir =setup['IRAFDIR'] 
     
     
    # No irafdir 
    if irafdir is None:
        raise ValueError('IRAFDIR NOT INPUT')
    # No scriptsdir 
    if scriptsdir is None:
        raise ValueError('SCRIPTSDIR NOT INPUT')    
     
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if iraf.check(irafdir=irafdir)==False:
        raise ValueError('IRAF TEST FAILED')

    # Check if the scripts exist in the current directory 
    scripts = ['getpsf.sh','photo.opt','apcor.opt','lstfilter.py','goodpsf.pro','allframe.opt',
               'default.sex','default.param','default.nnw','default.conv'] 
    nscripts = len(scripts) 
    # Loop through the scripts 
    for i in range(nscripts): 
        exists = os.path.exists(scriptsdir+'/'+scripts[i])
        if exits:
            size = os.path.size(scriptsdir+'/'+scripts[i])
        else:
            size = 0
        curexists = os.path.exists(scripts[i])
        if curexists:
            cursize = os.path.size(scripts[i])
        else:
            cursize = 0
         
        # No file 
        if exists == False or size == 0: 
            raise ValueError(scriptsdir+'/'+scripts[i],' NOT FOUND or EMPTY')
         
        # Check if the two files are the same size, if not copy it 
        if isize != cursize:
            if os.path.exists(scripts[i]): os.remove(scripts[i])
            shutil.copyfile(scriptsdir+'/'+scripts[i],scripts[i])

     
    logger.info('Combining images in '+filename) 
     
    # FILENAME 
    mchfile = os.path.basename(filename) 
    mchdir = os.path.dirname(filename) 
    mchbase = os.path.splitext(os.path.basename(filename))[0]
     
     
    # CD to the directory
    curdir = os.getcwd()
    os.chdir(mchdir)
     
    # Check that the mch, als, and opt files exist 
    if os.path.exists(mchfile)==False:
        raise ValueError(mchfile+' NOT FOUND')
     
    # Checking RAW file 
    if os.path.exists(mchbase+'.raw')==False:
        raise ValueError(mchbase+'.raw NOT FOUND')
     
     
    ############################################ 
    # CHECK NECESSARY FILES 
     
    # Load the MCH file
    files,trans,magoffset = io.readmch(mchfile)
     
    # Check that the fits, als, opt, and psf files exist 
    nfiles = len(files) 
    for i in range(nfiles): 
        base = os.path.splitext(os.path.basename(files[i]))[0]
        for e in ['.fits','.opt','.als.opt','.ap','.als']:
            if os.path.exists(base+e)==False:
                raise ValueError(base+e+' NOT FOUND')
     
    # FAKE, check that we have all the files that we need 
    if fake: 
        # weights, scale, zero, comb_psf, _shift.mch
        chkfiles = [mchbase+e for e in ['.weights','.scale','.zero','_comb.psf','_shift.mch']]
        exists = [os.path.exists(f) for f in chkfiles]
        nbdfiles = np.sum(np.array(exists)==False)
        if nbdfiles > 0:
            raise ValueError('FAKE.  Some necessary files not found. '+' '.join(chkfiles[bdfiles]))
     
     
    ############################################ 
    # STEP 1: IMALIGN PREP 
     
    logger.info('')
    logger.info('Step A: Getting Weights')
    logger.info('-----------------------')
     
    #----------------------------------- 
    # Computs Weights 
    if fake==False:
        weights,scales,sky = getweights(mchfile)
        #ALLFRAME_GETWEIGHTS,mchfile,weights,scales,sky#,raw2=raw2 
        invscales = 1.0/scales 
        bdscale, = np.where((scales < 1e-5) | (invscales > 900))
        if len(bdscale) > 0: 
            scales[bdscale] = 1.0 
            invscales[bdscale] = 1.0 
            weights[bdscale] = 0.0 
        weightfile = mchbase+'.weights'
        weights = dln.fread(weightfile,'F10.6')
        #WRITECOL,weightfile,weights,fmt='(F10.6)' 
        scalefile = mchbase+'.scale'
        invscales = dln.fread(scalefile,'F10.6')
        #WRITECOL,scalefile,invscales,fmt='(F10.6)'# want to scale it UP 
        zerofile = mchbase+'.zero'
        sky = -dln.fread(zerofile,'F12.4')
        #WRITECOL,zerofile,-sky,fmt='(F12.4)'# want to remove the background, set to 1st frame 

    # FAKE, use existing ones 
    else: 
        weightfile = mchbase+'.weights' 
        scalefile = mchbase+'.scale' 
        zerofile = mchbase+'.zero'
        weights = dln.fread(weightfile,'F10.6')
        invscales = dln.fread(scalefile,'F10.6')
        #READCOL,weightfile,weights,format='F',/silent 
        #READCOL,scalefile,invscales,format='F',/silent 
        scales = 1.0/invscales
        sky = -dln.fread(zerofile,'F12.4')
        #READCOL,zerofile,sky,format='F',/silent 
     
     
    #--------------------------------------- 
    # Get X/Y translations using DAOMASTER 
    #  NO ROTATION ONLY SHIFTS 
    #  Need these shifts for IMSHIFT 
    shiftmch = mchbase+'_shift'
    if fake==False:
        logger.info('Measuring X/Y shifts')
        if os.path.exists(mchbase+'.mch'): os.remove(mchbase+'.mch')
        shutil.copyfile(mchbase+'.mch',shiftmch+'.mch')
        # Make the DAOMASTER script 
        cmdlines = []
        cmdlines += ['#!/bin/csh']
        cmdlines += ['set input=${1}']
        cmdlines += ['daomaster <<DONE']
        cmdlines += ['${input}.mch']
        cmdlines += ['1,1,1']
        cmdlines += ['99.']
        cmdlines += ['2']
        cmdlines += ['10'] 
        cmdlines += ['5']
        cmdlines += ['4']
        cmdlines += ['3']
        cmdlines += ['2']
        cmdlines += ['1']
        cmdlines += ['0']
        cmdlines += ['n']
        cmdlines += ['n']
        cmdlines += ['n']
        cmdlines += ['n']
        cmdlines += ['y']
        cmdlines += ['']
        cmdlines += ['']
        cmdlines += ['n'] 
        cmdlines += ['n']
        cmdlines += ['n']
        cmdlines += ['DONE']
        tid,tempscript = tempfile.mkstemp(prefix="daomaster")  # absolute filename
        dln.writelines(tempscript,cmdlines)
        os.chmod(tempscript,0o755)
        # Run DAOMASTER 
        cmd2 = tempscript+' '+shiftmch 
        out2 = subprocess.check_output(cmd2)
        # Remove temporary DAOMASTER script
        if os.path.exists(tempscript): os.remove(tempscript)
    files2,trans2,magoffset2 = io.readmch(shiftmch+'.mch')
     
    xshift = trans2[:,0] 
    yshift = trans2[:,1]
    xyshifts = [[xshift],[yshift]] 
    logger.info('Image shifts')
    for i in range(nfiles): 
        logger.info(files[i],xshift[i],yshift[i])
     
     
    #----------------------------------- 
    # Create imalign prep files 
    # This is done by preimalign_k.sh 
    # Need an input list of fits files 
    # Need an output list of fits files 
    # Shift file 
    base = [os.path.splitext(os.path.basename(f))[0] for f in files]
    fitsfiles = [b+'.fits' for b in base]
    outfiles = [b+'.shft.fits' for b in base]
    infile = mchbase+'.inlist'
    outfile = mchbase+'.outlist'
    #WRITELINE,infile,fitsfiles   ; this is done below now with the temp files 
    dln.writelines(outfile,outfiles)  
    # Remove outfiles
    for f in outfiles:
        if os.path.exists(f): os.remove(f)
    # shift list 
    shiftfile = mchbase+'.shift'
    dln.writecol(shiftfile,xshift,yshift,fmt='(2F15.4)')
     
     
    # Make temporary files for bad pixel fixing and combining 
    #  FIXIMAGES doesn't work properly on the shifted images 
    #  because the interpolation can bring the bad pixel values down 
    #  below the saturation threshold, and we don't want to touch 
    #  the original images. 
    tempfits = base+'.temp.fits'
    for i in range(len(fitsfiles)):
        if os.path.exists(tempfits[i]): os.remove(tempfits[i])
        shutil.copyfile(fitsfiles[i],tempfits[i])
    dln.writelines(infile,tempfits)
     
     
    ############################################ 
    # STEP B: FIX BAD PIXELS 
    logger.info('')
    logger.info('Step B: Fixing bad pixels')
    logger.info('-------------------------')
    fiximages('@'+infile,satlevel=satlevel)  #6e4 
    # This also makes the FILE.mask.fits files for each image 
     
    # Find the maximum saturation level 
    satlevelarr = fltarr(nfiles) 
    for i in range(nfiles): 
        #head = headfits(base[i]+'.fits')
        im,head = io.readfile(base[i]+'.fits')
        saturate = head.get('SATURATE')
        if saturate is None:
            saturate = np.max(im)-1000. 
        satlevelarr[i] = saturate 
    maxsatlevel = np.max(satlevelarr) 
     
     
    ############################################ 
    # STEP C: IMALIGN 
    #  This figures out the X/Y-shifts between the images 
    #  and creates the shifted images (".shft.fits") 
    logger.info('')
    logger.info('Step C: IMALIGN')
    logger.info('---------------')
    reffile = mchbase+'.fits' 
     
    # IMALIGN basically is a script that runs: 
    #  IMCENTROID - to compute the shifts 
    #  IMSHIFT - shifts the images 
    #  IMCOPY - trims the images 
     
    # First, shift the images 
    logger.info('Shifting the images')
    imshift('@'+infile,'@'+outfile,shifts_file=shiftfile,interp_type='linear',
            boundary_type='constant',constant=0,irafdir=irafdir)
     
    # Trim the images 
    if trimcomb: 
        # Calculate the trim section 
        hd = io.readfile(reffile,header=True) 
        xsize = lonarr(nfiles)+hd['NAXIS1']
        ysize = lonarr(nfiles)+hd['NAXIS2']
        trimsection = iatrim(xshift,yshift,xsize,ysize)
        xoff = trimsection[0]-1 
        yoff = trimsection[2]-1 
         
        # Trim the shifted images 
        logger.info('Trimming the shifted images')
        xstart = trimsection[0]-1 
        xstop = trimsection[1]-1 
        ystart = trimsection[2]-1 
        ystop = trimsection[3]-1 
         
        for i in range(nfiles):
            im,head = io.readfile(outfiles[i])
            newim = im[xstart:xstop,ystart:ystop]
            fits.PrimaryHDU(newim,head).writeto(outfiles[i],overwrite=True)
         
    # Don't trim the images 
    else: 
        hd = io.readfile(reffile,header=True) 
        xsize = hd['NAXIS1']
        ysize = hd['NAXIS2']
        trimsection = [1,xsize,1,ysize] 
        xoff = 0 
        yoff = 0 
     
    # Delete the temporary FITS files
    if os.path.exists(tempfits): os.remove(tempfits)
     
     
     
    ############################################ 
    # STEP D: MAKE BAD PIXEL/WEIGHT MAP 
    # in how many images does the pixel need to be bad?? 
    # 1. shift the masks (created by FIXIMAGES.PRO) 
    #     using the shifts from IRAF_IMALIGN 
    # 2. trim the masks 
    # 3. combine the masks 
    logger.info('')
    logger.info('Step D: Making Bad Pixel Mask')
    logger.info('-----------------------------')
     
    # Make lists 
    maskfiles = os.path.dirname(tempfits)+'/'+os.path.basename(tempfits,'.fits')+'.mask.fits' 
    outmaskfiles = os.path.dirname(tempfits)+'/'+os.path.basename(tempfits,'.fits')+'.mask.shft.fits' 
    maskinfile = mchbase+'.maskinlist' 
    maskoutfile = mchbase+'.maskoutlist' 
    maskshiftsfile = mchbase+'.maskshifts'
    dln.writelines(maskinfile,maskfiles)
    dln.writelines(maskoutfile,outmaskfiles)
    for f in outmaskfiles:
        if os.path.exists(f): os.remove(f)
    strxyshifts = np.zeros(nfiles,(np.str,500))
    for i in range(nfiles): 
        strxyshifts[i] = '  '.join(xyshifts[i,:])
    dln.writelines(maskshiftsfile,strxyshifts)
              
    # Run IMSHIFT 
    #  set boundary to 0=BAD 
    logger.info('Shifting masks')
    iraflines = []
    iraflines += ['print("")']  # first line will be ignored 
    iraflines += ['cd '+curdir]
    iraflines += ['images']
    iraflines += ['imgeom']
    iraflines += ['imshift("@'+maskinfile+'","@'+maskoutfile+'",shifts_file="'+maskshiftsfile+'",'+               'interp_type="linear",boundary_typ="constant",constant=0)']
    iraflines += ['logout']
    imshiftscript = curdir+'/'+mchbase+'.imshift' 
    dln.writelines(imshiftscript,iraflines)
    out = run(imshiftscript,irafdir,verbose=False)
     
    # Trim 
    if trimcomb: 
        logger.info('Trimming masks') 
        xstart = trimsection[0]-1  # should be same as xoff 
        xstop = trimsection[1]-1 
        ystart = trimsection[2]-1  # should be same as yoff 
        ystop = trimsection[3]-1 
         
        for i in range(nfiles):
            im,head = io.readfile(outmaskfiles[i])
            shape = im.shape
            newim = im[xstart:xstop,ystart:ystop]
            # Add LTV1/LTV2 to the header 
            #  these are IRAF keywords to convert from logical to physical coords 
            ltv1 = head.get('LTV1')
            if ltv1 is None: ltv1=0
            ltv2 = head.get('LTV2')
            if ltv2 is None: ltv2=0
            head['LTV1'] = ltv1-xstart
            head['LTV2'] = ltv2-ystart
            fits.PrimaryHDU(newim,head).writeto(outmaskfiles[i],overwrite=True)
     
    # Combining masks 
    logger.info('Combining masks')
    bpm = None
    for i in range(nfiles):
        im,head = io.readfile(outmaskfiles[i])
        if i==0: 
            bpm = np.copy(im) 
            whead = np.copy(head)
        else: 
            bpm += im 
    #bpm = bpm/float(nfiles) 
     
    # masks have 0-bad, 1-good. 
    # anything with less than 1.0 is considered bad 
    # weight map, -1 is bad, +1 is good 
    #weightmap = -2.0*float(bpm lt nfiles) + 1. 
    #combweightfile = mchbase+'_comb.mask.fits' 
    #FITS_WRITE,combweightfile,weightmap,whead 
    # 
    # THIS IS NOW DONE BELOW AFTER THE IMAGE IS COMBINED 
    # DEPENDING ON IF THE IMAGES ARE SCALED OR NOT!!! 
     
     
    ############################################ 
    # STEP E: COMBINE IMAGES 
    logger.info('')
    logger.info('Step E: IMCOMBINE')
    logger.info('-----------------')
     
    # SCALE the images for combining 
    #-------------------------------
    if nocmbimscale==False:         
        # Put BPM mask names in file headers 
        #  these will be used by IMCOMBINE 
        for i in range(nfiles):
            head = os.readfile(outfiles[i],header=True)
            head['BPM'] = outmaskfiles[i]
            MODFITS,outfiles[i],0,head 
         
        # Combine the frames WITH scaling/offset/masking, for the bright stars 
        #logger.info('Creating SCALED image' 
        combfile = mchbase+'_comb.fits'
        for f in [combfile,mchbase+'_comb.bpm.pl']:
            if os.path.exists(f): os.remove(f)
        iraf.imcombine('@'+outfile,combfile,combine='average',reject='avsigclip',
                       weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',
                       irafdir=irafdir,scale='@'+scalefile,zero='@'+zerofile,
                       masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm')

         
        # Convert BPM mask from PL to FITS
        if os.path.exists(mchbase+'_comb.bpm.fits'): os.remove(mchbase+'_comb.bpm.fits')
        lines = []
        curdir = os.getcwd()
        lines += ['print("")']  # first line will be ignored 
        lines += ['cd '+curdir]
        lines += ['imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits']
        lines += ['logout']
        tmpfile = mktemp('tiraf') 
        dln.writelines(tmpfile,lines)
        out = iraf.run(tmpfile,irafdir,verbose=verbose)
         
        # Delete temporary scripts and PL file
        for f in [tmpfile,mchbase+'_comb.bpm.pl']:
            if os.path.exists(f): os.remove(f)
         
         
        # Fix the rdnoise and background/sky level and saturate 
        #  the bad pixels for DAOPHOT 
        #------------------------------------------------------ 
         
        # 10/02/12 
        # THE IMAGES ARE (1) ZERO-SUBTRACTED, (2) SCALED, AND (3) WEIGHT AVERAGED 
        # The algorithm is: 
        # 1.) add zero-level correction.  im = im+zero 
        # 2.) scale the images.  im = im*scale 
        # 3.) take weighted average.  combim=total(weight*im) 
        #      there is also clipping that takes place during the averaging 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2)) 
        # A gain that changes from frame to frame could be problematic, 
        # but this shouldn't happen since it's the same chip from the same night. 
         
        # IMCOMBINE wants rdnoise in electrons and gain in electrons/DN. 
        # DAOPHOT expects rdnoise in DN.  That's why mkopt converts 
        #  it with the gain.  So we are fine.  The header should have 
        #  rdnoise in ELECTRONS. 
         
        # page 63-65 of daophot2.pdf shows how you need to modify rdnoise/gain 
        # when averaging/summing frames. in observing/mosaic/. 
         
        # Load the IMCOMBINE output combined file and BPM
        combim,combhead = io.readfile(combfile)
        badmask,maskhead = io.readfile(mchbase+'_comb.bpm.fits')   # 0-good, 1-bad 
         
        # Fix the gain 
        # For N averaged frames gain(N)=N*gain(1) 
        # Leave the gain as is!  We are scaling everything to the reference 
        # and using its gain.  It's nearly impossible to figure out the real 
        # gain since we are scaling the images and then taking a weighted 
        # average with outlier rejection.  Find a gain that properly 
        # describes/follows Poisson errors for the final combined image is 
        # difficult/impossible.  But that's okay.  This is just for source 
        # detection and DAOPHOT FIND just cares about the noise in the 
        # background.  We just need to ensure that the sky and rdnoise 
        # are correct. 
         
        # Fix the rdnoise 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2)) 
        rdnoisearr = fltarr(nfiles) 
        for i in range(nfiles): 
            rdnoisearr[i] = utils.getrdnoise(base[i]+'.fits') 
        #  the "scales" array here is actually 1/scales used by IMCOMBINE. 
        rdnoise = np.sqrt(np.sum((weights*rdnoisearr/scales)**2)) 
        rdnoise = np.maximum(rdnoise, 0.01)  # must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
        dummy,rdnoisekey = utils.getrdnoise(combfile)  # get keyword 
        combhead[rdnoisekey] = rdnoise 
         
        # Fix the sky 
        # DAOPHOT FIND computes the random error per pixel in ADU as 
        # noise = sqrt( sky level/gain + rdnoise^2) 
        # So it assumes that the noise in the background is sqrt(sky/gain) 
        # in ADU.  We need to set the sky level so this is correct. 
        # The final noise should be 
        # final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
        # So the final sky level should be 
        # final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
        gain,gainkey = utils.getgain(combfile) 
        comb_sky = np.sum((weights*np.sqrt(np.maximum(sky,0)/gain)/scales)**2)*gain 
        # the "scales" array here is actually 1/scale 
        combim += float(comb_sky)# keep it float 
         
         
        # set the maximum to a "reasonable" level 
        # Rescale the image and increase the gain 
        if max(combim) > 50000: 
            rescale = 50000./np.max(combim) 
            combim = combim*rescale 
            combhead[gainkey] = gain/rescale 
            # rdnoise does NOT get modified since it's in electrons 
            # we just need to modify the gain which takes you from ADUs to electrons 
         
        maskdatalevel = np.max(combim) + 10000# set "bad" data level above the highest "good" value 
        combim2 = combim*(1-badmask) + maskdatalevel*badmask# set bad pixels to maskdatalevel
        fits.PrimaryHDU(combim2,combhead).writeto(combfile,overwrite=True) # fits_write can create an empty PDU 
         
        # Create the weight map for Sextractor using the BPM output by IMCOMBINE 
        #  bad only if bad in ALL images 
        weightmap = -2.0*float(badmask == 1) + 1.0 
        combweightfile = mchbase+'_comb.mask.fits'
        fits.PrimaryHDU(weightmap,whead).writeto(combweightfile,overwrite=True)
         
    # NO SCALING of the images for combining 
    #--------------------------------------- 
    else: 
        combfile = mchbase+'_comb.fits'
        if os.path.exists(combfile): os.remove(combfile)
        iraf.imcombine('@'+outfile,combfile,combine='average',reject='avsigclip',
                  weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',irafdir=irafdir)
         
        # Fix the rdnoise and background/sky level and saturate 
        #  the bad pixels for DAOPHOT 
        #------------------------------------------------------ 
         
        # See the explanations for all these steps above!! 
         
        # Load the IMCOMBINE output combined file and BPM 
        combim,combhead = io.readfile(combfile)
        badmask,maskhead = io.readfile(mchbase+'_comb.bpm.fits')  # 0-good, 1-bad 
         
        # Fix the rdnoise 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2)) 
        rdnoisearr = np.zeros(nfiles,float) 
        for i in range(nfiles): 
            rdnoisearr[i] = utils.getrdnoise(base[i]+'.fits') 
        rdnoise = np.sqrt(np.sum((weights*rdnoisearr)**2)) 
        rdnoise = np.maximum(rdnoise > 0.01)  # must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
        dummy,rdnoisekey = utils.getrdnoise(combfile)  # get keyword 
        combhead[rdnoisekey] = rdnoise 
         
        # Fix the sky 
        # So the final sky level should be 
        # final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
        gain,gainkey = utils.getgain(combfile) 
        comb_sky = np.sum((weights*np.sqrt(sky/gain))**2)*gain 
        # the "scales" array here is actually 1/scale 
        combim += comb_sky 
         
         
        # set the maximum to a "reasonable" level 
        # Rescale the image and increase the gain 
        if np.max(combim) > 50000: 
            rescale = 50000./np.max(combim) 
            combim = combim*rescale 
            combhead[gainkey] = gain/rescale 
            # rdnoise does NOT get modified since it's in electrons 
            # we just need to modify the gain which takes you from ADUs to electrons 
         
         
        # Making Sextractor "weight" map file 
        #------------------------------------ 
        # masks have 0-bad, 1-good. 
        # anything with less than 1.0 is considered bad 
        # weight map, -1 is bad, +1 is good 
        # "bpm" is the SUM of the bad pixel masks 
        # consider a pixel bad that is bad in ANY image 
        weightmap = -2.0*float(bpm < nfiles) + 1. 
        combweightfile = mchbase+'_comb.mask.fits' 
        fits.PrimaryHDU(weightmap,whead).writeto(combweightfile,overwrite=True)
         
        #--------------------------------------------- 
        # SATURATE BAD pixels in the COMBINED IMAGE 
        # DAOPHOT needs to have the bad pixels "saturated", 
        # SExtractor will know which pixels are bad from the "weight" map. 
        # 
        # We could skip the fiximage.pro step but we still need the 
        # individual bpm masks and setting the bad pixels to the background 
        # probably helps in the IMALIGN/IMCOMBINE steps. 
        logger.info('')
        logger.info('"Saturating" bad pixels in the COMBINED image')
        logger.info('') 
         
        badmask = float(weightmap < 0.5) 
        maskdatalevel = max(combim) + 10000# set "bad" data level above the highest "good" value 
        combim2 = combim*(1.0-badmask) + maskdatalevel*badmask# set bad pixels to 100,000
        fits.PrimaryHDU(combim2,combhead).writeto(combfile,overwrite=True)
     
    # Add TILETYPE to the combined image
    combhead = io.readfile(combfile,header=True)
    combhead['AFTILTYP'] = 'ORIG' 
    MODFITS,combfile,0,combhead 
     
    # Delete the shifted images
    shiftedfiles = dln.readlines(outfile)
    for f in shiftedfiles:
        if os.path.exists(f): os.remove(f)
     
    # Delete mask files
    for f in [maskfiles,outmaskfiles,maskinfile,maskoutfile,maskshiftsfile,imshiftscript]:
        if os.path.exists(f): os.remove(f)
     
    # Copy the original MCH file to COMB.MCH
    if os.path.exists(mchdir+'/'+mchbase+'.comb.mch'): os.remove(mchdir+'/'+mchbase+'.comb.mch')
    shutil.copyfile(filename,mchdir+'/'+mchbase+'.comb.mch')
     
    # Using CMN.LST of reference frame if it exists 
    if os.path.exists(mchbase+'.cmn.lst') and usecmn: 
        logger.info('Using reference image COMMON SOURCE file')
        if os.path.exists(mchbase+'_comb.cmn.lst'): os.remove(mchbase+'_comb.cmn.lst')
        shutil.copyfile(mchbase+'.cmn.lst',mchbase+'_comb.cmn.lst')
     
    # CD back to the original directory
    os.chdir(curdir)
 

    return maskdatalevel,fileinfo,xoff,yoff
