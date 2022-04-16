#!/usr/bin/env python

import os
import time
import numpy as np
import glob as glob
from astropy.io import fits
from astropy.wcs import WCS
from dlnpyutils import utils as dln
from . import io

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

 
def allframe_calcweights(mag,err,fwhm,rdnoise,medsky,actweight,scales):
    """
    Helper function to calculate weights
    """
     
    sz = size(mag) 
    ngd = sz[2] 
    nfiles = sz[1] 
     
    # Make 2D arrays for fwhm, rdnoise and medsky 
    fwhm2 = fwhm#(fltarr(ngd)+1.0) 
    rdnoise2 = rdnoise#(fltarr(ngd)+1.0) 
    medsky2 = medsky#(fltarr(ngd)+1.0) 
     
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
    weight = (intensity/(!dpi*fwhm2**2.0))/ (sqrt(rdnoise2**2.0 + medsky2 + intensity/(!dpi*fwhm2**2.0) ) ) 
    #weight(i,n)=(intensity/(Pi*FWHM(i)**2)) / (sqrt(RDNOISE(i)**2+ MEDSKY(i) + intensity/(Pi*FWHM(i)**2) ) ) 
     
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
     
    maxw = np.max(weight,dim=1) 
    maxw2 = (np.ones(nfiles))#maxw 
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
     
    #avgweight = total(weight2,2)/ngd 
    avgweight = np.sum(nweight,2) 
    actweight = avgweight/np.sum(avgweight) 
    # They sum to 1.0 
     
     
    # Compute the scaling for each frame.  scale=1 for the first frame 
    scales = fltarr(nfiles) 
    for i in range(nfiles): 
        ratio = reform(intensity[i,:]/intensity[0,:]) 
        ratio_error = reform(sqrt(err[i,:]**2 + err[0,:]**2))# mag errors are ~fractional 
        med_ratio = np.median(ratio) 
        #wmeanerr,ratio,ratio_error,xmean,xsigma 
        scales[i] = med_ratio 
     
    return scales


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

     
    mag = transpose(str.mag) 
    err = transpose(str.err) 
    fwhm = str.fwhm 
    rdnoise = str.rdnoise 
    medsky = str.medsky 
     
    if size(mag,/n_dim) == 2: 
        nstars = len(mag[0,:]) 
    else: 
        nstars=len(mag) 
    if size(mag,/n_dim) == 2: 
        nfiles = len(mag[:,0]) 
    else: 
        nfiles=1 
     
    # Getting the reference sources 
    totstars = np.sum(mag < 50,1) 
    si = np.flip(np.argsort(totstars))# get the stars with the most detections 
    #gdrefstars = si[0:(99<(nstars-1))] 
    gdrefstars = si[0:(49<(nstars-1))] 
    nrefstars = len(gdrefstars) 
    # Getting the "good" frames 
    totframe = np.sum(mag[:,gdrefstars] < 50,2)# # of ref stars detected per frame 
    gdframe , = np.where(totframe == nrefstars,ngdframe,comp=bdframe,ncomp=nbdframe) 
    # No good frames, lower number of reference stars 
    if ngdframe == 0: 
        gdrefstars = si[0:(29<(nstars-1))] 
        nrefstars = len(gdrefstars) 
        totframe = np.sum(mag[:,gdrefstars] < 50,2) 
        gdframe , = np.where(totframe == nrefstars,ngdframe,comp=bdframe,ncomp=nbdframe) 
    # No good frames again, say the frame with the most detections is "good" 
    #   get weights relative to that one for the others 
    if ngdframe == 0: 
        # say the frame with the most detections is "good" and 
        #  the rest are bad, get weights relative to this one frame 
        totstarsframe = np.sum(mag < 50,2) 
        gdframe = first_el(maxloc(totstarsframe)) 
        bdframe = lindgen(nfiles,1) 
        remove,gdframe,bdframe 
        nbdframe = len(bdframe) 
        # Get stars that are good in this frames and in ALL of the others 
        gdrefstars , = np.where(reform(mag[gdframe,:]) < 50 and totstars == nfiles,nrefstars) 
        #  lower threshold, at least half 
        if nrefstars == 0 : 
            gdrefstars , = np.where(reform(mag[gdframe,:]) < 50 and totstars > 0.5*nfiles,nrefstars) 
        #  just the good ones 
        if nrefstars == 0: 
            gdrefstars , = np.where(reform(mag[gdframe,:]) < 50,nrefstars) 
            si = np.flip(np.argsort(totstars[gdrefstars]))# order by how many other frames they are detected in 
            gdrefstars = gdrefstars[si[0:(49>(nrefstars-1))]]# only want 50 
            nrefstars = len(gdrefstars) 
     
    # Calculate the weights 
    weights = fltarr(nfiles) 
    scales = fltarr(nfiles) 
    mag2 = mag[gdframe,:] & mag2 = mag2[:,gdrefstars] 
    err2 = err[gdframe,:] & err2 = err2[:,gdrefstars] 
    ALLFRAME_CALCWEIGHTS,mag2,err2,fwhm[gdframe],rdnoise[gdframe],medsky[gdframe],         weights1,scales1 
    weights[gdframe] = weights1 
    scales[gdframe] = scales1 
     
    # If there are some "bad" frames calculate their weights 
    #  relative to some of the "good" ones 
    for i in range(nbdframe): 
         
        iframe = bdframe[i] 
         
        # First round of "good" stars 
        totgdstars = np.sum(mag[gdframe,:] < 50,1)# stars in good frames 
        igdrefstars1 , = np.where(mag[iframe,:] < 50 and totgdstars >= 5,nigdrefstars1) 
        if nigdrefstars1 < 3 : 
            igdrefstars1 , = np.where(mag[iframe,:] < 50 and totgdstars >= 1,nigdrefstars1) 
        if nigdrefstars1 < 2 : 
            goto,BOMB 
         
        totgdstars1 = np.sum( (mag[gdframe,:])[:,igdrefstars1] < 50,1) 
        si1 = np.flip(np.argsort(totgdstars1))# get the best ones 
        igdrefstars = igdrefstars1[si1[0:(49<(nigdrefstars1-1))]] 
        nirefstars = len(igdrefstars) 
         
        itotframe = np.sum( (mag[gdframe,:])[:,igdrefstars] < 50,2) 
        igdframe1 , = np.where( itotframe == nirefstars,nigdframe1) 
         
        # Whittle down to the best stars/frames 
        if nigdframe1 == 0: 
             
            igdframe = gdframe 
            igdrefstars = igdrefstars 
             
            # whittle down to best stars and frames 
            count = 0 
= 0 
        while endflag == 0: 
             
            # sum over frames 
            itot1 = np.sum( (mag[igdframe,:])[:,igdrefstars] < 50,1) 
            si1 = np.argsort(itot1) 
            # remove worst 5% stars 
            if len(igdrefstars) > 3: 
                bd1 = si1[0:int(np.round(len(igdrefstars)*0.05)] 
                remove,bd1,igdrefstars 
             
            # sum over stars 
            itot2 = np.sum( (mag[igdframe,:])[:,igdrefstars] < 50,2) 
            si2 = np.argsort(itot2) 
            # remove worst 5% frames 
            if len(igdframe) > 1: 
                bd2 = si2[0:int(np.round(len(igdframe)*0.05)] 
                remove,bd2,igdframe 
             
            # Testing if we should end 
            itotframe = np.sum( (mag[igdframe,:])[:,igdrefstars] < 50,2) 
            igdframe1 , = np.where( itotframe == len(igdrefstars),nigdframe1) 
            if nigdframe1 > 0 : 
 
         
        if endflag == 1 : 
            igdframe=igdframe[igdframe1] 
        if count > 50 :# not converging, end it 
            goto,BOMB 
         
        count++ 
     
 else igdframe=gdframe[igdframe1] 
     
    # Get weights relative to some "good" frames 
    indframes = [igdframe,iframe] 
    mag3 = (mag[indframes,:])[:,igdrefstars] 
    err3 = (err[indframes,:])[:,igdrefstars] 
    ALLFRAME_CALCWEIGHTS,mag3,err3,fwhm[indframes],rdnoise[indframes],medsky[indframes],                       weights3,scales3 
     
    # Scale these relative to the original ones 
    weights3a = weights[igdframe]# original 
    weights3b = weights3[0:nigdframe1-1]# new 
    wtfrac = np.median(weights3a/weights3b) 
    scales3a = scales[igdframe]# original 
    scales3b = scales3[0:nigdframe1-1]# new 
    sclfrac = np.median(scales3a/scales3b) 
    new_weights = weights3[nigdframe1] * wtfrac 
    new_scale = scales3[nigdframe1] * sclfrac 
    weights[iframe] = new_weights 
    scales[iframe] = new_scale 
     
    #print,iframe,new_weights,new_scale 
     
    BOMB: 
     
 
# Fix images with bad weights likely due to no overlap 
#  use the FLUX values to get a weight 
bdweights , = np.where(weights <= 0.0,nbdweights,comp=gdweights,ncomp=ngdweights) 
if nbdweights > 0: 
    # Use fluxrate10 to get weights and flux10 to get scales 
    # relative to "good" frame 
    if ngdweights > 0: 
        weights[bdweights] = str[bdweights].fluxrate10 * np.median([weights[gdweights]/str[gdweights].fluxrate10]) 
        scales[bdweights] = str[bdweights].flux10 * np.median([scales[gdweights]/str[gdweights].flux10]) 
        # all bad 
    else: 
        weights = str.fluxrate10 
        scales = str.flux10 
 
# Normalize the weights 
weights /= np.sum(weights) 
 
# Rescale SCALES so the reference frames has SCALE=1.0 
if scales[0] > 0.0: 
    scales /= scales[0] 
else: 
    scales/=max(scales) 
 
# Create the output structure 
schema = str[0] 
struct_assign,{dum:''},schema 
if tag_exist(schema,'weight') == 0 : 
    schema = create_struct(schema,'weight',0.0) 
if tag_exist(schema,'scale') == 0 : 
    schema = create_struct(schema,'scale',0.0) 
outstr = replicate(schema,nfiles) 
struct_assign,str,outstr 
    outstr['weight'] = weights 
    outstr['scale'] = scales 


    return outstr
 
def getweights(mchfile,imager=None,logfile=None,silent=False):
    """
    This calculates weights, scales and sky values from DAOPHOT 
    photometry files for images combination.  Used by ALLFRAME.PRO and 
    other programs. 
 
    You need to have run daophot, allstar, daomatch and daomaster 
    already.  There need als, mch and raw files. 
 
    Parameters
    ----------
    mchfile     The MCH filename 
    imager     Imager structure with basic information 
    logfile    A logfile to print to output to. 
    silent     Don't print anything to the screen. 
 
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
    
    global photred,setup
     
    # OUTPUTS: 
    #  actweight  The weight for each frame 
    #  scales     The scale for each frame 
    #  medsky     The sky value for each frame 
     
    tilesep = '+' 
    #tilesep = '.' 
    btilesep = int(byte(tilesep)) 
     
    # MCH file not found 
    if os.path.exists(mchfile) == 0: 
        print(mchfile,' NOT FOUND' 
        return 
     
    if len(logfile) == 0 : 
        logfile=-1 
     
    # Load the MCH file 
    LOADMCH,mchfile,files,trans 
    nfiles = len(files) 
     
    #----------------------------------- 
    # Get the information that we need 
     
    # Load the opt files 
    info = replicate({name:'',filter:'',exptime:0.0,fwhm:0.0,rdnoise:0.0,mnsky:0.0,medsky:0.0,
                      mag10:99.99,flux10:0.0,fluxrate10:0.0,weight:0.0,scale:0.0},nfiles) 
    info.name = files 
    for i in range(nfiles): 
        dir = os.path.dirname(mchfile) 
        base = os.path.basename(files[i],'.als') 
        optfile = dir+'/'+base+'.opt' 
        logfile1 = dir+'/'+base+'.log' 
         
        READCOL,optfile,name,dum,value,format='A,A,F',/silent 
        name = str(strupcase(name),2) 
         
        ind_re , = np.where(name == 'RE',nind_re) 
        info[i].rdnoise = value[ind_re[0]] 
        ind_fw , = np.where(name == 'FW',nind_fw) 
        info[i].fwhm = value[ind_fw[0]] 
         
        SPAWN,['grep','Clipped',logfile1],out,errout,/noshell 
        #              Clipped mean and median =  187.442  187.215 
         
        # daophot.sh log files are clipped on Tortoise for some reason 
        #  Get mean/median sky level 
        if len(out) == 1 and out[0] == '': 
            print('Getting mean/median sky levels for ',base 
            if os.path.exists('daophot.opt') == 0 : 
                file_copy,base+'.opt','daophot.opt'
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
            tempscript = MKTEMP('dfind')# absolute filename 
            dln.writelines(tempscript,cmdlines)
            os.chmod(tempscripts,0o755)
            os.remove(base+'.find.temp',/allow)
            os.remove(base+'.find.log',/allow)
            # Run DAOPHOT/FIND 
            out1 = subprocess.check_output(tempscript+' '+base)
             
            logfile2 = base+'.find.log'
            out = subprocess.check_output(['grep','Clipped',logfile2],shell=False)
            #              Clipped mean and median =  187.442  187.215 
             
            # Delete temporary files 
            os.remove(base+'.find.temp',/allow)
            os.remove(tempscript,/allow_
         
        arr = strsplit(out[0],' ',/extract) 
        info[i].mnsky = float(arr[5]) 
        info[i].medsky = float(arr[6]) 
         
        # Get exptime and filter 
        fitsfile = base+'.fits' 
        if os.path.exists(fitsfile) == 0 : 
            fitsfile+='.fz' 
        info[i].exptime = PHOTRED_GETEXPTIME(fitsfile) 
        info[i].filter = PHOTRED_GETFILTER(fitsfile) 

     
    # Only ONE file, return 1s 
    if nfiles == 1: 
        actweight = 1.0 
        scales = 1.0 
        medsky = info[0].medsky 
        return 
     
     
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
    dir = os.path.dirname(mchfile) 
    base = os.path.basename(mchfile,'.mch') 
    rawfile = dir+'/'+base+'.raw' 
    LOADRAW,rawfile,raw,rawhead 
     
    # Making mag and err arrays 
    nstars = len(raw) 
    tags = tag_names(raw) 
    magind , = np.where(stregex(tags,'MAG',/boolean) == 1,nmagind) 
    errind , = np.where(stregex(tags,'ERR',/boolean) == 1,nerrind) 
    mag = fltarr(nfiles,nstars) 
    err = fltarr(nfiles,nstars) 
    for i in range(nfiles): 
        mag[i,:] = raw.(magind[i]) 
    for i in range(nfiles): 
        err[i,:] = raw.(errind[i]) 
     
    # Calculate the magnitude and flux at the 10sigma magnitude 
    for i in range(nfiles): 
        gd , = np.where(mag[i,:] < 50,ngd) 
        if ngd > 0: 
            mag1 = mag[i,gd] 
            snr1 = 1.087/err[i,gd] 
            gdsnr , = np.where(abs(snr1-10.0) < 1,ngdsnr) 
            if ngdsnr < 5 : 
                gdsnr , = np.where(abs(snr1-10.0) < 2,ngdsnr) 
            if ngdsnr < 5: 
                si = np.argsort(abs(snr1-10.0)) 
                gdsnr = si[0:99<(ngd-1)] 
                ngdsnr = len(gdsnr) 
            mag10 = np.median([mag1[gdsnr]]) 
            info[i].mag10 = mag10 
            info[i].flux10 = 10.0**( (mag10-25.0)/(-2.5) )# total counts 
            info[i].fluxrate10 = info[i].flux10 / info[i].exptime# flux rate = counts / sec 
     
     
    # Using TILES 
    #-------------- 
    # We are using TILES and have multiple chips/amps 
    #   F1-00507800_39+T2.als, '+T' and two dots 
    if len(imager) > 0: 
        namps = imager['namps']
    else: 
        namps = 1 
    if np.sum(stregex(files,'\'+tilesep+'T',/boolean)) == nfiles and    np.sum(int(byte(files[0])) == btilesep) >= 2 and namps > 1: 
        usetiles = 1 
         
        # Number of unique exposures 
        expname = strarr(nfiles) 
        chip = strarr(nfiles) 
        for i in range(nfiles): 
            base1 = os.path.basename(files[i],'.als')# F1-00507800_39+T2 
            field1 = (strsplit(base1,'-',/extract))[0]# F1 
            expchptile = (strsplit(base1,'-',/extract))[1]# 00507800_39+T2 
            expchp = (strsplit(expchptile,tilesep,/extract))[0]# 00507800_39 
            expname[i] = (strsplit(expchp,imager.separator,/extract))[0] 
            chip[i] = (strsplit(expchp,imager.separator,/extract))[1] 
        # Unique exposures 
        uiexp = np.uniq(expname,np.argsort(expname)) 
        uexpname = expname[uiexp] 
        nexp = len(uexpname) 
         
        # Combine all the catalogs for a given exposure 
        # Calculate the weights 
        expstr = replicate({mag:fltarr(nstars),err:fltarr(nstars),nfiles:0L,index:lonarr(nfiles),
                      exptime:0.0,filter:'',fwhm:0.0,rdnoise:0.0,medsky:0.0,mag10:99.99,flux10:0.0,fluxrate10:0.0},nexp) 
        expmag = fltarr(nexp,nstars) 
        experr = fltarr(nexp,nstars) 
        for i in range(nexp): 
            ind , = np.where(expname == uexpname[i],nind) 
            expstr[i].nfiles = nind 
            expstr[i].index[0:nind-1] = ind 
            expstr[i].filter = info[ind[0]].filter 
            expstr[i].exptime = info[ind[0]].exptime 
            expstr[i].fwhm = np.median([info[ind].fwhm]) 
            expstr[i].rdnoise = np.median([info[ind].rdnoise]) 
            expstr[i].medsky = np.median([info[ind].medsky]) 
            expstr[i].mag10 = np.median([info[ind].mag10]) 
            expstr[i].flux10 = np.median([info[ind].flux10]) 
            expstr[i].fluxrate10 = np.median([info[ind].fluxrate10]) 
            # Combine the photometry 
            mag1 = mag[ind,:] 
            err1 = err[ind,:] 
            # Multiple chips 
            #   they shouldn't overlap, so just use the mean/median 
            #   and that should pick up the detections 
            if nind > 1: 
                bd , = np.where(mag1 > 50,nbd) 
                if nbd > 0 : 
                    mag1[bd]=!values.f_nan 
                if nbd > 0 : 
                    err1[bd]=!values.f_nan 
                expmag[i,:] = np.median(mag1,dim=1) 
                experr[i,:] = np.median(err1,dim=1) 
            else: 
                expmag[i,:] = mag1 
                experr[i,:] = err1 
        # Replace NANs with 99.9999 
        bdmag , = np.where(finite(expmag) == 0,nbdmag) 
        expmag[bdmag] = 99.99 
        experr[bdmag] = 9.99 
        expstr.mag = transpose(expmag) 
        expstr.err = transpose(experr) 
        # Perform the weights and scales calculations 
        outexpstr = getweights_raw(expstr)
        # Give each chip the weight of its respective exposure 
        for i in range(nexp): 
            info[expstr[i].index[0:expstr[i].nfiles-1]].weight = outexpstr[i].weight 
            info[expstr[i].index[0:expstr[i].nfiles-1]].scale = outexpstr[i].scale 
         
    # REGULAR Method 
    #--------------- 
    else: 
        tab = replicate({mag:fltarr(nstars),err:fltarr(nstars),
                         exptime:0.0,filter:'',fwhm:0.0,rdnoise:0.0,
                         medsky:0.0,flux10:0.0,fluxrate10:0.0},nfiles) 
        struct_assign,info,str 
        tab['mag'] = transpose(mag) 
        tab['err'] = transpose(err) 
        # Perform the weights and scales calculations 
        outstr = getweights_raw(tab)
        info['weight'] = outstr['weight'] 
        info['scale'] = outstr['scale']
     
     
    # Print out the information 
    if silent==False:
        logger.info('        FILE         FILTER EXPTIME FWHM RDNOISE MEDSKY WEIGHT SCALE')
        for i in range(nfiles):
            fmt = '(A-23,A4,F6.1,F6.2,F6.2,F8.1,F7.3,F7.3)'
            data = (info[i].name,info[i].filter,info[i].exptime,info[i].fwhm,info[i].rdnoise,
                    info[i].medsky,info[i].weight,info[i].scale)
            logger.info(fmt % data)
     
    # Final output 
    actweight = info['weight']
    scales = info['scale']
    medsky = info['medsky']

    return actweight,scales,medsky

                      
def combine(filename,tile=None,scriptsdir=None,logfile=None,irafdir=None,
            satlevel=6e4,nocmbimscale=False,fake=False,usecmn=False,imager=None):
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
    logfile : str
       A logfile to print to output to. 
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
    filestr : table
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

    global photred,setup 

    # Logfile 
    if keyword_set(logfile): 
        logf=logfile 
    else: 
        logf=-1 
     
    # Getting scripts directory and iraf directory 
    nsetup = len(setup) 
    if nsetup > 0: 
        scriptsdir = setup['SCRIPTSDIR']
        irafdir = setup['IRAFDIR']
     
    # No irafdir 
    if irafdir is None:
        raise ValueError('IRAFDIR NOT INPUT')
    # No irafdir 
    if scriptsdir is None:
        raise ValueError('SCRIPTSDIR NOT INPUT')    
     
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if check_iraf(iraftest,irafdir=irafdir)==False:
        raise ValueError('IRAF TEST FAILED.  EXITING')

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
     
    # FILENAME 
    mchfile = os.path.basename(filename) 
    mchdir = os.path.dirname(filename) 
    mchbase = os.path.splitext(filename)[0]
     
    # CD to the directory
    curdir = os.path.abspath(os.getcwd())
    os.chdir(mchdir)



    
    # Check that the mch, als, and opt files exist 
    mchtest = os.path.exists(mchfile) 
    if mchtest == 0: 
        logger.info(mchfile,' NOT FOUND' 
        return 
     
    # Checking RAW file 
    rawtest = os.path.exists(mchbase+'.raw') 
    if rawtest == 0: 
        logger.info(mchbase+'.raw NOT FOUND' 
        return 
     
     
    ############################################ 
    # CHECK NECESSARY FILES 
     
    # Load the MCH file
    files,trans,magoff = io.loadmch(mchfile)
    nfiles = len(files) 
     
    # FAKE, check that we have all the files that we need 
    if fake: 
        # weights, scale, zero, comb_psf, _shift.mch 
        chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_comb.mch'] 
        bdfiles , = np.where(os.path.exists(chkfiles) == 0,nbdfiles) 
        if nbdfiles > 0: 
            error = 'FAKE.  Some necessary files not found. '+strjoin(chkfiles[bdfiles],' ') 
            logger.info(error 
            return 
     
     
    # Final files 
    base = os.path.basename(files,'.als') 
    fitsfiles = base+'.fits' 
    outfiles = base+'.shft.fits' 
    outmaskfiles = base+'.mask.shft.fits' 
     
     
    # Gather information on all of the files 
    # photred_gatherfileinfo.pro can do most of this 
    logger.info('Gathering file information' 
    ntrans = len(trans[0,:]) 
    filestr = replicate({fitsfile:'',catfile:'',nx:0L,ny:0L,trans:dblarr(ntrans),magoff:fltarr(2),head:ptr_new(),
                         vertices_ra:dblarr(4),vertices_dec:dblarr(4),pixscale:0.0,saturate:0.0,
                         background:0.0,comb_zero:0.0,comb_scale:0.0,comb_weights:0.0,
                         resampfile:'',resampmask:'',resamptrans:dblarr(ntrans),resamptransrms:0.0},nfiles) 
    filestr['fitsfile'] = fitsfiles 
    filestr['catfile'] = files 
    filestr['trans'] = transpose(trans) 
    filestr['magoff'] = magoff 
    filestr['resampfile'] = outfiles 
    filestr['resampmask'] = outmaskfiles 
    for i in range(nfiles): 
        im1 = PHOTRED_READFILE(filestr[i].fitsfile,head1) 
        filestr[i].head = ptr_new(head1) 
        filestr[i].nx = sxpar(head1,'NAXIS1') 
        filestr[i].ny = sxpar(head1,'NAXIS2') 
        HEAD_XYAD,head1,[0,filestr[i].nx-1,filestr[i].nx-1,0],[0,0,filestr[i].ny-1,filestr[i].ny-1],vra,vdec,/degree 
        filestr[i].vertices_ra = vra 
        filestr[i].vertices_dec = vdec 
        GETPIXSCALE,'',pixscale,head=head1 
        filestr[i].pixscale = pixscale 
        saturate = sxpar(head1,'SATURATE',count=nsaturate,/silent) 
        if nsaturate == 0 : 
            saturate=50000L 
        filestr[i].saturate = saturate 
        gdpix , = np.where(im1 < saturate,ngdpix,ncomp=nbdpix) 
        background = np.median(im1[gdpix]) 
        filestr[i].background = background 
     
     
    ################################################## 
    # Create default reference frame if TILE not input 
    if len(tileinp) > 0: 
        tile = tileinp 
    else: 
        tile = {'type':'WCS'} 
    if tile.type == 'WCS' and n_tags(tile) == 1: 
        logger.info('Creating TILE projection' 
        # The default projection is a tangent plane centered 
        # at halfway between the ra/dec min/max of all of 
        # the images.  The mean pixel scale is used. 
        #  near RA=0 line 
        if range(filestr.vertices_ra) > 180: 
            vertices_ra = filestr.vertices_ra 
            over , = np.where(vertices_ra > 180,nover,comp=under,ncomp=nunder) 
            if nover > 0 : 
                vertices_ra[over]-=360 
            rar = minmax(vertices_ra) 
            cenra = np.mean(rar) 
        else: 
            rar = dln.minmax(filestr['vertices_ra']) 
            cenra = np.mean(rar) 
        decr = dln.minmax(filestr['vertices_dec']) 
        cendec = np.mean(decr) 
        pixscale = np.mean(filestr.pixscale) 
        # Set up the tangent plane projection 
        step = pixscale/3600.0d0 
        delta_dec = dln.valrange(decr) 
        delta_ra = dln.valrange(rar)*np.cos(np.deg2rad(cendec))
        nx = np.ceil(delta_ra*1.01/step) 
        ny = np.ceil(delta_dec*1.01/step) 
        xref = nx/2 
        yref = ny/2
        tilehead = fits.Header()
        tilehead['NAXIS1'] = nx 
        tilehead['CDELT1'] = step 
        tilehead['CRPIX1'] = xref+1L 
        tilehead['CRVAL1'] = cenra 
        tilehead['CTYPE1'] = 'RA---TAN' 
        tilehead['NAXIS2'] = ny 
        tilehead['CDELT2'] = step 
        tilehead['CRPIX2'] = yref+1L 
        tilehead['CRVAL2'] = cendec 
        tilehead['CTYPE2'] = 'DEC--TAN' 
        tilewcs = WCS(tilehead)
        tileast['equinox'] = 2000 
         
        logger.info('RA range = [',str(rar[0]),',',str(rar[1]),'] deg')
        logger.info('DEC range = [',str(decr[0]),',',str(decr[1]),'] deg')
        logger.info('Central RA = ',str(cenra))
        logger.info('Central DEC = ',str(cendec)) 
        logger.info('NX = ',str(nx))
        logger.info('NY = ',str(ny))
         
        # Create the TILE structure 
        tile = {'type':'WCS','naxis':int([nx,ny]),'cdelt':double([step,step]),'crpix':double([xref+1L,yref+1L]),
                'crval':double([cenra,cendec]),'ctype':['RA--TAN','DEC--TAN'],
                'head':tilehead,'wcs':tilewcs,'xrange':[0,nx-1],'yrange':[0,ny-1],'nx':nx,'ny':ny} 
     
    # Check that the TILE is valid 
    if allframe_validtile(tile,error=tilerror) == 0: 
        error = tilerror 
        if not keyword_set(silent) : 
            logger.info(error 
        return 
     
    # Add HEAD/AST to TILE if needed 
    if tag_exist(tile,'HEAD') == 0: 
        MKHDR,tilehead,fltarr(5,5) 
        tilehead['NAXIS1'] = tile.naxis[0] 
        tilehead['CDELT1'] = tile.cdelt[0] 
        tilehead['CRPIX1'] = tile.crpix[0] 
        tilehead['CRVAL1'] = tile.crval[0] 
        tilehead['CTYPE1'] = tile.ctype[0] 
        tilehead['NAXIS2'] = tile.naxis[1] 
        tilehead['CDELT2'] = tile.cdelt[1] 
        tilehead['CRPIX2'] = tile.crpix[1] 
        tilehead['CRVAL2'] = tile.crval[1] 
        tilehead['CTYPE2'] = tile.ctype[1] 
        if 'CD' in tile: 
            tilehead['CD1_1'] = tile.cd[0,0] 
            tilehead['CD1_2'] = tile.cd[0,1] 
            tilehead['CD2_1'] = tile.cd[1,0] 
            tilehead['CD2_2'] = tile.cd[1,1] 
        EXTAST,tilehead,tileast 
        tileast.equinox = 2000 
        tile = CREATE_STRUCT(tile,'HEAD',tilehead,'AST',tileast) 
    if 'AST' not in tile: 
        EXTAST,tile.head,tileast 
        tileast.equinox = 2000 
        tile = CREATE_STRUCT(tile,'AST',ast) 
    # Add XRANGE/YRANGE 
    if 'XRANGE' not in tile: 
        tile = CREATE_STRUCT(tile,'XRANGE',int([0,tile.naxis[0]-1]),'NX',tile.naxis[0]) 
    if 'YRANGE' not in tile: 
        tile = CREATE_STRUCT(tile,'YRANGE',int([0,tile.naxis[1]-1]),'NY',tile.naxis[1]) 
     
     
    ############################################ 
    # STEP 1: IMALIGN PREP 
     
    logger.info('------------------------')
    logger.info('STEP 1: GETTING WEIGHTS')
    logger.info('------------------------')
     
    #----------------------------------- 
    # Computs Weights 
    if fake==False:
        weights,scales,sky = getweights(mchfile,imager=imager,logfile=logf) #,raw2=raw2 
        invscales = 1.0/scales 
        bdscale, = np.where((scales < 1e-5) | (invscales > 900))
        if len(bdscale) > 0: 
            scales[bdscale] = 1.0 
            invscales[bdscale] = 1.0 
            weights[bdscale] = 0.0 
        weightfile = mchbase+'.weights' 
        WRITECOL,weightfile,weights,fmt='(F10.6)' 
        scalefile = mchbase+'.scale' 
        WRITECOL,scalefile,invscales,fmt='(F10.6)'# want to scale it UP 
        zerofile = mchbase+'.zero' 
        WRITECOL,zerofile,-sky,fmt='(F12.4)'# want to remove the background, set to 1st frame 
         
    # FAKE, use existing ones 
    else: 
        weightfile = mchbase+'.weights' 
        scalefile = mchbase+'.scale' 
        zerofile = mchbase+'.zero' 
        READCOL,weightfile,weights,format='F',/silent 
        READCOL,scalefile,invscales,format='F',/silent 
        scales = 1.0/invscales 
        READCOL,zerofile,sky,format='F',/silent 
        sky = -sky 
    # Stuff the information into the FILEINFO structure 
    filestr['comb_weights'] = weights 
    filestr['comb_scale'] = invscales 
    filestr['comb_zero'] = -sky 
     
     
    ############################################ 
    # STEP 2: Resample the images 
     
    logger.info('---------------------------' 
    logger.info('STEP 2: Resample the images' 
    logger.info('---------------------------' 
     
    CASE tile.type of 
         
        # --- WCS projection on the sky --- 
        'WCS': begin 
         
        # Convert x/y and ra/dec grid for entire tile image 
        #  that we'll be using, ~15sec 
        xb = (lindgen(tile.nx)+tile.xrange[0])#replicate(1,tile.ny) 
        yb = replicate(1,tile.nx)#(lindgen(tile.ny)+tile.yrange[0]) 
        HEAD_XYAD,tile.head,xb,yb,rab,decb,/deg 
         
        # Loop through the files 
        for i in range(nfiles): 
            im1 = PHOTRED_READFILE(filestr[i].fitsfile,head1) 
             
            # Make the mask 
            mask = bytarr(filestr[i].nx,filestr[i].ny)+1 
            gdpix , = np.where(im1 < filestr[i].saturate,ngdpix,comp=bdpix,ncomp=nbdpix) 
            # 0-bad, 1-good 
            if nbdpix > 0: 
                mask[bdpix] = 0 
                im1[bdpix] = filestr[i].background 
             
            # Get X/Y range for this image in the final coordinate system 
            HEAD_ADXY,tile.head,filestr[i].vertices_ra,filestr[i].vertices_dec,vx,vy,/deg 
            xout = [floor(min(vx))-2 > tile.xrange[0], ceil(max(vx))+2 < tile.xrange[1]] 
            xoutrel = xout-tile.xrange[0]# relative to xrange[0] 
            nxout = xout[1]-xout[0]+1 
            yout = [floor(min(vy))-2 > tile.yrange[0], ceil(max(vy))+2 < tile.yrange[1]] 
            youtrel = yout-tile.yrange[0]# relative to yrange[0] 
            nyout = yout[1]-yout[0]+1 
            rr = rab[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] 
            dd = decb[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] 
            ALLFRAME_ADXYINTERP,head1,rr,dd,xx,yy,nstep=10 
             
            # The x/y position to bilinear need to be in the original system, ~1sec 
            rim = BILINEAR(im1,xx,yy,missing=filestr[i].background) 
            rmask = BILINEAR(mask,xx,yy,missing=0) 
             
            # Contruct final image 
            fim = fltarr(tile.nx,tile.ny)+filestr[i].saturate 
            fim[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] = rim 
            fmask = bytarr(tile.nx,tile.ny) 
            fmask[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] = rmask 
             
            # Contruct the final header 
            fhead = head1 
            # Delete any previous WCS keywords 
            sxdelpar,fhead,['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2','CTYPE1','CTYPE2'] 
            sxdelpar,fhead,['CD1_1','CD1_2','CD2_1','CD2_2'] 
            pvind , = np.where(stregex(strmid(fhead,0,5),'PV[0-9]_+[0-9]',/boolean) == 1,npvind) 
            if npvind > 0 : 
                remove,pvind,fhead 
            # Add the new WCS 
            PUTAST,fhead,tile.ast 
            sxaddpar,fhead,'NAXIS1',tile.nx 
            sxaddpar,fhead,'NAXIS2',tile.ny 
            sxaddpar,fhead,'BPM',filestr[i].resampmask 
            logger.info(filestr[i].resampfile,' ['+str(xout[0],2)+':'+str(xout[1],2)+','+                  str(yout[0],2)+':'+str(yout[1],2)+']' 
            MWRFITS,fim,filestr[i].resampfile,fhead,/create 
            mhead = fhead 
            sxaddpar,mhead,'BITPIX',8 
            sxdelpar,mhead,'BPM' 
            MWRFITS,fmask,filestr[i].resampmask,mhead,/create 
            # this takes about ~37-50 sec for a 2kx4k image. 
            #  now it takes ~5 sec for a 2kx4k image 
 
     
    # --- Pixel based --- 
    'PIXEL': begin 
     
    # Expand the images to sizes that will allow all of the shifts 
    hd1 = PHOTRED_READFILE(fitsfiles[0],/header) 
    nx = sxpar(hd1,'NAXIS1') 
    ny = sxpar(hd1,'NAXIS2') 
    #  left, down, right, up 
    pix_expand = [ abs(floor(min(xshift)) < 0), abs(floor(min(yshift)) < 0),                  ceil(max(xshift)) > 0 , ceil(max(yshift)) > 0 ] 
    outmaskfiles = os.path.dirname(tempfits)+'/'+os.path.basename(tempfits,'.fits')+'.mask.shft.fits' 
    logger.info('Expanding images by [',strjoin(str(pix_expand,2),','),'] pixels' 
    LOADMCH,shiftmch+'.mch',files2,trans2 
    nxf = nx+pix_expand[0]+pix_expand[2] 
    nyf = ny+pix_expand[1]+pix_expand[3] 
    xx1 = lindgen(nx)#replicate(1,ny) 
    yy1 = replicate(1,nx)#lindgen(ny) 
    for i in range(nfiles): 
        # Image 
        tim = PHOTRED_READFILE(tempfits[i],thead) 
        background = np.median(tim) 
        out = trans_coo(xx1[:],yy1[:],reform(trans[i,:])) 
        xx2 = xx1*0. 
        yy2 = yy1*0. 
        xx2[:] = reform(out[0,:]) + pix_expand[0]# shift to expanded grid 
        yy2[:] = reform(out[1,:]) + pix_expand[1]# shift to expanded grid 
        triangulate,xx2,yy2,tr,b# triangulate 
        xout = lindgen(nxf) 
        yout = lindgen(nyf) 
        tim2 = TRIGRID(xx2,yy2,tim, tr, XOUT = xout, YOUT = yout, missing=background) 
        #tim2 = fltarr(nxf,nyf)+background  ; set out of bounds pixels to background 
        #tim2[pix_expand[0]:pix_expand[0]+nx-1,pix_expand[1]:pix_expand[1]+ny-1]=tim 
        thead2 = thead 
        sxaddpar,thead2,'NAXIS1',nxf 
        sxaddpar,thead2,'NAXIS2',nyf 
        sxaddpar,thead2,'CRPIX1',sxpar(thead2,'CRPIX1')+pix_expand[0] 
        sxaddpar,thead2,'CRPIX2',sxpar(thead2,'CRPIX2')+pix_expand[1] 
        MWRFITS,tim2,outfiles[i],thead2,/create 
        # Mask 
        mfile = os.path.basename(tempfits[i],'.fits')+'.mask.fits' 
        mim = PHOTRED_READFILE(mfile,mhead) 
        mim2 = TRIGRID(xx2,yy2,mim, tr, XOUT = xout, YOUT = yout, missing=background) 
        #mim2 = fltarr(nxf,nyf)   ; out of bounds pixels set to 0=bad 
        #mim2[pix_expand[0]:pix_expand[0]+nx-1,pix_expand[1]:pix_expand[1]+ny-1]=mim 
        mhead2 = mhead 
        sxaddpar,mhead2,'NAXIS1',nxf 
        sxaddpar,mhead2,'NAXIS2',nyf 
        sxaddpar,mhead2,'CRPIX1',sxpar(mhead2,'CRPIX1')+pix_expand[0] 
        sxaddpar,mhead2,'CRPIX2',sxpar(mhead2,'CRPIX2')+pix_expand[1] 
        #outmaskfile = FILE_DIRNAME(outfiles[i])+'/'+FILE_BASENAME(outfiles[i],'.fits')+'.mask.shft.fits' 
        MWRFITS,mim2,outmaskfiles[i],mhead2,/create 
 
else: import pdb; pdb.set_trace(),tile.type+' not implemented yet' 
 
 
 
# Creating new MCH file for the combined file 
if not keyword_set(fake): 
print('Deriving new transformation equations for the resampled coordinate system' 
for i in range(nfiles): 
     
    # Convert X/Y of this system into the combined reference frame 
    #  The pixel values are 1-indexed like DAOPHOT uses. 
    ngridbin = 50 
    nxgrid = filestr[i].nx / ngridbin 
    nygrid = filestr[i].ny / ngridbin 
    xgrid = (lindgen(nxgrid)*ngridbin+1)#replicate(1,nygrid) 
    ygrid = replicate(1,nxgrid)#(lindgen(nygrid)*ngridbin+1) 
    HEAD_XYAD,(*filestr[i].head),xgrid-1,ygrid-1,ragrid,decgrid,/deg 
    HEAD_ADXY,tile.head,ragrid,decgrid,refxgrid,refygrid,/deg 
    refxgrid += 1# convert 0-indexed to 1-indexed 
    refygrid += 1 
     
    # Now fit the transformation 
    xdiff = refxgrid-xgrid 
    ydiff = refygrid-ygrid 
    xmed = np.median([xdiff],/even) 
    ymed = np.median([ydiff],/even) 
    # Fit rotation with linear fits if enough points 
    coef1 = robust_poly_fitq(ygrid,xdiff,1)# fit rotation term 
    coef1b = dln_poly_fit(ygrid,xdiff,1,measure_errors=xdiff*0+0.1,sigma=coef1err,/bootstrap) 
    coef2 = robust_poly_fitq(xgrid,ydiff,1)# fit rotation term 
    coef2b = dln_poly_fit(xgrid,ydiff,1,measure_errors=ydiff*0+0.1,sigma=coef2err,/bootstrap) 
    #theta = mean([-coef1[1],coef2[1]]) 
    WMEANERR,[-coef1[1],coef2[1]],[coef1err[1],coef2err[1]],theta,thetaerr 
     
    # [xoff, yoff, cos(th), sin(th), -sin(th), cos(th)] 
    #trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0] 
    trans = [xmed, ymed, 1.0-theta**2, theta, -theta, 1.0-theta**2] 
    # Adjust Xoff, Yoff with this transformation 
    xyout = trans_coo(xgrid,ygrid,trans) 
    trans[0] += np.median([refxgrid - xyout[0,:]],/even) 
    trans[1] += np.median([refygrid - xyout[1,:]],/even) 
     
    # Fit full six parameters if there are enough stars 
    fa = {x1:(refxgrid)(*),y1:(refygrid)(*),x2:(xgrid)(*),y2:(ygrid)(*)} 
    initpar = trans 
    fpar = MPFIT('trans_coo_dev',initpar,functargs=fa, perror=perror, niter=iter, status=status,                 bestnorm=chisq,:f=dof, autoderivative=1, /quiet) 
    trans = fpar 
     
    diff = trans_coo_dev(fpar,x1=refxgrid,y1=refygrid,x2=xgrid,y2=ygrid) 
    rms = sqrt(np.mean(diff**2.)) 
    filestr[i].resamptrans = trans 
    filestr[i].resamptransrms = rms 
     
    # The output is: 
    # filename, xshift, yshift, 4 trans, mag offset, magoff sigma 
    format = '(A2,A-30,A1,2A10,4A12,F9.3,F8.4)' 
    # In daomaster.f the translations are 10 digits with at most 4 
    # decimal places (with a leading space), the transformation 
    # coefficients are 12 digits with at most 9 decimal places. 
    # Need a leading space to separate the numbers. 
    strans = ' '+[str(string(trans[0:1],format='(F30.4)'),2),                 str(string(trans[2:5],format='(F30.9)'),2)] 
    newline = STRING("'",filestr[i].catfile,"'", strans, filestr[i].magoff[0], rms, format=format) 
    PUSH,mchfinal,newline 
     
    # Printing the transformation 
    logger.info(format='(A-20,2A10,4A12,F9.3,F8.4)',filestr[i].catfile,strans,filestr[i].magoff[0],rms 
# Write to the new MCH file 
combmch = mchbase+'_comb.mch' 
WRITELINE,combmch,mchfinal 
 
# FAKE, use existing one 
else: 
combmch = mchbase+'_comb.mch' 
# don't need to load the information 
 
 
############################################ 
# STEP 5: COMBINE IMAGES 
logger.info('-------------------' 
logger.info('STEP 5: IMCOMBINE' 
logger.info('-------------------' 
 
# The imcombine input file 
resampfile = mchbase+'.resamp' 
WRITELINE,resampfile,filestr.resampfile 
 
# SCALE the images for combining 
#------------------------------- 
if not keyword_set(nocmbimscale): 
 
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
os.remove(combfile,/allow 
os.remove(mchbase+'_comb.bpm.pl',/allow 
IRAF_IMCOMBINE,'@'+resampfile,combfile,combine='average',reject='avsigclip',                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',                 irafdir=irafdir,error=imcombineerror2,scale='@'+scalefile,zero='@'+zerofile,                 masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm' 
 
if len(imcombineerror2) != 0: 
    logger.info('ERROR in IRAF_IMCOMBINE' 
    logger.info(imcombineerror2 
    error = imcombineerror2 
    return 
 
# Convert BPM mask from PL to FITS 
os.remove(mchbase+'_comb.bpm.fits',/allow 
undefine,lines 
cd,current=curdir 
push,lines,'print("")'# first line will be ignored 
push,lines,'cd '+curdir 
push,lines,'imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits' 
push,lines,'logout' 
tempfile = mktemp('tiraf') 
WRITELINE,tempfile,lines 
IRAF_RUN,tempfile,irafdir,silent=silent,out=out,error=error 
 
# Delete temporary scripts and PL file 
os.remove([tempfile,mchbase+'_comb.bpm.pl'],/allow 
 
 
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
combim = PHOTRED_READFILE(combfile,combhead) 
badmask = PHOTRED_READFILE(mchbase+'_comb.bpm.fits',maskhead)# 0-good, 1-bad 
 
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
    rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits') 
#  the "scales" array here is actually 1/scales used by IMCOMBINE. 
rdnoise = sqrt(np.sum((weights*rdnoisearr/scales)**2)) 
rdnoise = rdnoise > 0.01# must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey)# get keyword 
sxaddpar,combhead,rdnoisekey,rdnoise 
 
# Fix the sky 
# DAOPHOT FIND computes the random error per pixel in ADU as 
# noise = sqrt( sky level/gain + rdnoise^2) 
# So it assumes that the noise in the background is sqrt(sky/gain) 
# in ADU.  We need to set the sky level so this is correct. 
# The final noise should be 
# final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
# So the final sky level should be 
# final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
gain = PHOTRED_GETGAIN(combfile,keyword=gainkey) 
comb_sky = np.sum((weights*sqrt((sky>0)/gain)/scales)**2)*gain 
# the "scales" array here is actually 1/scale 
combim += float(comb_sky)# keep it float 
 
 
# set the maximum to a "reasonable" level 
# Rescale the image and increase the gain 
if max(combim) > 50000: 
    rescale = 50000./max(combim) 
    combim = combim*rescale 
    sxaddpar,combhead,gainkey,gain/rescale 
    # rdnoise does NOT get modified since it's in electrons 
    # we just need to modify the gain which takes you from ADUs to electrons 
 
maskdatalevel = max(combim) + 10000# set "bad" data level above the highest "good" value 
combim2 = combim*(1-badmask) + maskdatalevel*badmask# set bad pixels to maskdatalevel 
sxaddpar,combhead,'SATURATE',maskdatalevel 
MWRFITS,combim2,combfile,combhead,/create# fits_write can create an empty PDU 
 
# Create the weight map for Sextractor using the BPM output by IMCOMBINE 
#  bad only if bad in ALL images 
weightmap = -2.0*float(badmask == 1) + 1.0 
combweightfile = mchbase+'_comb.mask.fits' 
MWRFITS,weightmap,combweightfile,whead,/create 
 
# NO SCALING of the images for combining 
#--------------------------------------- 
else: 
 
combfile = mchbase+'_comb.fits' 
os.remove(combfile,/allow 
IRAF_IMCOMBINE,'@'+resampfile,combfile,combine='average',reject='avsigclip',                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',                 irafdir=irafdir,error=imcombineerror2 
 
if len(imcombineerror2) != 0: 
    logger.info('ERROR in IRAF_IMCOMBINE' 
    logger.info(imcombineerror2 
    error = imcombineerror2 
    return 
 
# Fix the rdnoise and background/sky level and saturate 
#  the bad pixels for DAOPHOT 
#------------------------------------------------------ 
 
# See the explanations for all these steps above!! 
 
# Load the IMCOMBINE output combined file and BPM 
combim = PHOTRED_READFILE(combfile,combhead) 
badmask = PHOTRED_READFILE(mchbase+'_comb.bpm.fits',maskhead)# 0-good, 1-bad 
 
# Fix the rdnoise 
# The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2)) 
rdnoisearr = fltarr(nfiles) 
for i in range(nfiles): 
    rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits') 
rdnoise = sqrt(np.sum((weights*rdnoisearr)**2)) 
rdnoise = rdnoise > 0.01# must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey)# get keyword 
sxaddpar,combhead,rdnoisekey,rdnoise 
 
# Fix the sky 
# So the final sky level should be 
# final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
gain = PHOTRED_GETGAIN(combfile,keyword=gainkey) 
comb_sky = np.sum((weights*sqrt(sky/gain))**2)*gain 
# the "scales" array here is actually 1/scale 
combim += comb_sky 
 
 
# set the maximum to a "reasonable" level 
# Rescale the image and increase the gain 
if max(combim) > 50000: 
    rescale = 50000./max(combim) 
    combim = combim*rescale 
    sxaddpar,combhead,gainkey,gain/rescale 
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
FITS_WRITE,combweightfile,weightmap,whead 
 
#--------------------------------------------- 
# SATURATE BAD pixels in the COMBINED IMAGE 
# DAOPHOT needs to have the bad pixels "saturated", 
# SExtractor will know which pixels are bad from the "weight" map. 
# 
# We could skip the fiximage.pro step but we still need the 
# individual bpm masks and setting the bad pixels to the background 
# probably helps in the IMALIGN/IMCOMBINE steps. 
logger.info('' 
logger.info('"Saturating" bad pixels in the COMBINED image' 
logger.info('' 
 
 
badmask = float(weightmap < 0.5) 
maskdatalevel = max(combim) + 10000# set "bad" data level above the highest "good" value 
combim2 = combim*(1.0-badmask) + maskdatalevel*badmask# set bad pixels to 100,000 
sxaddpar,combhead,'SATURATE',maskdatalevel 
FITS_WRITE,combfile,combim2,combhead 
 
# no scaling of images for combining 
 
# Add TILETYPE to the combined image 
combhead = PHOTRED_READFILE(combfile,/header) 
sxaddpar,combhead,'AFTILTYP',tile.type 
MODFITS,combfile,0,combhead 
 
# Delete the resampled images 
os.remove(filestr.resampfile,/allow,/quiet 
os.remove(filestr.resampmask,/allow,/quiet 
 
# Make the common source file for the combined image 
#--------------------------------------------------- 
if keyword_set(usecmn): 
print('Combining COMMON SOURCE files for the combined image.' 
# Loop through the files and convert to coordinates to the comined file 
undefine,allcmn 
for i in range(nfiles): 
    cmnfile1 = os.path.basename(filestr[i].fitsfile,'.fits')+'.cmn.lst' 
    if os.path.exists(cmnfile1) == 1: 
        cmn1 = IMPORTASCII(cmnfile1,fieldnames=['id','x','y','mag','err','sky','skysig','sharp','round','round2'],                         skipline=3,/silent) 
        READLINE,cmnfile1,cmnlines1 
        coohead1 = cmnlines1[0:1] 
        ncmn1 = len(cmn1) 
        # Get coordinates on the resampled/combined image grid 
        out = trans_coo(cmn1.x,cmn1.y,filestr[i].resamptrans) 
        newx = reform(out[0,:]) 
        newy = reform(out[1,:]) 
        cmn1.x = newx 
        cmn1.y = newy 
        if len(allcmn) == 0: 
            allcmn = cmn1 
        else: 
            # Remove any duplicates 
            SRCMATCH,allcmn.x,allcmn.y,cmn1.x,cmn1.y,2.0,ind1,ind2,count=nmatch 
            if nmatch > 0: 
                if nmatch < ncmn1: 
                    remove,ind2,cmn1 
                else: 
                    undefine,cmn1 
            if len(cmn1) > 0 : 
                push,allcmn,cmn1 
if len(allcmn) > 0: 
    WRITECOL,mchbase+'_comb.cmn.lst',allcmn.id,allcmn.x,allcmn.y,allcmn.mag,allcmn.err,allcmn.sky,allcmn.skysig,             allcmn.sharp,allcmn.round,allcmn.round2,fmt='(I7,2F9.2,3F9.3,F9.2,3F9.3)' 
    WRITELINE,mchbase+'_comb.cmn.lst',[coohead1,''],/prepend# prepend the COO header 
 else logger.info('No common file to combine' 
     
    # Don't use common file 
else: 
    # Make sure it doesn't exist otherwise it will be used 
    os.remove(mchbase+'_comb.cmn.lst',/allow 
 
    BOMB: 


 
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
 
    Example
    -------

    maskdatalevel = combine_orig('ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/',finditer=2)
 
 
    By D.Nidever   February 2008 
    Automation of steps and scripts by J.Ostheimer and Rachael Beaton 
    Translated to Python by D. Nidever, April 2022
    """
 
              
    global photred,setup 

     
    # Logfile 
    if keyword_set(logfile): 
        logf=logfile 
    else: 
        logf=-1 
     
    # Getting scripts directory and iraf directory 
    nsetup = len(setup) 
    if nsetup > 0: 
        scriptsdir = READPAR(setup,'SCRIPTSDIR') 
        irafdir = READPAR(setup,'IRAFDIR') 
     
     
    # No irafdir 
    if len(scriptsdir) == 0: 
        logger.info('SCRIPTSDIR NOT INPUT' 
        error = 'SCRIPTSDIR NOT INPUT' 
        return 
     
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if iraf.check(irafdir=irafdir)==False:
        raise ValueError('IRAF TEST FAILED')
     
    # No scriptsdir 
    if len(scriptsdir) == 0: 
        logger.info('SCRIPTSDIR NOT INPUT' 
        error = 'SCRIPTSDIR NOT INPUT' 
        return 
    # Check if the scripts exist in the current directory 
    scripts = ['getpsf.sh','photo.opt','apcor.opt','lstfilter.py','goodpsf.pro','allframe.opt',
               'default.sex','default.param','default.nnw','default.conv'] 
    nscripts = len(scripts) 
    # Loop through the scripts 
    for i in range(nscripts): 
        info = FILE_INFO(scriptsdir+'/'+scripts[i]) 
        curinfo = FILE_INFO(scripts[i]) 
         
        # No file 
        if info.exists == 0 or info.size == 0: 
            logger.info(scriptsdir+'/'+scripts[i],' NOT FOUND or EMPTY' 
            error = scriptsdir+'/'+scripts[i]+' NOT FOUND or EMPTY' 
            return 
         
        # Check if the two files are the same size, if not copy it 
        if info.size != curinfo.size: 
            FILE_COPY,info.name,curinfo.name,/overwrite 
     
     
    logger.info('Combining images in ',filename) 
     
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
    files,trans,magoffset = io.loadmch(mchfile)
     
    # Check that the fits, als, opt, and psf files exist 
    nfiles = len(files) 
    for i in range(nfiles): 
        base = os.path.splitext(os.path.basename(files[i]))[0]
         
        # Checking FITS file 
        fitstest = os.path.exists(base+'.fits') 
        if fitstest == 0: 
            logger.info(base+'.fits NOT FOUND' 
            return 
         
        # Checking OPT file 
        opttest = os.path.exists(base+'.opt') 
        if opttest == 0: 
            logger.info(base+'.opt NOT FOUND' 
            return 
         
        # Checking ALS.OPT file 
        alsopttest = os.path.exists(base+'.als.opt') 
        if alsopttest == 0: 
            logger.info(base+'.als.opt NOT FOUND' 
            return 
         
        # Checking AP file 
        aptest = os.path.exists(base+'.ap') 
        if aptest == 0: 
            logger.info(base+'.ap NOT FOUND' 
            return 
         
        # Checking ALS file 
        alstest = os.path.exists(base+'.als') 
        if alstest == 0: 
            logger.info(base+'.als NOT FOUND' 
            return 
         
     
    # FAKE, check that we have all the files that we need 
    if keyword_set(fake): 
        # weights, scale, zero, comb_psf, _shift.mch 
        chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_shift.mch'] 
        bdfiles , = np.where(os.path.exists(chkfiles) == 0,nbdfiles) 
        if nbdfiles > 0: 
            error = 'FAKE.  Some necessary files not found. '+strjoin(chkfiles[bdfiles],' ') 
            logger.info(error 
            return 
     
     
     
    ############################################ 
    # STEP 1: IMALIGN PREP 
     
    logger.info('')
    logger.info('Step A: Getting Weights')
    logger.info('-----------------------')
     
    #----------------------------------- 
    # Computs Weights 
    if not keyword_set(fake): 
        ALLFRAME_GETWEIGHTS,mchfile,weights,scales,sky#,raw2=raw2 
        invscales = 1.0/scales 
        bdscale , = np.where(scales < 1e-5 or invscales > 900,nbdscale) 
        if nbdscale > 0: 
            scales[bdscale] = 1.0 
            invscales[bdscale] = 1.0 
            weights[bdscale] = 0.0 
        weightfile = mchbase+'.weights' 
        WRITECOL,weightfile,weights,fmt='(F10.6)' 
        scalefile = mchbase+'.scale' 
        WRITECOL,scalefile,invscales,fmt='(F10.6)'# want to scale it UP 
        zerofile = mchbase+'.zero' 
        WRITECOL,zerofile,-sky,fmt='(F12.4)'# want to remove the background, set to 1st frame 
         
    # FAKE, use existing ones 
    else: 
        weightfile = mchbase+'.weights' 
        scalefile = mchbase+'.scale' 
        zerofile = mchbase+'.zero' 
        READCOL,weightfile,weights,format='F',/silent 
        READCOL,scalefile,invscales,format='F',/silent 
        scales = 1.0/invscales 
        READCOL,zerofile,sky,format='F',/silent 
        sky = -sky 
     
     
    #--------------------------------------- 
    # Get X/Y translations using DAOMASTER 
    #  NO ROTATION ONLY SHIFTS 
    #  Need these shifts for IMSHIFT 
    shiftmch = mchbase+'_shift'
    if fake==False:
        logger.info('Measuring X/Y shifts')
        FILE_COPY,mchbase+'.mch',shiftmch+'.mch',/overwrite,/allow 
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
        tempscript = MKTEMP('daomaster')# absolute filename 
        dln.writelines(tempscript,cmdlines)
        os.chmod(tempscript,0o755)
        # Run DAOMASTER 
        cmd2 = tempscript+' '+shiftmch 
        out2 = subprocess.check_output(cmd2)
        # Remove temporary DAOMASTER script 
        os.remove(tempscript,/allow_non 
    files2,trans2,magoffset2 = io.loadmch(shiftmch+'.mch')
     
    xshift = reform(trans2[:,0]) 
    yshift = reform(trans2[:,1]) 
    xyshifts = [[xshift],[yshift]] 
    logger.info('Image shifts' 
    for i in range(nfiles): 
        logger.info(files[i],xshift[i],yshift[i])
     
     
    #----------------------------------- 
    # Create imalign prep files 
    # This is done by preimalign_k.sh 
    # Need an input list of fits files 
    # Need an output list of fits files 
    # Shift file 
    base = os.path.splitext(os.path.basename(files))[0]
    fitsfiles = base+'.fits' 
    outfiles = base+'.shft.fits' 
    infile = mchbase+'.inlist' 
    outfile = mchbase+'.outlist' 
    #WRITELINE,infile,fitsfiles   ; this is done below now with the temp files 
    dln.writelines(outfile,outfiles)  
    # Remove outfiles 
    os.remove(outfiles,/allow 
    # shift list 
    shiftfile = mchbase+'.shift' 
    WRITECOL,shiftfile,xshift,yshift,fmt='(2F15.4)' 
     
     
    # Make temporary files for bad pixel fixing and combining 
    #  FIXIMAGES doesn't work properly on the shifted images 
    #  because the interpolation can bring the bad pixel values down 
    #  below the saturation threshold, and we don't want to touch 
    #  the original images. 
    tempfits = base+'.temp.fits' 
    FILE_COPY,fitsfiles,tempfits,/overwrite,/allow 
    dln.writelines(infile,tempfits)
     
     
    ############################################ 
    # STEP B: FIX BAD PIXELS 
    logger.info('')
    logger.info('Step B: Fixing bad pixels')
    logger.info('-------------------------')
    FIXIMAGES,'@'+infile,satlevel=satlevel  #6e4 
    # This also makes the FILE.mask.fits files for each image 
     
    # Find the maximum saturation level 
    satlevelarr = fltarr(nfiles) 
    for i in range(nfiles): 
        #head = headfits(base[i]+'.fits') 
        im = PHOTRED_READFILE(base[i]+'.fits',head) 
        saturate = sxpar(head,'SATURATE',count=nsaturate,/silent) 
        if nsaturate == 0 : 
            saturate=max(im)-1000. 
        satlevelarr[i] = saturate 
    maxsatlevel = max(satlevelarr) 
     
     
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
    if len(imshifterror) != 0: 
        logger.info('ERROR in IRAF_IMSHIFT')
        logger.info(imshifterror)
        error = imshifterror)
        return 
     
    # Trim the images 
    if trimcomb: 
        # Calculate the trim section 
        hd = io.readfile(reffile,header=True) 
        xsize = lonarr(nfiles)+hd['NAXIS1']
        ysize = lonarr(nfiles)+hd['NAXIS2']
        IA_TRIM,xshift,yshift,xsize,ysize,trimsection 
        xoff = trimsection[0]-1 
        yoff = trimsection[2]-1 
         
        # Trim the shifted images 
        logger.info('Trimming the shifted images')
        xstart = trimsection[0]-1 
        ximport pdb; pdb.set_trace() = trimsection[1]-1 
        ystart = trimsection[2]-1 
        yimport pdb; pdb.set_trace() = trimsection[3]-1 
         
        for i in range(nfiles): 
            im = PHOTRED_READFILE(outfiles[i],head) 
            newim = im[xstart:ximport pdb; pdb.set_trace(),ystart:yimport pdb; pdb.set_trace()] 
            MWRFITS,newim,outfiles[i],head,/create,/silent 
        # could also use IRAF_IMCOPY here instead 
         
         
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
    out = run(imshiftscript,irafdir,silent=True)
    if len(iraferror) != 0: 
        logger.info('ERROR in running IMSHIFT with IRAF_RUN')
        logger.info(iraferror)
        error = iraferror 
        return 
     
    # Trim 
    if keyword_set(trimcomb): 
        logger.info('Trimming masks' 
        xstart = trimsection[0]-1# should be same as xoff 
        ximport pdb; pdb.set_trace() = trimsection[1]-1 
        ystart = trimsection[2]-1# should be same as yoff 
        yimport pdb; pdb.set_trace() = trimsection[3]-1 
         
        for i in range(nfiles): 
            im = PHOTRED_READFILE(outmaskfiles[i],head) 
            sz = size(im) 
            newim = im[xstart:ximport pdb; pdb.set_trace(),ystart:yimport pdb; pdb.set_trace()] 
            # Add LTV1/LTV2 to the header 
            #  these are IRAF keywords to convert from logical to physical coords 
            ltv1 = sxpar(head,'LTV1',/silent)# 0 if not found 
            ltv2 = sxpar(head,'LTV2',/silent)# 0 if not found 
            sxaddpar,head,'LTV1',ltv1-xstart 
            sxaddpar,head,'LTV2',ltv2-ystart 
            MWRFITS,newim,outmaskfiles[i],head,/create,/silent 
     
    # Combining masks 
    logger.info('Combining masks' 
    undefine,bpm 
    for i in range(nfiles): 
        im = PHOTRED_READFILE(outmaskfiles[i],head) 
        if i == 0: 
            bpm = im 
            whead = head 
        else: 
            bpm = bpm+im 
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
            head = PHOTRED_READFILE(outfiles[i],/header) 
            sxaddpar,head,'BPM',outmaskfiles[i] 
            MODFITS,outfiles[i],0,head 
         
        # Combine the frames WITH scaling/offset/masking, for the bright stars 
        #logger.info('Creating SCALED image' 
        combfile = mchbase+'_comb.fits' 
        os.remove(combfile,/allow 
        os.remove(mchbase+'_comb.bpm.pl',/allow 
        IRAF_IMCOMBINE,'@'+outfile,combfile,combine='average',reject='avsigclip',                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',                 irafdir=irafdir,error=imcombineerror2,scale='@'+scalefile,zero='@'+zerofile,                 masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm' 
         
        if len(imcombineerror2) != 0: 
            logger.info('ERROR in IRAF_IMCOMBINE' 
            logger.info(imcombineerror2 
            error = imcombineerror2 
            return 
         
        # Convert BPM mask from PL to FITS 
        os.remove(mchbase+'_comb.bpm.fits',/allow 
        undefine,lines 
        cd,current=curdir 
        push,lines,'print("")'# first line will be ignored 
        push,lines,'cd '+curdir 
        push,lines,'imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits' 
        push,lines,'logout' 
        tempfile = mktemp('tiraf') 
        WRITELINE,tempfile,lines 
        IRAF_RUN,tempfile,irafdir,silent=silent,out=out,error=error 
         
        # Delete temporary scripts and PL file 
        os.remove([tempfile,mchbase+'_comb.bpm.pl'],/allow 
         
         
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
        combim = PHOTRED_READFILE(combfile,combhead) 
        badmask = PHOTRED_READFILE(mchbase+'_comb.bpm.fits',maskhead)# 0-good, 1-bad 
         
         
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
            rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits') 
        #  the "scales" array here is actually 1/scales used by IMCOMBINE. 
        rdnoise = sqrt(np.sum((weights*rdnoisearr/scales)**2)) 
        rdnoise = rdnoise > 0.01# must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
        dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey)# get keyword 
        sxaddpar,combhead,rdnoisekey,rdnoise 
         
        # Fix the sky 
        # DAOPHOT FIND computes the random error per pixel in ADU as 
        # noise = sqrt( sky level/gain + rdnoise^2) 
        # So it assumes that the noise in the background is sqrt(sky/gain) 
        # in ADU.  We need to set the sky level so this is correct. 
        # The final noise should be 
        # final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
        # So the final sky level should be 
        # final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
        gain = PHOTRED_GETGAIN(combfile,keyword=gainkey) 
        comb_sky = np.sum((weights*sqrt((sky>0)/gain)/scales)**2)*gain 
        # the "scales" array here is actually 1/scale 
        combim += float(comb_sky)# keep it float 
         
         
        # set the maximum to a "reasonable" level 
        # Rescale the image and increase the gain 
        if max(combim) > 50000: 
            rescale = 50000./max(combim) 
            combim = combim*rescale 
            sxaddpar,combhead,gainkey,gain/rescale 
            # rdnoise does NOT get modified since it's in electrons 
            # we just need to modify the gain which takes you from ADUs to electrons 
         
        maskdatalevel = max(combim) + 10000# set "bad" data level above the highest "good" value 
        combim2 = combim*(1-badmask) + maskdatalevel*badmask# set bad pixels to maskdatalevel 
        MWRFITS,combim2,combfile,combhead,/create# fits_write can create an empty PDU 
         
        # Create the weight map for Sextractor using the BPM output by IMCOMBINE 
        #  bad only if bad in ALL images 
        weightmap = -2.0*float(badmask == 1) + 1.0 
        combweightfile = mchbase+'_comb.mask.fits' 
        MWRFITS,weightmap,combweightfile,whead,/create 
         
    # NO SCALING of the images for combining 
    #--------------------------------------- 
    else: 
         
        combfile = mchbase+'_comb.fits' 
        os.remove(combfile,/allow 
        IRAF_IMCOMBINE,'@'+outfile,combfile,combine='average',reject='avsigclip',                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',                 irafdir=irafdir,error=imcombineerror2 
         
        if len(imcombineerror2) != 0: 
            logger.info('ERROR in IRAF_IMCOMBINE' 
            logger.info(imcombineerror2 
            error = imcombineerror2 
            return 
         
        # Fix the rdnoise and background/sky level and saturate 
        #  the bad pixels for DAOPHOT 
        #------------------------------------------------------ 
         
        # See the explanations for all these steps above!! 
         
        # Load the IMCOMBINE output combined file and BPM 
        combim = PHOTRED_READFILE(combfile,combhead) 
        badmask = PHOTRED_READFILE(mchbase+'_comb.bpm.fits',maskhead)# 0-good, 1-bad 
         
        # Fix the rdnoise 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2)) 
        rdnoisearr = fltarr(nfiles) 
        for i in range(nfiles): 
            rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits') 
        rdnoise = sqrt(np.sum((weights*rdnoisearr)**2)) 
        rdnoise = rdnoise > 0.01# must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
        dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey)# get keyword 
        sxaddpar,combhead,rdnoisekey,rdnoise 
         
        # Fix the sky 
        # So the final sky level should be 
        # final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
        gain = PHOTRED_GETGAIN(combfile,keyword=gainkey) 
        comb_sky = np.sum((weights*sqrt(sky/gain))**2)*gain 
        # the "scales" array here is actually 1/scale 
        combim += comb_sky 
         
         
        # set the maximum to a "reasonable" level 
        # Rescale the image and increase the gain 
        if max(combim) > 50000: 
            rescale = 50000./max(combim) 
            combim = combim*rescale 
            sxaddpar,combhead,gainkey,gain/rescale 
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
        FITS_WRITE,combweightfile,weightmap,whead 
         
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
        FITS_WRITE,combfile,combim2,combhead 
         
# no scaling of images for combining 
     
    # Add TILETYPE to the combined image 
    combhead = PHOTRED_READFILE(combfile,/header) 
    sxaddpar,combhead,'AFTILTYP','ORIG' 
    MODFITS,combfile,0,combhead 
     
    # Delete the shifted images
    shiftedfiles = dln.readlines(outfile)
    os.remove(shiftedfiles,/allow,/quiet 
     
    # Delete mask files 
    os.remove([maskfiles,outmaskfiles],/allow 
    os.remove([maskinfile,maskoutfile,maskshiftsfile,imshiftscript],/allow 
     
    # Copy the original MCH file to COMB.MCH 
    FILE_COPY,file,mchdir+'/'+mchbase+'.comb.mch',/allow,/over 
     
    # Using CMN.LST of reference frame if it exists 
    if os.path.exists(mchbase+'.cmn.lst') and keyword_set(usecmn): 
        logger.info('Using reference image COMMON SOURCE file')
        FILE_COPY,mchbase+'.cmn.lst',mchbase+'_comb.cmn.lst',/over,/allow 
     
    # CD back to the original directory 
    cd,curdir 
    os.chdir(curdir)
              
    BOMB: 
     
 
