#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning
import warnings
from dlnpyutils import utils as dln
from . import utils,imfwhm,io

# Ignore these warnings
warnings.simplefilter('ignore', category=AstropyWarning)
#warnings.filterwarnings(action="ignore", message=r'FITSFixedWarning:*')


def mkopt(inpfiles,hilimit=6.4e4,va=2,fitradius_fwhm=None,
          inp_fwhm=None,verbose=False):
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
    nfiles = len(files)
     
    # Not enough inputs 
    if nfiles == 0: 
        raise ValueError('No files') 
     
    # More than one name input 
    if nfiles > 1: 
        fwhm = fltarr(nfiles) 
        for i in range(nfiles): 
            fwhm1 = mkopt(files[i],hilimit=hilimit,va=va,fitradius_fwhm=fitradius_fwhm,verbose=verbose)
            fwhm[i] = fwhm1
            if verbose:
                print('')
        return fwhm
     
    filename = str(files[0]).strip()
    im = None
     
    # Default settings
    VA = np.maximum(va,2)
    if fitradius_fwhm is None or fitradius_fwhm <= 0:
        fitradius_fwhm = 1.0
    fitradius_fwhm = np.maximum(fitradius_fwhm,1.0)
         
         
    # Processing ONE file 
    #-------------------- 
    if os.path.exists(filename)==False:
        raise ValueError(filename+' NOT FOUND')

    if verbose:
        print('Running MKOPT on '+str(filename))
         
    # Get the base
    base = utils.fitsext(filename,basename=True)
    if filename.endswith('.fz'):
        fpack = True
    else:
        fpack = False
    fdir = os.path.dirname(filename) 
    if fdir=='':
        fdir = '.'
    
             
    # Get the FITS header 
    if fpack:
        head = io.readfile(filename,exten=1,header=True)
    else: 
        head = io.readfile(filename,header=True)
         
    # We need GAIN, READNOISE, FWHM, and Hi-limit 
    #-------------------------------------------- 
         
    # Getting GAIN
    gain,gainkey = io.getgain(filename)
    # Getting READNOISE 
    rdnoise,rdnoisekey = io.getrdnoise(filename)
         
    # Run IMFWHM to get the FWHM
    if inp_fwhm is None:
        fwhm,ellip,gtab,ptab = imfwhm.imfwhm(filename)
        im = io.readfile(filename)
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
            return None
             
    else:
        fwhm = inp_fwhm 
             
    # Load the image
    if im is None:
        im,head = io.readfile(filename)
             
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
    HI = hi 
    
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
            form = '%-5s%8.2f\n'
            if anotarr[j] == 'HI = ':
                form = '%5s,%8d'
            f.write(form % (anotarr[j],outarr[j]))
        f.write('\n')  # we need an extra blank line at end             

    # Writing the ALLSTAR parameter file 
    #----------------------------------- 
             
    # FI    : fitting radius 
    # IS    : inner sky radius
    # OS    : outer sky radius
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
            form = '%-5s%8.2f\n'
            f.write(form % (anotarr2[j],outarr2[j]))
        f.write('\n')  # we need an extra blank line at end
             
    # Verbose output 
    if verbose:
        print('Created ',fdir+'/'+base+'.opt')
        print('Created ',fdir+'/'+base+'.als.opt') 

