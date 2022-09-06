#!/usr/bin/env python

import os
import time
import numpy as np
import subprocess 
import tempfile
import shutil
from glob import glob
from astropy.table import Table
from dlnpyutils import utils as dln
from . import io,utils

def mksexconfig(filename,configfile=None,catfile=None,flagfile=None,
                paramfile=None,wtfile=None,cattype=None):
    """
    This program creates a SExtractor configuration file 
    for an inputs FITS file. 
 
    Parameters
    ----------
    filename : str
       Filename of the FITS image 
    configfile: str, optional
       Filename for the SExtractor configuration file. 
         This is optional.  If not input then the 
         name will be FITSFILEBASE.sex 
    catfile : str, optional
       The output catalog name.  If not input then 
          FITSFILEBASE.cat will be used 
    flagfile : str, optional
       The flag/mask file to use. 
    paramfile : str, optional
       The output parameter/columns file to use. 
    wtfile : str, optional
       The weight image file to use. 
    cattype : str, optional
       The output type of the catalog: ASCII or FITS_LDAC 
 
    Returns
    -------
    The SExtractor config file will be written to CONFIGFILE. 
 
    Example
    -------
    mksexconfig('F1-89230023_05.fits','F1-89230023_05.sex')
 
    By D. Nidever    January 2019 
    Translated to python by D. Nidever  Sep 2022
    """

    base,_ = os.path.splitext(os.path.basename(filename))
     
    # Make sure file exists 
    if os.path.exists(filename) == 0: 
        raise ValueError(filename+' NOT FOUND')
    if os.path.exists(base+'.opt') == 0: 
        raise ValueError(base+'.opt NOT FOUND')
     
    # Read the .opt file 
    optlines = dln.readlines(base+'.opt')
    opt = io.readopt(base+'.opt')
    satlevel = opt['HI']
    gain = opt['GA']
    fwhm = opt['FW']

    # Get pixel scale 
    scale = io.getpixscale(filename)
    if scale is None or scale>90:
        scale = 0.5 

    # Make customized SEXTRACTOR file 
    #-------------------------------- 
    if configfile is None:
        configfile = base+'.sex' 
    if catfile is None:
        catfile = base+'.cat' 
    if os.path.exists('default.sex') == False:
        out = subprocess.run(['sex','-d'],shell=False)        
        dln.writelines('default.sex',out)
    shutil.copyfile('default.sex',configfile)
    #FILE_COPY,'default.sex',configfile,/overwrite,/allow 
    sexlines = dln.readlines(configfile)
    sexlines2 = np.char.array(sexlines).copy()
    # CATALOG_NAME 
    g, = np.where(sexlines2.find('CATALOG_NAME') > -1) 
    sexlines2[g[0]] = 'CATALOG_NAME    '+catfile+'  # name of the output catalog' 
    # SATUR_LEVEL 
    g, = np.where(sexlines2.find('SATUR_LEVEL') > -1)
    sexlines2[g[0]] = 'SATUR_LEVEL     '+str(satlevel)+'  # level (in ADUs) at which arises saturation' 
    # GAIN 
    g, = np.where(sexlines2.find('GAIN') > -1)
    sexlines2[g[0]] = 'GAIN            '+str(gain)+'  # detector gain in e-/ADU.' 
    # PIXEL_SCALE 
    g, = np.where(sexlines2.find('PIXEL_SCALE') > -1)
    sexlines2[g[0]] = 'PIXEL_SCALE     '+str(scale)+'  # size of pixel in arcsec (0=use FITS WCS info).' 
    # SEEING_FWHM 
    g, = np.where(sexlines2.find('SEEING_FWHM') > -1)
    fwhmas = float(fwhm)*float(scale) 
    sexlines2[g[0]] = 'SEEING_FWHM     '+str(fwhmas)+'  # stellar FWHM in arcsec' 
    # DETECT_MINAREA 
    g, = np.where(sexlines2.find('DETECT_MINAREA') > -1)
    sexlines2[g[0]] = 'DETECT_MINAREA     4  # minimum number of pixels above threshold' 
    # Catalog type 
    if cattype is not None:
        g, = np.where(sexlines2.find('CATALOG_TYPE') > -1)
        if len(g) > 0: 
            sexlines2[g[0]] = 'CATALOG_TYPE    '+cattype 
        else: 
            sexlines2 = np.append(sexlines2,'CATALOG_TYPE    '+cattype)
    # Parameters file 
    if paramfile is not None:
        g, = np.where(sexlines2.find('PARAMETERS_NAME') > -1)
        if len(g) > 0: 
            sexlines2[g[0]] = 'PARAMETERS_NAME  '+paramfile 
        else: 
            sexlines2 = np.append(sexlines2,'PARAMETERS_NAME  '+paramfile)
    # FLAG/MASK file 
    if flagfile is not None:
        sexlines2 = np.append(sexlines2,'FLAG_IMAGE  '+flagfile)
        sexlines2 = np.append(sexlines2,'FLAG_TYPE   OR')
    # WEIGHT file 
    if wtfile is not None:
        sexlines2 = np.append(sexlines2,'WEIGHT_IMAGE  '+wtfile)
        sexlines2 = np.append(sexlines2,'WEIGHT_TYPE   MAP_WEIGHT')
    # Write the file 
    dln.writelines(configfile,sexlines2)


def sex2daophot(catfile,fitsfile,daofile):
    """
    This program converts a SExtractor output file 
    to DAOPHOT format.  Currently only coo. 
 
    Parameters
    ----------
    catfile : str
       The SExtractor catalog filename. 
    fitsfile : str
       The name of the associated FITS file. 
    daofile : str,
       The name of the output DAOPHOT file. 

    Returns
    -------
    The catalog is written to DAOFILE. 
 
    Example
    -------
    sex2daophot('F1-12340056_01.cat','F1-12340056_01.fits','F1-12340056_01.coo')
 
    By D. Nidever  Jan 2019 
    Translated to python by D. Nidever  Sep 2022
    """

    # Check that the needed files exist 
    if os.path.exists(catfile) == False:
        raise ValueError(catfile+' NOT FOUND')
    if os.path.exists(fitsfile) == False:
        raise ValueError(fitsfile+' NOT FOUND')

    #------------------------------------- 
    # Load sextractor output file 
    # default.param specifies the output columns 
    if utils.file_isfits(catfile) == False: 
        #  fields = ['ID','X','Y','MAG','ERR','FLAGS','STAR'] 
        fields = ['NUMBER','X_IMAGE','Y_IMAGE','MAG_APER','MAGERR_APER','FLAGS','CLASS_STAR'] 
        sex = Table.read(catfile,format='ascii',names=fields) 
    else: 
        sex = Table.read(catfile,1) 
        if len(sex.colnames)==0:
            sex = Table.read(catfile,2) 
    nsex = len(sex) 

    #------------------------------------- 
    # Get meta-data from the FITS file 
    head = io.readfile(fitsfile,header=True) 
    naxis1 = head['NAXIS1'] 
    naxis2 = head['NAXIS2']
    saturate = head['SATURATE']
    rdnoise,rdnoisekey = io.getrdnoise(fitsfile) 
    gain,gainkey = io.getgain(fitsfile) 
    lowbad = 1.0 
    thresh = 20.0 
     
    # Header values:  this information comes from daophot2.pdf pg.69 
    # NL: Originally meant "number of lines" but not anymore 
    # NX: size of X-dimension of image in pixels 
    # NY: size of Y-dimension of image in pixels 
    # LOWBAD: lower good data limit, calculated by FIND 
    # HIGHBAD: upper good data limit, specified in option file 
    # THRESH: threshold calculated by FIND 
    # AP1: radius (pixels) of the first aperture used by PHOTOMETRY 
    # PH/ADU: gain in photons/ADU used when running FIND 
    # RDNOISE: rdnoise (ADU) used when running FIND 
    # FRAD: value of fitting radius 
     
    # Making ALS structure for new SEX sources 
    dt = [('ID',int),('X',float),('Y',float),('MAG',float),('ERR',float),
          ('SKY',float),('ITER',float),('CHI',float),('SHARP',float)]
    dao = np.zeros(nsex,dtype=np.dtype(dt))
    dao['ID'] = sex['NUMBER']
    dao['X'] = sex['X_IMAGE']
    dao['Y'] = sex['Y_IMAGE']
    dao['MAG'] = sex['MAG_APER'] 
    dao['ERR'] = sex['MAGERR_APER'] 
    if 'BACKGROUND' in sex.columns:
        dao['SKY'] = sex['BACKGROUND']
    else: 
        dao['SKY'] = 0.0 
    dao['ITER'] = 1 
    dao['CHI'] = 1.0 
    dao['SHARP'] = 0.0 
    ndao = nsex 
     
    #NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD 
    #  1  2046  4094  1472.8 38652.0   80.94    0.00    3.91    1.55    3.90 
    # 
    #      1  1434.67    15.59   -0.045    0.313    0.873    1.218 
    #      2   233.85    18.42   -0.018    0.218   -0.781    1.433 
    #    ID      X         Y       MAG     SHARP    ROUND    ROUND2 
    with open(daofile,'w') as f:
        f.write(' NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD\n')
        f.write('%3d%6d%6d%8.1f%8.1f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % (1,naxis1,naxis2,lowbad,saturate,thresh,3.0,gain,rdnoise/gain,3.9))
        f.write(' \n')
        # Write the data 
        for i in range(ndao): 
            f.write('%7d%9.2f%9.2f%9.3f%9.3f%9.3f%9.3f\n' % (dao['ID'][i],dao['X'][i],dao['Y'][i],dao['MAG'][i],0.6,0.0,0.0))


def getpsf(base,fake=False,verbose=False,logger=None):
    """
    This runs the DAOPHOT routines to create the PSF for an image. 
 
    Parameters
    ----------
    base : str
       The base name of the FITS and MCH files. 
    fake : boolean, optional
       Run for artificial star tests.  Default is False.
    verbose : boolean, optional
       Verbose output to the screen.  Default is False.
    logger : logging object, optional
       A logging object to use for output to the screen.
 
    Returns
    -------
    psffile : str
      The name of the DAOPHOT PSF file. 
 
    Example
    -------
    psffile = getpsf('F1-23430911_10')
 
    By D.Nidever  Feb 2019 
    Translated to python by D.Nidever  Sep 2022
    """

    if logger is None:
        logger = dln.basiclogger()
     
    # Make SExtractor output parameter file 
    secols = ['NUMBER','X_IMAGE','Y_IMAGE','MAG_APER(1)','MAGERR_APER(1)','MAG_AUTO','MAGERR_AUTO','BACKGROUND',
              'THRESHOLD','ISOAREA_IMAGE','A_IMAGE','B_IMAGE','THETA_IMAGE','ELLIPTICITY','FWHM_IMAGE','FLAGS',
              'IMAFLAGS_ISO(1)','NIMAFLAGS_ISO(1)','CLASS_STAR'] 
    dln.writelines(base+'.param',secols)
     
    # Create the SExtractor config file 
    mksexconfig(base+'.fits',configfile=base+'.sex',catfile=base+'.cat',flagfile=base+'.bpm.fits',                    
                paramfile=base+'.param',cattype='FITS_1.0')
     
    # Run SExtractor for detection 
    logger.info('Running Source Extractor')
    out = subprocess.check_output(['sex',base+'.fits','-c',base+'.sex'],shell=False,stderr=subprocess.STDOUT)
    for o in out.decode().split('\n'): logger.info(o)
    hd = io.readfile(base+'.cat',exten=1,header=True) 
    nsources = hd.get('NAXIS2')
    if nsources < 1: 
        raise ValueError('Only '+str(nsources)+' sources. Need at least 1 to create a PSF')
     
    # Apply cuts to get good stars 
    sex = Table.read(base+'.cat',1)
    if len(sex.colnames) == 1:
        sex = Table.read(base+'.cat',2) 
    nsex = len(sex) 
    si = np.argsort(sex['FWHM_IMAGE']) 
    fwhm80 = np.percentile(sex['FWHM_IMAGE'],80)
    bad = ((sex['IMAFLAGS_ISO'] > 0) | (sex['CLASS_STAR'] < 0.5) | (sex['ELLIPTICITY'] > 0.8) | 
           (sex['FWHM_IMAGE'] > fwhm80) | ( ((sex['FLAGS'] & 8) > 0) | ((sex['FLAGS'] & 16) > 0) | 
                                            (1.087/sex['MAGERR_APER'] <= 10)))
    nbad = np.sum(bad)
    ngood = np.sum(~bad)
    # Not enough stars, remove class_star cut and raise S/N cut 
    if ngood < 50: 
        bad = ((sex['IMAFLAGS_ISO'] > 0) | (sex['ELLIPTICITY'] > 0.8) | (sex['FWHM_IMAGE'] > fwhm80) |
               ((sex['FLAGS'] & 8) > 0) | ((sex['FLAGS'] & 16) > 0) | (1.087/sex['MAGERR_APER'] <= 7))
    nbad = np.sum(bad)
    ngood = np.sum(~bad)
    # No sources left 
    if ngood == 0: 
        logger.info('No sources left after stellar cuts')
        return 
    sex = sex[~bad] 
    nsources = ngood
    sex.write(base+'.cat',format='fits',overwrite=True)

    # Convert to DAOPHOT coo format 
    sex2daophot(base+'.cat',base+'.fits',base+'.coo')
     
    # Check that there are enough stars for our VA setting 
    #   VA=2 quadratic spatial variations, 6 stars minimum 
    #   VA=1 linear spatial variations, 3 stars minimum 
    #   VA=0 empirical corrections but no spatial variations, 1 star 
    #   minimum 
    #   VA=-1 analytical, 1 star minimum 
    opt = io.readopt(base+'.opt')
    optlines = dln.readlines(base+'.opt')
    optlines = np.char.array(optlines)
    vaval = opt['VA']
    minsources = [1,1,3,6] 
    if nsources < minsources[int(vaval)+1]: 
        # Get new VA value, largest allowed for this number of sources 
        if nsources >= 3: 
            newvaval = 1.0 
        else: 
            newvaval = 0.0 
        logger.info('Only '+str(nsources)+' sources and need '+str(minsources[int(vaval)+1])+' for VA='+str(vaval)+'. Lower to VA='+str(newvaval)) 
        newoptlines = optlines.copy()
        newoptlines[vaind] = ('VA = %8.2f' % newvaval)
        dln.writelines(base+'.opt',newoptlines)
     
    # Sometimes the filenames get too long for DAOPHOT 
    # use temporary files and symlinks
    tid,tbasefile = tempfile.mkstemp(prefix='cmb',dir='.')  # create base, leave so other processes won't take it   
    tbase = os.path.basename(tbasefile)
    tfits = tbase+'.fits'
    if os.path.exists(tfits): os.remove(tfits)
    os.symlink(base+'.fits',tfits)
    topt = tbase+'.opt'
    if os.path.exists(topt): os.remove(topt)
    os.symlink(base+'.opt',topt)
    taopt = tbase+'.als.opt'
    if os.path.exists(taopt): os.remove(taopt)
    os.symlink(base+'.als.opt',taopt)
    tcoo = tbase+'.coo'
    if os.path.exists(tcoo): os.remove(tcoo)
    os.symlink(base+'.coo',tcoo)
     
    # Get the PSF of the combined image 
    if os.path.exists('getpsfnofind.sh')==False:
        raise ValueError('getpsfnofind.sh script not found')
    os.chmod('getpsfnofind.sh',0o755)
    os.chmod('lstfilter.py',0o755)
    logger.info('Running DAOPHOT')
    out = subprocess.check_output(['./getpsfnofind.sh',tbase],shell=False,stderr=subprocess.STDOUT)
    for o in out.decode().split('\n'): logger.info(o)     

    # If getpsf failed or has NaN, change to VA=0 
    psfline1 = dln.readlines(tbase+'.psf',nreadline=1)[0]
    if os.path.exists(tbase+'.psf')==False or os.path.getsize(tbase+'.psf')==0 or psfline1.find('NaN') != -1:
        optlines = dln.readlines(base+'.opt')
        optlines = np.char.array(optlines)
        opt = io.readopt(base+'.opt')
        vaval = opt['VA']
        vaind, = np.where(optlines.find('VA') > -1)
        newoptlines = optlines.copy()
        newoptlines[vaind] = 'VA = %8.2f' % 0.0
        dln.writelines(base+'.opt',newoptlines)
        logger.info('getpsfnofind.sh failed or NaN.  Changing to VA=0.  Trying again.')
        out = subprocess.run(['getpsfnofind.sh',tbase],shell=False,stderr=subprocess.STDOUT)
        for o in out.decode().split('\n'): logger.info(o)     

    # Delete the temporary symlinks 
    for f in [tbase,tfits,topt,taopt,tcoo]:
        if os.path.exists(f): os.remove(f)

    # Getpsf succeeded, rename files 
    psfline1 = dln.readlines(tbase+'.psf',nreadline=1)[0]
    if os.path.exists(tbase+'.psf') and os.path.getsize(tbase+'.psf')>0 and psfline1.find('NaN') == -1:
        outfiles = glob(tbase+'*')
        noutfiles = len(outfiles)
        renamefiles = [o.replace(tbase,base) for o in outfiles]
        for i in range(len(outfiles)):
            if os.path.exists(renamefiles[i]): os.remove(renamefiles[i])
            shutil.move(outfiles[i],renamefiles[i])
     
    # No PSF file found 
    if os.path.exists(base+'.psf')==False or os.path.getsize(base+'.psf')==0 or psfline1.find('NaN') != -1:
        logger.info('Could not create PSF for '+base)
        return None
 
    return base+'.psf'
