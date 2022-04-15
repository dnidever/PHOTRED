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
from astropy.time import Time
from dlnpyutils import utils as dln
import struct
from itertools import zip_longest
from itertools import accumulate
from io import StringIO
from . import utils

def make_parser(fieldwidths,fieldtypes=None):
    # https://stackoverflow.com/questions/4914008/how-to-efficiently-parse-fixed-width-files
    cuts = tuple(cut for cut in accumulate(abs(fw) for fw in fieldwidths))
    pads = tuple(fw < 0 for fw in fieldwidths) # bool flags for padding fields
    if fieldtypes is None:
        flds = tuple(zip_longest(pads, (0,)+cuts, cuts))[:-1]  # ignore final one        
        slcs = ', '.join('line[{}:{}]'.format(i, j) for pad, i, j in flds if not pad)
    else:
        tdict = {'s':'str','d':'int','f':'float'}
        ftypes = [tdict[ft] for ft in fieldtypes]
        flds = tuple(zip_longest(ftypes,pads, (0,)+cuts, cuts))[:-1]  # ignore final one        
        slcs = ', '.join('{}(line[{}:{}])'.format(ftype, i, j) for ftype, pad, i, j in flds if not pad)        
    parse = eval('lambda line: ({})\n'.format(slcs))  # Create and compile source code.
    # Optional informational function attributes.
    parse.size = sum(abs(fw) for fw in fieldwidths)
    if fieldtypes is None:
        parse.fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                                   for fw in fieldwidths)
    else:
        parse.fmtstring = ' '.join('{}{}'.format(a[0],a[1]) for a in zip(fieldwidths,fieldtypes))
    return parse

def loadsetup(fake=False,setupdir=None,std=False):
    """

    Parameters
    ----------
    fake      Load the fakered.setup file. 
    setupdir  The directory in which to look for the setup file. 
 
    Returns
    -------
    setup      The setup file.  It is a string array with 
                 dimensions of 2xN_parameters.  READPAR can be 
                 used to read the parameters. 

    Example
    -------

    setup = loadsetup()
 
    By D.Nidever  March 2008 
    Translated to Python by D. Nidever,  April 2022
    """
 
    if setupdir is None:  # default setup directory 
        setupdir = '.' 
     
    # Type of setup file 
    setupfile = 'photred'
    if std:
        setupfile = 'stdred'
    if fake:
        setupfile = 'fakered' 
     
    # LOAD THE SETUP FILE 
    #-------------------- 
    # This is a 2xN array.  First colume are the keywords 
    # and the second column are the values. 
    # Use READPAR.PRO to read it 
    setupfiles = glob(setupdir+'/'+setupfile+'.*setup')
    nsetupfiles = len(setupfiles)
    if (nsetupfiles < 1): 
        raise ValueError('NO '+strupcase(setupfile)+' SETUP FILE')
    if (nsetupfiles > 1):
        raise ValueError('MORE THAN ONE '+strupcase(setupfile)+' SETUP FILE FOUND')
     
    # Read the setup file
    lines = dln.readlines(setupfiles[0],comment='#')
    nlines = len(files)
     
    # Process the lines
    setup = {}
    if nlines > 0:
        for l in lines:
            if l.strip()!='':
                arr = l.split()
                if len(arr)==1:
                    setup[arr[0]] = ''
                else:
                    setup[arr[0]] = arr[1]
         
    # No lines to process 
    else: 
        raise ValueError(setupfiles[0],' HAS NOT LINES TO PROCESS')

    return setup

def loadmch(mchfile):
    """
    This loads a DAOMATCH/DAOMASTER mch file 
     
    Parameters
    ----------
    mchfile  The MCH filename 
     
    Returns
    -------
    files    The list of files in the MCH file 
    trans    The transformation equations in a Nfilesx6 array. 
    magoff   The magnitude offsets and errors array. 
     
    Example
    -------

    files,trans,magoff = loadmch('ccd1001.mch')
     
    By D.Nidever   February 2008 
    Translated to Python by D. Nidever,  April 2022
    """ 
     
    # Test the file 
    if os.path.exists(mchfile)==False:
        raise ValueError(mchfile+' NOT FOUND')

    # Read in the file
    lines = dln.readlines(mchfile)

    # Creating the trans array
    lines2 = [re.sub("'","",l) for l in lines]
    nlines = len(lines)

    arr = lines2[0].split()
    ntrans = len(arr)-3  # first line is the samee, last two are mag offsets 
 
    # Getting the file names
    arr2 = [l.split() for l in lines2]
    files = [a[0] for a in arr2]
 
    # Initializing the array
    trans = np.zeros((nlines,ntrans),float)
    magoff = np.zeros((nlines,2),float)
    # Filling the aray 
    for i in range(nlines):
        arr = lines2[i].split()
        trans[i,:] = arr[1:ntrans+1]
        magoff[i,:] = arr[ntrans+1:]
 
    return files,trans,magoff

def loadmakemag(filename):
    """"
    This loads the .makemag file output in ALLFRAME.PRO. 
    This used to be handled by LOADRAW.PRO, but they diverged. 
    This is basically the old version of loadraw.pro. 
     
    Parameters
    ----------
    filename   The name of the .makemag file 
     
    Returns
    -------
    phot       A structure with the makemag data 
    head       A string array with the makemag header 
     
    Example
    ------
    
    phot,head = loadmakemag('temp.makemag')
     
    By D. Nidever   Feb. 2008 (basically a copy of loadals.pro 
    Translated to Python by D. Nidever,  April 2022
    """
     
    if os.path.exists(filename)==False:
        raise ValueError(filename+' DOES NOT EXIST')
     
    # Is this an ALS or PHOT file
    with open(filename,'r') as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
        
    # This is an ALS file
    arr1 = line1.split()
    if arr1[0]=='NL' and line3.strip()=='': 
        # Figure out the number of columns
        f = open(filename,'r')
        line1 = f.readline().replace("\n", "")
        line2 = f.readline().replace("\n", "")
        line3 = f.readline().replace("\n", "")
         
        # First line for the first star
        line4 = f.readline().replace("\n", "")
        instr = line4 

        # Check for continuation lines 
        endflag = 0 
        nstarline = 1 
        continuation = False
        line4 = '1'
        while (endflag != 1) and line4 != '':
            line4 = f.readline().replace("\n", "")
         
            # If there are too many frames/columns then these 
            # go on separate lines and lead with 27 blank spaces 
         
            # This is a continuation line 
            if line4[0:16].strip()=='':
                trial = line4[33:34]
                if trial == ' ': 
                    nspaces = 24 
                else: 
                    nspaces = 25
                instr1 = line4[nspaces:]
                instr += instr1 
                nstarline += 1 
                continuation = True 
            else:
                endflag = 1 
        f.close()

        
        # Now parse the long line for a single star with a formatted read 
        #fmt = '(I7,2A9,'+strtrim(2*nmag+2,2)+'F9.4)' 
        nmag = int(  (len(instr)-(7+2*9+2*9)) / 9 / 2 )
        ncol = nmag*2+5 
         
        # ID  X  Y  MAG1  ERR1  MAG2  ERR2 ...  CHI SHARP 
        #nmag = (ncol-5)/2 
         
        # Stars in this file
        numstar = int((dln.numlines(filename)-3 )/nstarline)
         
        # LOADING THE DATA 
        #------------------              
        # mastable is where everything is stored, id, x, y, unsolved magnitudes, chi, sharp
        dt = [('id',int),('x',float),('y',float)]
        for m in np.arange(nmag): dt += [('mag'+str(m+1),float),('err'+str(m+1),float)]
        dt += [('chi',float),('sharp',float)]
        phot = np.zeros(numstar,dtype=np.dtype(dt))
            
        # Reading in the magnitude file
        f = open(filename,'r')
        line1 = f.readline().replace("\n", "")
        line2 = f.readline().replace("\n", "")
        line3 = f.readline().replace("\n", "")
        head = [line1,line2] 

        # WRITE (1,111) IDMAST(IMASTR), POSIT1, POSIT2, 
        # .            ((DATUM(J,I), J=1,2), I=1,NFRM), SUMCHI, SUMSHP 
        #        111       FORMAT (I7, 2A9, 12F9.4:/ (25X, 12F9.4)) 
        fieldwidths = tuple([7,9,9]+(2*nmag+2)*[9])
        fieldtypes = tuple(['d','f','f']+(2*nmag+2)*['f'])
        parser = make_parser(fieldwidths,fieldtypes)
        
        # Loop through the stars 
        for j in np.arange(numstar): 
            instr = '' 
            instr1 = ' '
            inline = np.zeros(2*nmag+5,float)
                 
            # Loop through the lines per star 
            for k in np.arange(nstarline):
                instr1 = f.readline().replace("\n", "")
                if k > 0: 
                    # There are leading spaces, 24 or 25 
                    # Use the first character AFTER the first column to figure out 
                    #   how many spaces we need to strip off 
                    # The MAG/ERR columns are F9.4, so if there are 25 spaces then 
                    #   the 34th character (index=33) should be the last digit of MAG1 (25+9=34) 
                    #   and be a number.  If there are 24 leading space then the 
                    #   34th character will right after MAG1 and be a space.
                    trial = instr1[33:34]
                    if trial == ' ': 
                        nspaces = 24 
                    else: 
                        nspaces = 25 
                    instr1 = instr1[nspaces:]
                #if k gt 0 then instr1=strmid(instr1,25) ; 2nd and later lines have 25 leading spaces 
                instr += instr1 
                
            out = parser(instr)
            phot[j] = out
        f.close()

        # Convert to astropy table
        phot = Table(phot)
        
    # Not a raw/makemag file 
    else: 
        print('This is NOT a MAKEMAG file')
        return None,None

    return phot,head

def readopt(optfile):
    """
    Read in a DAOPHOT option file.

    Parameters
    ----------
    optfile : str
       Filename of optional file.

    Returns
    -------
    opt : dict
       Dictionary of DAOPHOT option parameters.

    Example
    -------

    opt = readopt('filename.opt')

    Written by D. Nidever,  April 2022
    """

    if os.path.exists(optfile)==False:
        raise ValueError(optfile+' NOT FOUND')

    # Read in the file
    lines = dln.readlines(optfile)
    nlines = len(lines)
    opt = {}
    # Parse the values
    for i in range(nlines):
        line = lines[i]
        if len(line)>0 and line[0] != '#':  # skip comment lines
            arr = line.split()
            if len(arr)==3 and arr[1]=='=':
                key = arr[0]
                val = arr[2]
                if val.isnumeric():
                    val = int(val)
                elif val.find('.') or val.find('e') or val.find('E'):
                    if utils.isfloat(val):
                        val = float(val)
                opt[key] = val
    return opt

def readpar(array,keyword):
    """
    This allows you to get a parameter value from 
    a 2xN array of keyword/value pairs. 
    This is similar to getting keyword values from 
    headers with SXPAR.PRO. 
    
    Parameters
    ----------
    array    A 2xN array of keyword-value pairs 
    keyword  A keyword string for which to return the value. 
            Case insensitive. 
 
    Returns
    -------
    value    The value corresponding to the input keyword 
              is output.  If none is found then '0' is returned. 
              If the keyword exists in array but has not value 
              then an empty string '' is returned. 
    count   The number of parameters found by READPAR.  count=0 
              if the parameter was not found. 

    Example
    -------

    value = readpar(setup,'MOSAIC') 
 
    By D. Nidever    Oct. 2007 
    Translated to Python by D. Nidever,  April 2022
    """
    
    #sz = size(array)
    #if sz[0] != 2 or sz[1] != 2: # must be 2xN 
    #    return None
    #if size(array,/type) != 7: # must be string 
    #    return None 
     
    # Looking for keyword 
    keyword2 = strlowcase(str(keyword[0],2)) 
    keys = reform(array[0,:]) 
    values = reform(array[1,:])
    gd, = dln.grep(keys==keyword2)
    if len(gd)==0:
        return None
    value = str(values[gd[0]])  # returning the first value 
    return value

def mkopt(base=None,meta=None,VA=1,LO=7.0,TH=3.5,LS=0.2,HS=1.0,LR=-1.0,HR=1.0,
          WA=-2,AN=-6,EX=5,PE=0.75,PR=5.0,CR=2.5,CE=6.0,MA=50.0,RED=1.0,WA2=0.0,
          fitradius_fwhm=1.0,HI=None,RD=None,GA=None,FW=None,logger=None):
    """
    Create the DAOPHOT and ALLSTAR option files (.opt and .als.opt) for an exposure.

    Parameters
    ----------
    base : str
         The base name to use for the option files.  The DAOPHOT option file will
         be called `base`.opt and the ALLSTAR option file `base`.als.opt
    meta : astropy dictionary
         The metal-data dictionary for the image.    
    VA : int, default = 1
       The variable type of PSF to use.
       -1: Analytic PSF only
        0: Analytic PSF and look-up table of empirical corrections
        1: linear variations across the field
        2: quadratic variations across the field
    LO : float, default = 7.0
       Low good datum (7. works fine on most imags).
    TH : float, default = 3.5
       Threshold in sigma above the background (3.5 works fine).
    LS : float, default = 0.2
       Lower sharpness cutoff.
    HS : float, default = 1.0
       High sharpness cutoff.
    LR : float, default = -1.0
       Lower roundness cutoff.
    HR : float, default = 1.0
       High roundness cutoff.
    WA : int, default = -2
       Watch progress for DAOPHOT.  Determines what output is displayed.
    AN : int, default = -6
       Analytic model PSF.
        1: Gaussian (3 pararameters)
        2: Moffat function (3 parameters), beta=1.5
        3: Moffat function (3 parameters), beta=2.5
        4: Lorentz function (3 parameters)
        5: Penny function, Gauss+Lorentz (4 parameters), G+L are parallel
        6: Penny function (5 parameters), G and L can be in different directions
        A negative sign in front means to try all functions up to X and pick the best one.
    EX : int, default = 5
       Extra PSF cleaning passes.
    PE : float, default = 0.75
       Percent error due to the uncertainty in the fine-scale structure of the flat field.
    PR : float, default = 5.0
       Profile error due to the incompleteness of the PSF model.
    CR : float, default = 2.5
       Clipping range.  Used to remove outlier pixels. Parameter "a" in the formula given in
       Stetson 1987, PASP, 99, 191, section III.D.2.d "Resisting bad data".
    CE : float, default = 6.0
       Clipping exponent.  Parameter b in above clipping formula.
    MA : float, default = 50.0
       Maximum group size
    RED : float, default = 1.0
        Redetermine centroid (0 = no, 1 = yes).
    WA2 : float, default = 0.0
        Watch progress for ALLSTAR.      
    fitradius_fwhm : float, default = 1.0
        The fitting radius size in units of the seeing FWHM for the area to be fit.
    HI : float, optional
       High good datum.  Normally set by `saturate` from `meta`.
    RD : float, optional
       The read noise in electrons. Normally set by `rdnoise` from `meta`.
    GA : float, optional
       The gain in electrons/ADU. Normally set by `gain` from `meta`.
    FW : float, optional
       The seeing FWHM in pixels.  Normally set by `fwhm`/`pixscale` from `meta`.
    logger : logger object, optional
           The Logger to use for logging output.

    Returns
    -------
    Nothing is returned.  The DAOPHOT option file is written to `base`.opt and the ALLSTAR
    option file to `base`.als.opt.

    Example
    -------

    .. code-block:: python

        mkopt("image",meta)

    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % MAKING THE OPT FILES
    #
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
    #
    # (2) ALLSTAR parameters
    # 
    # CR    : Clipping range (leave it)
    # CE    : Clipping exponent (leave it)
    # MA    : Maximum group size
    # RED   : Redetermine centroid (0 = no, 1 = yes)
    #
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
    # VA  defined above
    #AN = -6     # It will try all PSF models (#1-6) and use the one with the lowest chi value
    #EX =  5     # extra PSF passes

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary

    optfile = base+".opt"
    alsoptfile = base+".als.opt"

    # Get frame specific parameters from meta if necessary
    if GA is None: GA = meta['gain']
    if RD is None: RD = meta['rdnoise']
    if FW is None: FW = meta['fwhm'] / meta['pixscale']
    if HI is None: HI = meta['saturate']


    # Calculating some things
    FW = np.min([ FW , 20 ])            # daophot won't accept anything higher than 20
    RE = RD/GA
    FI = np.min([ fitradius_fwhm*FW , 51 ])                  # daophot won't accept anything higher than 51
    PS = np.min([ (4.0*FW) , 51 ])       # daophot won't accept anything higher than 51
    IS = np.min([ (FI - 1.0) , 35 ])     # daophot won't accept anything higher than 35
    OS = np.min([ (PS + 1.0) , 100 ])    # daophot won't accept anything higher than 100

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
    nanot = len(anotarr)

    # Delete file if it exists
    if os.path.exists(optfile):
        os.remove(optfile)
    # Write opt file
    f = open(optfile,'w')
    for j in range(len(outarr)):
        if anotarr[j] == "HI":
            f.write("%2s = %8d\n" % (anotarr[j], outarr[j]))
        else:
            f.write("%2s = %8.2f\n" % (anotarr[j], outarr[j]))
    f.close()

    # Writing the ALLSTAR parameter file
    #-----------------------------------
    #
    # FI    : fitting radius
    # IS    :  ??
    # OS    :  ??
    # RED   : Redetermine centroid (0 = no, 1 = yes)
    # WA2   : Watch progress
    # PE    : Percent error
    # PR    : Profile error
    # CR    : Clipping range (leave it)
    # CE    : Clipping exponent (leave it)
    # MA    : Maximum group size

    outarr2 = [FI,IS,OS,RED,WA2,PE,PR,CR,CE,MA]
    anotarr2 = ['FI','IS','OS','RE','WA','PE','PR','CR','CE','MA']
    nanot2 = len(anotarr2)
    form = '(A5,F8.2)'

    # Delete file if it exists
    if os.path.exists(alsoptfile):
        os.remove(alsoptfile)
    # Write opt file
    f = open(alsoptfile,'w')
    for j in range(len(outarr2)):
        f.write("%2s = %8.2f\n" % (anotarr2[j], outarr2[j]))
    f.close()

    logger.info("Created "+optfile+" and "+alsoptfile)

def loadals(filename,silent=False):
    """
    This loads the ALLSTAR photometry file (.als).
 
    Parameters
    ----------
    filename : str
       The name of the ALS file.
    silent : boolean, optional
       Don't print anything to the screen.  Default is False.
 
    Returns
    -------
    phot : astropy table
       A table with the ALS data.
    head : list
       A string list with the ALS header.
 
    Example
    -------
    
    phot,head = loadals('obj2153_1.als')
 
    By D. Nidever   January 2007 
    Translated to Python by D. Nidever,  April 2022
    """
     
    if os.path.exists(filename)==False:
        raise ValueError(filename+' NOT FOUND')
     
    nlines = dln.numlines(filename) 
    if nlines < 4:
        if silent==False:
            print('No sources in '+filename)
        return None,None
     
    # Is this an ALS or PHOT file
    with open(filename,'r') as f:
        line1 = f.readline().replace("\n","")
        line2 = f.readline().replace("\n","")
        line3 = f.readline().replace("\n","")
        line4 = f.readline().replace("\n","")        
     
    # This is an ALS file
    arr1 = line1.split()
    if arr1[0] == 'NL' and line3.strip() == '':
        numstars = dln.numlines(filename)-3        
        dt = [('id',int),('x',float),('y',float),('mag',float),('err',float),('sky',float),
              ('niter',int),('chi',float),('sharp',float)]
        phot = np.zeros(numstars,dtype=np.dtype(dt))

        fieldwidths = tuple([7,9,9,9,9,9,9,9,9])
        fieldtypes = tuple(['d','f','f','f','f','f','f','f','f'])
        parser = make_parser(fieldwidths,fieldtypes)

        f = open(filename,'r')
        line1 = f.readline().replace("\n","")
        line2 = f.readline().replace("\n","")
        line3 = f.readline().replace("\n","")
        head = [line1,line2] 
        
        # Loop through the stars 
        for i in np.arange(numstars):
            line = f.readline().replace("\n","")
            out = parser(line)
            phot[i] = out
        f.close()
        phot = Table(phot)  # convert to astropy table

    # This is a PHOT file 
    else: 
        raise ValueError('This is NOT an ALLSTAR output file')

    return phot,head
     

def loadcoo(filename,silent=False):
    """
    This loads the DAOPHOT coordinates file (.coo).

    Parameters
    ----------
    filename : str
       The name of the coo file 
    silent : boolean, optional
       Don't print anything to the screen.  Default is False.
 
    Returns
    -------
    phot : astropy table
       A table with the coo data.
    head : list
       A string list with the coo header.
 
    Example
    -------
    
    phot,head = loadcoo('obj2153_1.coo')
 
    By D. Nidever   January 2007 
    Translated to Python by D. Nidever,  April 2022
    """
     
    if os.path.exists(filename)==False:
        raise ValueError(filename+' NOT FOUND')
     
    nlines = dln.numlines(filename) 
    if nlines < 4:
        if silent==False:
            print('No sources in '+filename)
        return None,None
     
    # Is this an DAOPHOT file
    with open(filename,'r') as f:
        line1 = f.readline().replace("\n","")
        line2 = f.readline().replace("\n","")
        line3 = f.readline().replace("\n","")
        line4 = f.readline().replace("\n","")        
     
    # This is a DAOPHOT file
    arr1 = line1.split()
    if arr1[0] == 'NL' and line3.strip() == '':
        numstars = dln.numlines(filename)-3        
        dt = [('id',int),('x',float),('y',float),('mag',float),('sharp',float),('round',float),('round2',float)]
        tab = np.zeros(numstars,dtype=np.dtype(dt))

        fieldwidths = tuple([7,9,9,9,9,9,9])
        fieldtypes = tuple(['d','f','f','f','f','f','f'])
        parser = make_parser(fieldwidths,fieldtypes)

        f = open(filename,'r')
        line1 = f.readline().replace("\n","")
        line2 = f.readline().replace("\n","")
        line3 = f.readline().replace("\n","")
        head = [line1,line2] 
        
        # Loop through the stars 
        for i in np.arange(numstars):
            line = f.readline().replace("\n","")
            out = parser(line)
            tab[i] = out
        f.close()
        tab = Table(tab)  # convert to astropy table

    # This is not a DAOPHOT file 
    else: 
        raise ValueError('This is NOT an DAOPHOT output file')

    return tab,head
     

def loadaper(filename,silent=False):
    """
    This loads the DAOPHOT aperture photometry file (.ap).
 
    Parameters
    ----------
    filename : str
       The name of the aperture photometry file.
    silent : boolean, optional
       Don't print anything to the screen.  Default is False.
 
    Returns
    -------
    phot : astropy table
       A table with the aper data.
    head : list
       A string list with the aper header.
 
    Example
    -------
    
    phot,head = loadaper('obj2153_1.ap')
 
    By D. Nidever   January 2007 
    Translated to Python by D. Nidever,  April 2022
    """
     
    if os.path.exists(filename)==False:
        raise ValueError(filename+' NOT FOUND')
     
    nlines = dln.numlines(filename) 
    if nlines < 4:
        if silent==False:
            print('No sources in '+filename)
        return None,None
     
    # Is this an aperture photometry file
    with open(filename,'r') as f:
        line1 = f.readline().replace("\n","")
        line2 = f.readline().replace("\n","")
        line3 = f.readline().replace("\n","")
        line4 = f.readline().replace("\n","")
        line5 = f.readline().replace("\n","")        
     
    # This is an ALS file
    arr1 = line1.split()
    arr2 = line2.split()    
    if arr1[0] == 'NL' and arr2[0].strip()[0]=='2' and line3.strip()=='' and line4.strip()=='':
        
        # This is what a Daophot aperture photometry file looks like 
        # NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD 
        #  2  2040  2047  1837.1 26000.0  114.56    4.00    2.20    2.86    1.98 
        # 
        # 
        #      1    3.000    1.000   97.999 
        #      2055.039 41.29  0.00  9.9999 
        # 
        #      2  397.000    1.000   97.999 
        #      2066.848 37.25  0.12  9.9999 
        # 
        #  The columns are: ID, X, Y, Mag1, Mag2, etc.. 
        #                   Sky, St.Dev. of sky, skew of sky, Mag1err, Mag2err, etc. 
        # 
        
        # Figure out the number of columns
        # First line for the first star
        ncol = (len(line5)-7)//9 + 1
        nmag = ncol-3
         
        # ID  X  Y  MAG1  ERR1  MAG2  ERR2 ...  CHI SHARP 
        #nmag = (ncol-5)/2 
         
        # Stars in this file
        #  3 header lines
        #  3 lines per star except the last one
        numstar = int((dln.numlines(filename)-3 )/3)
         
        # LOADING THE DATA 
        #------------------              
        dt = [('id',int),('x',float),('y',float)]
        for m in np.arange(nmag): dt += [('mag'+str(m+1),float)]
        dt += [('sky',float),('skysig',float),('skyskew',float)]
        for m in np.arange(nmag): dt += [('err'+str(m+1),float)]        
        phot = np.zeros(numstar,dtype=np.dtype(dt))
            
        # Reading in the magnitude file
        f = open(filename,'r')
        line1 = f.readline().replace("\n", "")
        line2 = f.readline().replace("\n", "")
        line3 = f.readline().replace("\n", "")
        head = [line1,line2] 

        # WRITE (1,111) IDMAST(IMASTR), POSIT1, POSIT2, 
        # .            ((DATUM(J,I), J=1,2), I=1,NFRM), SUMCHI, SUMSHP 
        #        111       FORMAT (I7, 2A9, 12F9.4:/ (25X, 12F9.4))
        #      1   39.000   19.740   99.999   99.999   99.999   99.999   99.999   99.999   99.999   99.999   99.999   99.999   99.999   99.999
        #      1080.900 17.42  0.00  9.9999   9.9999   9.9999   9.9999   9.9999   9.9999   9.9999   9.9999   9.9999   9.9999   9.9999   9.9999
        fieldwidths = tuple([7]+(nmag+2)*[9]+[14,6,6]+nmag*[9])
        fieldtypes = tuple(['d']+(2*nmag+5)*['f'])
        parser = make_parser(fieldwidths,fieldtypes)
        
        # Loop through the stars 
        for i in np.arange(numstar):
            line1 = f.readline().replace("\n", "")   # blank
            line2 = f.readline().replace("\n", "")   # line 1: ID, X, Y, MAGs
            line3 = f.readline().replace("\n", "")   # line 2: sky, skysig, skyskew, ERRs
            line = line2 + line3
            out = parser(line)
            phot[i] = out
        f.close()

        # Convert to astropy table
        phot = Table(phot)
        
    # Not a DAOPHOT aperture photometry file
    else:
        raise ValueError('This is NOT a DAOPHOT aperture photometry file')

    return phot,head


def writeals(outfile,phot,head):
    """
    This writes a photometry structure in the ALLSTAR output
    format.

    Parameters
    ----------
    outfile : str
       The name of the ALS output file
    phot : astropy table
       The input photometry structure.  This should have
         the following data (in the same order):
          ID, X, Y, MAG, ERR, SKY, NITER, CHI, SHARP
    head : list
       The ALS header as list.

    Returns
    -------
    The ALS file to "outfile"

    Example
    -------

    writeals('temp.als',phot,head)

    By D.Nidever  September 2007
    Translated to Python by D. Nidever,  April 2022
    """

    if type(outfile) is not str:
        raise ValueError('outfile must be a string')
    if type(head) is not list:
        raise ValueError('head must be a list')
    if np.array(head).size != 2:
        raise valueError('head must be a 2-element list')

    # Opening the output file
    if os.path.exists(outfile): os.remove(outfile)
    f = open(outfile,'w')

    # Print header
    f.write(head[0]+'\n')
    f.write(head[1]+'\n')
    f.write(' \n')

    # NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
    #  1  2046  4094   114.2 38652.0   11.70    3.00    3.91    1.55    7.02
    # 
    #     10 1476.102   34.512   16.313   0.0318  161.060       4.    0.954   -0.244
    #     24  461.950   55.882   15.043   0.0127  160.980       4.    1.048   -0.372
    
    # Setting up the command
    # From allstar.f
    #            STRNG1 = RNDOFF(XC(I), 9, 3)
    #            STRNG2 = RNDOFF(YC(I), 9, 3)
    #            STRNG3 = RNDOFF(SKY(I), 9, 3)
    #            WRITE (1,321) ID(I), STRNG1, STRNG2, MAG(I), ERR, STRNG3,
    #     .           FLOAT(NITER), CHI(I), SHARP
    #  321       FORMAT (I7, 2A9, F9.3, F9.4, A9, F9.0, 2F9.3)
    numstars = len(phot)
    
    fmt = '%7d%9.3f%9.3f%9.3f%9.4f%9.3f%9s%9.3f%9.3f\n'
    for i in range(numstars):
        data = (phot['id'][i],phot['x'][i],phot['y'][i],phot['mag'][i],phot['err'][i],phot['sky'][i],
                str(phot['niter'][i])+'.',phot['chi'][i],phot['sharp'][i])
        f.write(fmt % data)

    # Closing the file
    f.close()

def loadtrans(transfile):
    pass
 
def loadtfr(tfrfile):
    """
    This loads a DAOMATCH/DAOMASTER tfr file. 
 
    Parameters
    ----------
    tfrfile : str
       The TFR filename 
 
    Returns
    -------
    files : list
      The list of files in the TFR file 
    tab : astropy table
       The structure with final ID, X, Y and 
         the index array.  The indices are 
           FORTRAN/IRAF 1-index based. 
 
    Example
    -------

    files,tab = loadtfr('ccd1001.tr')
 
    By D.Nidever   August 2016 
    Translated to Python by D. Nidever,  April @022
    """
     
    # Test the file 
    if os.path.exists(tfrfile)==False:
        raise ValueError(tfrfile+' NOT FOUND')
     
    # A TFR file looks like this: 
    # F1-00423034_01.als              99.9999   9.9999 
    # F1-00423033_01.als              99.9999   9.9999 
    # F1-00423035_01.als              99.9999   9.9999 
    # F1-00423036_01.als              99.9999   9.9999 
    # F1-00423037_01.als              99.9999   9.9999 
    # F1-00423038_01.als              99.9999   9.9999 
    # F1-00423039_01.als              99.9999   9.9999 
    # F1-00423040_01.als              99.9999   9.9999 
    # F1-00423041_01.als              99.9999   9.9999 
    # F1-00423042_01.als              99.9999   9.9999 
    # F1-00423043_01.als              99.9999   9.9999 
    # F1-00423044_01.als              99.9999   9.9999 
    # F1-00423045_01.als              99.9999   9.9999 
    # F1-00423046_01.als              99.9999   9.9999 
    # F1-00423047_01.als              99.9999   9.9999 
    # F1-00423048_01.als              99.9999   9.9999 
    # ============================== 
    #      1 -1037.98 2452.949      0      0      0      0      0      0    340      0      0      0      0      0      0      0      0      0 
    #      2 -1037.67 3505.380      0      0      0      0      0      0   1222      0      0      0      0      0      0      0      0      0 
    #      3 -1036.54 3174.594      0      0      0      0      0      0    448      0      0      0      0      0      0      0      0      0 
    #      4 -1035.85 4321.116      0      0      0      0      0      0   1263      0      0      0      0      0      0      0      0      0 
    #      5 -1035.28 5458.115      0      0      0      0      0      0    729      0      0      0      0      0      0      0      0      0 
    #      6 -1033.22 2134.540      0      0      0      0      0      0    838      0      0      0      0      0      0      0      0      0 
    #      7 -1032.40 3823.881      0      0      0      0      0      0   1126      0      0      0      0      0      0      0      0      0 
    #      8 -1031.18 5946.214      0      0      0      0      0      0   1075      0      0      0      0      0      0      0      0      0 
    #      9 -1030.97 5823.931      0      0      0      0      0      0      0      0   1773      0      0      0      0      0      0      0 
    #     10 -1030.16 5403.574      0      0      0      0      0      0    725      0   2157      0      0      0      0      0      0      0 
    #     11 -1029.83 4989.178      0      0      0      0      0      0      0      0   2110      0      0      0      0      0      0      0 
    #     12 -1029.31 5322.905      0      0      0      0      0      0      0      0    700      0      0      0      0      0      0      0 
    #     13 -1029.17 3798.451      0      0      0      0      0      0      0      0    377      0      0      0      0      0      0      0 
     
    # Read in the information
    tfrlines = dln.readlines(tfrfile)
     
    # Find the break in in the list 
    brkind = dln.grep(tfrlines,'====',index=True)
    if len(brkind) == 0:
        raise ValueError('ERROR. No ===== break line')
     
    # Filenames
    filelines = tfrlines[0:brkind[0]]
    files = [l.split()[0] for l in filelines]
    nfiles = len(files) 
     
    # Transformation info 
    tlines = tfrlines[brkind[0]+1:] 
    ntlines = len(tlines)
    #fieldwidths = tuple([7,9,9]+(2*nmag+2)*[9])
    #fieldtypes = tuple(['d','f','f']+(2*nmag+2)*['f'])
    #parser = make_parser(fieldwidths,fieldtypes)
    dum = tlines[0].split()
    ncol = len(dum)
    allid = np.zeros(ntlines,int)
    allx = np.zeros(ntlines,float)
    ally = np.zeros(ntlines,float)
    index = np.zeros((ncol-3,ntlines),int)
    nindexcol = index.shape[0]
    for i,l in enumerate(tlines):
        arr = l.split()
        allid[i] = arr[0]
        allx[i] = arr[1]
        ally[i] = arr[2]
        index[:,i] = arr[3:]
     
    # Get unique star IDs and X/Y values
    dum,ui = np.unique(allid,return_index=True)
    ui = np.sort(ui)
    uid = allid[ui]
    ux = allx[ui] 
    uy = ally[ui] 
    nstars = len(uid) 
    # Create structure and stuff in ID, X, Y 
    dt = [('id',int),('x',float),('y',float),('index',int,nfiles)]
    tab = np.zeros(nstars,dtype=np.dtype(dt))
    tab['id'] = uid 
    tab['x'] = ux 
    tab['y'] = uy 
     
    # Load the INDEX information 
    nlinesperstar = ntlines // nstars 
    
    # Load the index information 
    uindex = np.zeros((nfiles,nstars),int)
    for i in range(nstars):
        # Pull out the Nlinesperstar from INDEX for this star
        index1 = index[:,i*nlinesperstar:(i+1)*nlinesperstar]
        # Reformat from 2D to 1D
        index2 = index1.reshape(nlinesperstar*nindexcol)
        # Take only the lines we need 
        #   there might be extra blank elements
        uindex[:,i] = index2[0:nfiles].astype(int)
    tab['index'] = uindex.T

    return files,tab


def daophot_imprep(fluxfile,maskfile,header=False):
    """
    Use the CP Instcal flux and mask files to prepare 
    an image for DAOPHOT 
 
    Parameters
    ----------
    fluxfile : str 
      The filename for the CP Instcal flux file. 
    maskfile : str
       The filename for the CP Instcal mask file. 
    header : boolean, optional
       Return the header only. 
 
    Returns
    -------
    im : numpy array
      The DAOPHOT-prepared image array 
    meta : header
      The header/metadata for IM. 
 
    Example
    -------
    
    im,meta = daophot_imprep(fluxfile,maskfile)
 
    By D. Nidever  Feb 2019 
    Translated to Python by D. Nidever,  April 2022
    """

    nodiffmaskflag = None
 
    # Create the DAOPHOT file 
    #   taken from smashred_imprep_single.pro 
    if header==False:
        fim,fhead = fits.getdata(fluxfile,header=True)
        mim,mhead = fits.getdata(maskfile,header=True)
    else:
        fhead = fits.getheader(fluxfile)
    ccdnum = fhead['CCDNUM']
    
    # --- Prepare the header --- 
         
    # add gain, rdnoise, saturation 
    meta = fhead.copy()
    if 'XTENSION' in meta:
        meta[0]='SIMPLE  =                    T / file does conform to FITS standard             ' 
         
    #gain = (arr[ccd-1].gaina+arr[ccd-1].gainb)*0.5 
    #rdnoise = (arr[ccd-1].rdnoisea+arr[ccd-1].rdnoiseb)*0.5 
    #gain = sxpar(fhead,'ARAWGAIN')
    gainA = fhead['GAINA']
    gainB = fhead['GAINB']
    gain = (gainA+gainB)*0.5
    rdnoiseA = fhead['RDNOISEA']
    rdnoiseB = fhead['RDNOISEB']
    rdnoise = (rdnoiseA+rdnoiseB)*0.5
    meta['gain'] = gain
    meta['rdnoise'] = rdnoise

    # REMOVE DUPLICATE KEYWORDS!!  They cause lots of annoying errors 
    # EXTVER, CHECKSUM, DATASUM
    if 'extver' in meta:
        del meta[bd[1:]]
    if 'checksum' in meta:
        del meta[bd[1:]]
    if 'datasum' in meta:
        del meta[bd[1:]]
         
    # Add "COMMENT " before "BEGIN EXTENSION HEADER ---", it causes problems in daophot
    if 'BEGIN' in meta:
        beg = meta['BEGIN']
        del meta['BEGIN']
        meta['COMMENT'] = beg
         
    if header:
        return header 

    # --- Prepare the image --- 
    med1 = np.median(fim,axis=1) 
    med1slp = np.diff(med1) 
         
    im = np.copy(fim)
         
    # DAOPHOT cannot handle DOUBLE arrays (BITPIX=-64)
    if meta['bitpix'] == -64:
        meta['bitpix'] = -32
        im = im.astype(np.float32)
         
    # Check for differences in amp background levels 
    med1 = np.median(im[:,800:1023]) 
    med2 = np.median(im[:,1024:1200]) 
    err1 = dln.mad(im[:,800:1023])/np.sqrt(len(im[:,800:1023])) 
    err2 = dln.mad(im[:,1024:1200])/np.sqrt(len(im[:,1024:1200])) 
    err = np.sqrt(err1**2 + err2**2) 
         
    # Set bad pixels to saturation value 
    # --DESDM bit masks (from Gruendl): 
    # BADPIX_BPM 1          /* set in bpm (hot/dead pixel/column)        */ 
    # BADPIX_SATURATE 2     /* saturated pixel                           */ 
    # BADPIX_INTERP 4 
    #     /* interpolated pixel                        */ 
    # BADPIX_LOW     8      /* too little signal- i.e. poor read         */ 
    # BADPIX_CRAY   16      /* cosmic ray pixel                          */ 
    # BADPIX_STAR   32      /* bright star pixel                         */ 
    # BADPIX_TRAIL  64      /* bleed trail pixel                         */ 
    # BADPIX_EDGEBLEED 128  /* edge bleed pixel                          */ 
    # BADPIX_SSXTALK 256    /* pixel potentially effected by xtalk from super-saturated source */ 
    # BADPIX_EDGE   512     /* pixel flagged to exclude CCD glowing edges */ 
    # BADPIX_STREAK 1024    /* pixel associated with satellite (airplane/meteor) streak     */ 
    # BADPIX_FIX    2048    /* a bad pixel that was fixed                */ 
    # --CP bit masks, Pre-V3.5.0 (PLVER) 
    # Bit   DQ Type  PROCTYPE 
    # 1  detector bad pixel          InstCal 
    # 1  detector bad pixel/no data  Resampled 
    # 1  No data                     Stacked 
    # 2  saturated                   InstCal/Resampled 
    # 4  interpolated                InstCal/Resampled 
    # 16  single exposure cosmic ray InstCal/Resampled 
    # 64  bleed trail                InstCal/Resampled 
    # 128  multi-exposure transient  InstCal/Resampled 
    # --CP bit masks, V3.5.0 on (after ~10/28/2014), integer masks 
    #  1 = bad (in static bad pixel mask) 
    #  2 = no value (for stacks) 
    #  3 = saturated 
    #  4 = bleed mask 
    #  5 = cosmic ray 
    #  6 = low weight 
    #  7 = diff detect 
    # You can't have combinations but the precedence as in the order 
    # of the list (which is also the order in which the processing 
    # discovers them).  So a pixel marked as "bad" (1) won't ever be 
    # flagged as "diff detect" (7) later on in the processing. 
    # 
    # "Turn off" the "difference image masking", clear the 8th bit 
    # 128 for Pre-V3.5.0 images and set 7 values to zero for V3.5.0 or later. 
    if nodiffmaskflag is None:  # set by default 
        nodiffmaskflag = 1 
    if nodiffmaskflag:
        plver = fhead.get('plver')  # DESDM doesn't have this 
        if plver is not None and plver[0:3] != 'DES': # CP data 
            # DES, didn't have this flag, so skip
                 
            # V3.5.0 and on, Integer masks
            versnum = np.array(plver[1:].split('.')).astype(int)
            if versnum[0] > 3 or (versnum[0] == 3 and versnum[1] >= 5): 
                bdpix = (mim == 7)
                if np.sum(bdpix) > 0: 
                    mim[bdpix] = 0 
            # Pre-V3.5.0, Bitmasks 
            else: 
                bdpix = ( (mim & 2**7) > 0)
                if np.sum(bdpix) > 0:  # clear 128 
                    mim[bdpix] -= 128 
         
    # Add background back in for DES SV data
    skybrite = meta.get('skybrite')
    skysigma = meta.get('skybrite')
    bunit = meta.get('bunit')
    if bunit is None:
        bunit = 'adu'
    else:
        bunit = bunit.strip()
    if bunit == 'electrons' and skybrite is not None:
        saturate = meta.get('saturate')
        if saturate is None:
            saturate = 65000.0 
        gdpix = ((mim==0) & (im<saturate))
        medim = np.median(im[gdpix]) 
        if medim < 100: 
            # Add sky background back in so DAOPHOT can create a proper 
            # noise model 
            # SKYBRITE is given in electrons as well 
            if skybrite != 0: 
                im[gdpix] += skybrite 
            else: 
                # Sometimes SKYBRITE=0, in that case use SKYSIGMA^2 
                #  measure sigma ourselves if SKYSIGMA is not given 
                if skysigma is not None and skysigma > 0: 
                    im[gdpix] += skysigma**2 
                else: 
                    im[gdpix] += dln.mad(im[gdpix])**2 
         
    # Mask bad half of DECam chip 31
    dateobs = meta.get('DATE-OBS')
    if dateobs is not None:  # assume bad if no date 
        t = Time(dateobs, format='fits')
        mjd = t.mjd
    else: 
        mjd = 58000
    if meta.get('INSTRUME')=='DECam' and meta.get('CCDNUM')==31 and mjd>56660:
        print('Masking bad half of DECam chip 31')
        # X: 1-1000 okay 
        # X: 1000-2049 bad 
        fim[:,1000:] = 1e6 
         
    # Set saturated pixels to 65000.0 
    if bunit != 'electrons':
        saturate = meta.get('saturate')
        if saturate is not None:   # set it slightly lower than 65000 for DAOPHOT 
            saturate = np.minimum(saturate,64500.0)
        else: 
            saturate = 64500.0
        meta['saturate'] = saturate
        bdpix = ((mim > 0.0) | (fim > 65000.0))
        if np.sum(bdpix) > 0: 
            im[bdpix] = 65000.0 
    # Allow saturation value to be larger for DES SV data in electrons 
    else: 
        saturate = meta.get('saturate')
        if saturate is None:
            saturate = 64500.0   # set it slightly lower than 65000 for DAOPHOT 
            bdpix = ((mim > 0.0) | (fim > 65000.0))
            if np.sum(bdpix) > 0: 
                im[bdpix] = 65000.0
            meta['saturate'] = saturate
        else: 
            bdpix = ((mim > 0.0) | (fim > saturate))
            if np.sum(bdpix) > 0:  # set slightly higher for DAOPHOT 
                im[bdpix] = saturate*1.01 
            meta['saturate'] = saturate 
         
    return im,meta


def loadresource(rfile):
    """ Load resource file"""
    rlines = dln.readlines(rfile,comment='#',noblank=True)
    arr = [l.split('=') for l in rlines]
    names = [a[0].strip() for a in arr]
    vals = [a[1].strip() for a in arr]
    rstr = dict(zip(names,vals))
    return rstr
 
def fits_read_resource(filename,header=False,nowrite=False):
    """
    Read a FITS file by using the resource file. 
 
    Parameters
    ----------
    filename : str
      The FITS filename to read.  This is the actual file, not 
        the resource file. 
    header : boolean, optional
      Return the header only.  The header will be returned 
        as the output of the function. 
    nowrite : boolean, optional
      Only return the data and don't write the file to the disk. 
 
    Returns
    -------
    im : numpy array
      The 2D FITS image. 
    meta : header
      The header array for IM. 
 
    Example
    -------

    im,meta = fits_read_resource(file,header=header,nowrite=nowrite) 
 
    By D. Nidever  Feb 2019 
    Translated to Python by D. Nidever,  April 2022
    """
     
    t0 = time.time() 
     
    fdir = os.path.abspath(os.path.dirname(filename))
    base = os.path.basename(filename) 
    rfile = fdir+'/.'+base 
     
    # Check if there's a lock file
    lockexists = os.path.exists(filename+'.lock')
    if os.path.exists(filename+'.lock') and (t0-os.path.getmtime(filename+'.lock'))<100:
        print('Lock file.  Waiting')
        while (os.path.exists(filename+'.lock') and time.time()-t0 <100.): time.sleep(5)
         
         
    # If the file exists and is not empty then just use it's data 
    #  Might need to use resource header though 
    #============================================================ 
    if os.path.exists(filename) and os.path.getsize(filename)>1:
        # Check if there's a resouce header we should use 
        # It has a resource file
        rmeta = []
        if os.path.exists(rfile): 
            # Load the resource file
            rstr = loadresource(rfile)
            # There is a resource header 
            if 'header' in rstr.keys():
                if rstr['header'][0] == '/': 
                    hfile = rstr['header']
                else: 
                    hfile = fdir+'/'+rstr['header']
                if os.path.exists(hfile)==False:
                    raise ValueError('Local header file '+hfile+' NOT FOUND')
                else: 
                    rmeta = dln.readlines(hfile)
        # Load the regular FITS file header
        meta = fits.getheader(filename)
             
        # Decide which header to use 
        if len(rmeta) > 0: 
            #print,'Using resource header ',hfile 
            meta = rmeta 
             
        # Only return the header
        if header:
            return meta 
        # Load the data in the fits file 
        im = fits.getdata(filename)
        return im,meta
         
    # If the file is empty then use the resource information 
    #======================================================= 
         
    # Create lock file
    if header==False:
        dln.touch(filename+'.lock')
         
    # Load the resource file
    rstr = loadresource(rfile)
         
    # Only the header, local 
    # get it from the resource file or a stand-alone file
    if header and 'header' in rstr.keys():
        if rstr['header'][0] == '/': 
            hfile = rstr['header']
        else: 
            hfile = fdir+'/'+rstr['header']
        if os.path.exists(hfile) == False:
            raise ValueError('Local header file '+hfile+' NOT FOUND')
        else:
            meta = dln.readlines(hfile)
            return meta
         
         
    # Create a temporary directory to work in
    if header==False:
        tmpdir = tempfile.mkdtemp(prefix="rsrc")
        os.chmod(tmpdir,0o755)
         
    # Get the flux file 
    #   /mss1/archive/pipe/20170328/ct4m/2014B-0404/c4d_170329_043921_ooi_g_v1.fits.fz[39] 
    lo = rstr['fluxfile'].find('[') 
    hi = rstr['fluxfile'].find(']') 
    fluxfile = rstr['fluxfile'][0:lo]
    fext = rstr['fluxfile'][lo+1:hi]
    if header==False:
        tfluxfile = tmpdir+'/flux.fits'
        utils.file_wait(fluxfile)
        out = subprocess.check_output(['funpack','-E',fext,'-O',tfluxfile,fluxfile],shell=False)
         
    # Construct the header from the extension and the main headers: 
    #                       <required keywords for the extension: 
    #                       XTENSION, BITPIX, 
    #                               NAXIS, ...> 
    #                       BEGIN MAIN HEADER 
    #                       -------------------------------- 
    #                       <PDU header keyword and history less required 
    #                       keywords: 
    #                               SIMPLE, BITPIX, NAXIS, ...> 
    #                       BEGIN EXTENSION HEADER 
    #                       --------------------------- 
    #                       <extension header less required keywords that 
    #                       were 
    #                               placed at the beginning of the header. 
    #                       END 
    # Need PDU header with exposure information
    mhead0 = fits.getheader(fluxfile,0)
    if header:
        ehead0 = fits.getheader(fluxfile,fext)
    else:
        ehead0 = fits.getheader(tfluxfile,0)
    # Required keywords 
    #XTENSION= 'IMAGE   '           /extension type 
    #   SIMPE = T 
    #BITPIX  =                  -32 /bits per data value 
    #NAXIS   =                    2 /number of axes 
    #NAXIS1  =                 2046 / 
    #NAXIS2  =                 4094 / 
    #PCOUNT  =                    0 /Number of group parameters 
    #GCOUNT  =                    1 /Number of groups 
    
    # Start the final header 
    head = fits.Header()
    head['SIMPLE'] = 'T','file does conform to FITS standard'
    head['BITPIX'] = ehead0['BITPIX'],'bits per data value' 
    head['NAXIS'] = ehead0['NAXIS'],'number of data axes' 
    head['NAXIS1'] = ehead0['NAXIS1'],'length of data axis 1' 
    head['NAXIS2'] = ehead0['NAXIS2'],'length of data axis 2' 
    head['PCOUNT'] = 0,'Number of group parameters' 
    head['GCOUNT'] = 1,'Number of groups' 
    # Remove required keywords from the main header 
    mhead = mhead0 
    todel = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','EXTEND','DATASUM','CHECKSUM','END'] 
    for j in range(len(todel)):
        if todel[j] in mhead:
            del mhead[todel[j]]
    # Remove required keywords from the extension header 
    ehead = ehead0 
    todel = ['SIMPLE','XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','EXTEND','DATASUM','CHECKSUM','END'] 
    for j in range(len(todel)):
        if todel[j] in ehead:
            del ehead[todel[j]]
    # Combine them all 
    head['COMMENT'] = ' BEGIN MAIN HEADER ---------------------------------                     '
    head.extend(mhead)
    head['COMMENT'] = 'BEGIN EXTENSION HEADER ----------------------------                     ' 
    head.extend(ehead)
    
    # Fix/remove FPACK keywords 
    if header:
        head['BITPIX'] = head['ZBITPIX']
        head['NAXIS1'] = head['ZNAXIS1']
        head['NAXIS2'] = head['ZNAXIS2']
        todel = ['ZNAXIS','ZNAXIS1','ZNAXIS2','ZIMAGE','ZTILE1','ZTILE2','ZCMPTYPE','ZNAME1','ZVAL1',
                 'ZNAME2','ZVAL2','ZPCOUNt','ZGCOUNT','ZQUANTIZ','ZDITHER0','ZTENSION','TFIELDS',
                 'TTYPE1','TFORM1','TTYPE2','TFORM2','TTYPE3','TFORM3'] 
        for j in range(len(todel)): 
            if todel[j] in head:
                del head[todel[j]] 
             
        # DAOPHOT_IMPREP adds GAIN and RDNOISE 
        if header:
            gainA = head['GAINA'] 
            gainB = head['GAINB']
            gain = (gainA+gainB)*0.5 
            rdnoiseA = head['RDNOISEA'] 
            rdnoiseB = head['RDNOISEB']
            rdnoise = (rdnoiseA+rdnoiseB)*0.5 
            head['GAIN'] = gain 
            head['RDNOISE'] = rdnoise 
         
    # This is done in daophot_imprep.pro but if we are only returning 
    # the header then do it here 
    if header:
        saturate = head.get('saturate')
        if saturate is not None:  # set it slightly lower than 65000 for DAOPHOT 
            saturate = np.minumum(saturate,64500.0)
        else: 
            saturate = 64500.0 
        head['saturate'] = saturate 
         
    # Returning only the header
    if header:
        return head
         
         
    # Now update the fluxfile with this new header 
    tflux = fits.getdata(tfluxfile,0)
    os.remove(tfluxfile)
    fits.HDUList(fits.PrimaryHDU(tflux,head)).writeto(tfluxfile,overwrite=True)
         
    # Get the mask file 
    lo = rstr['maskfile'].find('[') 
    hi = rstr['maskfile'].find(')')
    maskfile = rstr['maskfile'][0:lo]
    mext = rstr['maskfile'][lo+1:hi]
    tmaskfile = tmpdir+'/mask.fits'
    utils.file_wait(maskfile)
    out = subprocess.check_output(['funpack','-E',mext,'-O',tmaskfile,maskfile],shell=False)
         
         
    # Create the DAOPHOT file 
    #   taken from smashred_imprep_single.pro
    im,meta = daophot_imprep(tfluxfile,tmaskfile,header=header)         
         
    # Use the local header 
    #   to write to the file and return 
    if header in rstr.keys():
        if rstr['header'][0]=='/':
            hfile = rstr['header']
        else: 
            hfile = fdir+'/'+rstr['header']
        if os.path.exists(hfile)==False:
            raise ValueError('Local header file '+hfile+' NOT FOUND')
        else:
            meta = dln.readlines(hfile)
              
    # Write new image
    if nowrite==False and header==False:
        fits.HDUList(fits.PrimaryHDU(im,meta)).writeto(filename,overwrite=True)
         
    # Delete the lock file
    if os.path.exists(filename+'.lock'): os.remove(filename+'.lock')
                  
    # Clean up
    for f in [tfluxfile,tmaskfile]:
        if os.path.exists(f): os.remove(f)
         
    dt = time.time()-t0 
    
    # Return header only 
    if header:
        return meta 
         
    return im,meta


def readfile(filename,exten=None,header=False,nowrite=False):
    """
    Generic file reading program for PHOTRED 
 
    Parameters
    ----------
    filename : str
       The name of the file to read from.  The extension 
         and type of the file is used to figure out how 
         to read it. 
    exten : int, optional
       The extension to read 
    header : boolean, optional
       Only return the header/metadata. 
    nowrite : boolean, optional
       Don't write the file if using the resource file. 
 
    Returns
    -------
    results : numpy array or table
       The primary data in the file 
    meta : header or dict
       Meta-data in the 1st extension if it exists. 

    Example
    -------

    cat = photred_readfile('F4-20440011_09.als',meta) 
 
    By D. Nidever, Jan 2019 
    Translated to Python by D. Nidever,  April 2022
    """
     
    # Check that file exists 
    if os.path.exists(filename)==False:
        raise ValueError(filename+' NOT FOUND')

    meta = None
     
    # Is this a FITS file 
    isfits = utils.file_isfits(filename) 
     
    # Break down the file into its components
    exists = os.path.exists(filename)
    fdir = os.path.abspath(os.path.dirname(filename))
    base = os.path.basename(filename)
    ext = os.path.splitext(base)[1][1:]
    if isfits and ext=='gz': # gzipped FITS file 
        ext = 'fits' 
    if isfits and ext=='fz': # fpacked FITS file 
        ext = 'fits' 
    if isfits: # fits file with different extension 
        ext = 'fits' 
     
    # Go through the various options
    if ext=='opt':
        return readopt(filename)
    elif ext=='coo':
        return loadcoo(filename)
    elif ext=='ap':
        return loadaper(filename)
    elif ext=='mch':
        return loadmch(filename)
    elif ext=='raw':
        return loadraw(filename)    # STILL NEED THIS!!
    elif ext=='tfr':
        return loadtfr(filename)
    elif ext=='als' or ext=='alf':
        return loadals(filename)
    elif ext=='mag':
        if isfits: 
            if header:  # only read header
                return fits.getheader(filename)
            phot = Table.read(filename,format='ascii')
            hdu = fits.open(filename)
            nhdu = len(hdu)
            hdu.close()
            if nhdu >= 2 : 
                meta = fits.getdata(filename,2)
            return phot,meta
        else:
            if header: #only read header
                return dln.readlines(filename,nreadline=1)
            phot = loadmag(filename)    # STILL NEED THIS!!!!
            return phot,meta  
    elif ext=='ast':
        if isfits:
            hdu = fitsopen(filename)
            nhdu = len(hdu)
            hdu.close()
            if nhdu >= 2 : 
                meta = fits.getdata(filename,2)
            if header:  # only read header 
                if meta is None:
                    meta = fits.getheader(filename)
                return meta
            phot = Table.read(filename,format='ascii')
            return phot,meta
        else: 
            if header:  # only read header
                meta = dln.readlines(filename,nreadline=1)
                return meta
            phot = Table.read(filename,format='ascii')
            return phot,meta
    elif ext=='trans':
        return read_trans(filename)   # STILL NEED THIS
    elif ext=='phot':
        if isfits:
            hdu = fits.open(filename)
            nhdu = len(hdu)
            hdu.close()
            if nhdu>=2:
                meta = fits.getdata(filename,2)
            if header:  # only read header 
                if meta is None:
                    meta = HEADFITS(filename) 
                return meta 
            phot = Table.read(filename,format='ascii')
            return phot,meta
        else: 
            if header: # only read header
                return dln.readlines(filename,nreadline=1)
            phot = Table.read(filename,format='ascii')
            # We need ID to be STRING not LONG
            if 'ID' in phot.keys():
                phot['ID'] = phot['ID'].astype(str)
            else:
                phot['id'] = phot['id'].astype(str)
            return phot,meta
                
    elif ext=='cmb':
        if isfits:
            hdu = fits.open(filename)
            nhdu = len(hdu)
            hdu.close()
            if nhdu>=2:
                meta = fits.getdata(filename,2)
                if header:  # only read header 
                    if meta is None:
                        meta = fits.getheader(filename) 
                    return meta 
                phot = Table.read(fileame,format='ascii')
                return phot,meta
        else: 
            if os.path.exists(filename+'.meta'):
                meta = Table.read(filename+'.meta',format='ascii')
            if header:  # only read header
                if meta is None:
                    meta = dln.readlines(filename,nreadline=1)
                return meta 
            phot = Table.read(filename,format='ascii')
            return phot,meta
                    
    elif ext=='dered':
        if isfits:
            hdu = fits.open(filename)
            nhdu = len(hdu)
            hdu.close()
            if nhdu>=2:
                meta = fits.getdata(filename,2)
            if header:  # only read header 
                if meta is None:
                    meta = fits.getheader(filename) 
                return meta 
            phot = Table.read(filename,format='ascii')
            return phot,meta
        else: 
            if os.path.exists(filename+'.meta'): 
                meta = Table.read(filename+'.meta',format='ascii')
            if header:  # only read header 
                if meta is None:
                    meta = dln.readlines(filename,nreadline=1)
                return meta 
            phot = Table.read(filename,format='ascii')
            return phot,meta
 
    elif ext=='final':
        if isfits:
            hdu = fits.open(filename)
            nhdu = len(hdu)
            hdu.close()
            if nhdu >= 2 : 
                meta = fits.getdata(filename,2)
            if header: # only read header 
                if meta is None:
                    meta = fits.getheader(filename) 
                return meta 
            phot = Table.read(filename,1,format='ascii')
            return phot,meta
        else: 
            if os.path.exists(filename+'.meta'): 
                meta = Table.read(filename+'.meta',format='ascii')
            if header:  # only read header 
                if meta is None:
                    meta = dln.readlines(filename,nreadline=1)
                return meta 
            phot = Table.read(filename,format='ascii')
            return phot,meta
 
    elif ext=='fits':
        # Load using resource file 
        rfilename = fdir+'/.'+base
        if os.path.exists(rfilename):
            return fits_read_resource(filename,header=header,nowrite=nowrite) 
        # Regular FITS file 
        if exten is not None:
            if header: 
                return fits.getheader(filename,exten) 
            return fits.getdata(filename,exten,header=True)
        else: 
            # Figure out whether to read HDU0 or HDU1
            hdu = fits.open(filename)
            nhdu = len(hdu)
            hd0 = hdu[0].header
            data0 = hdu[0].data
            if nhdu > 1:
                hd1 = hdu[1].header
                data1 = hdu[1].data
            hdu.close()
            # Only HDU0 exists
            if nhdu==1:
                if header:
                    return hd0
                return data0,hd0
            # Both exist, HDU0 has no data but HDU1 does
            if hd0['NAXIS']==0 and hd1['NAXIS']>0:
                if header: 
                    return hd1
                return data1,hd1
            # Both exist, both HDU0 and HDU1 have data (meta-data in HDU1)
            if hd0['NAXIS']>0 and hd1['NAXIS']>0:
                meta = data1
                # HDU1 data is NOT meta-data, use header0
                if meta.dtype != np.str:
                    meta = hd0 
                if header: 
                    return meta
                result = data0
                return result,meta
                
    else:
        raise ValueError('Extension '+ext+' not supported')