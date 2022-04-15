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
from dlnpyutils import utils as dln
import struct
from itertools import zip_longest
from itertools import accumulate
from io import StringIO

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
                    if isfloat(val):
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
