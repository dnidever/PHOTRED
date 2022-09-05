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
from astropy.wcs import WCS
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning
import warnings
from dlnpyutils import utils as dln
import struct
from itertools import zip_longest
from itertools import accumulate
from io import StringIO
from . import utils,imfwhm

# Ignore these warnings
warnings.simplefilter('ignore', category=AstropyWarning)
#warnings.filterwarnings(action="ignore", message=r'FITSFixedWarning:*')

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

def fileinfo(files):
    """
    Gather information about files like their corner coordinates, etc. 
    This is used by the various tiling-related programs. 
 
    Parameters
    ----------
    files : str or list
      The FITS file names. 
 
    Returns
    -------
    info : table
       The structure with information for each file. 
 
    Example
    -------

    info = fileinfo(files)
 
    By D.Nidever  Jan. 2017 
    Translated to Python by D. Nidever,  April 2022
    """

    nfiles = np.array(files).size
    if nfiles==1 and (type(files)==str or type(files)==np.str or type(files)==np.str_):
        files = [files]

    # Create structure 
    dt = [('file',(np.str,300)),('dir',(np.str,300)),('base',(np.str,100)),('ext',(np.str,10)),('exists',bool),
          ('size',int),('nx',int),('ny',int),('filter',(np.str,50)),('exptime',float),
          ('dateobs',(np.str,30)),('mjd',float),('pixscale',float),('cenra',float),('cendec',float),
          ('vertices_ra',(float,4)),('vertices_dec',(float,4))]
    info = np.zeros(nfiles,dtype=np.dtype(dt))
    # File loop 
    for i in range(nfiles): 
        info['file'][i] = files[i]
        info['dir'][i] = os.path.dirname(os.path.abspath(files[i]))
        info['base'][i] = utils.fitsext(files[i],basename=True)
        info['ext'][i] = utils.fitsext(files[i])
        info['exists'][i] = os.path.exists(files[i])
        if info['exists'][i]==False:
            continue
        if files[i][-7:]=='fits.fz':
            head = readfile(files[i],exten=1,header=Tre)
            # Fix the NAXIS1/NAXIS2 in the header 
            head['NAXIS1'] = head['ZNAXIS1']
            head['NAXIS2'] = head['ZNAXIS2']
        else: 
            head = readfile(files[i],header=True)
        info['nx'][i] = head['NAXIS1'] 
        info['ny'][i] = head['NAXIS2']
        info['file'][i] = files[i] 
        try:
            info['filter'][i] = getfilter(files[i],noupdate=True,silent=True)
        except:
            info['filter'][i] = head['filter']
        info['exptime'][i] = head['exptime']
        info['dateobs'][i] = head['date-obs'] 
        info['mjd'][i] = utils.date2jd(info['dateobs'][i],mjd=True) 
        wcs = WCS(head)
        pcoo = wcs.pixel_to_world([info['nx'][i],info['nx'][i]+1],[info['ny'][i],info['ny'][i]+1])
        pixscale = 3600*pcoo[0].separation(pcoo[1]).deg
        info['pixscale'][i] = pixscale 
        coo = wcs.pixel_to_world(info['nx'][i]//2,info['ny'][i]//2)
        cenra1 = coo.ra.deg
        cendec1 = coo.dec.deg
        info['cenra'][i] = cenra1 
        info['cendec'][i] = cendec1 
        vcoo = wcs.pixel_to_world([0,info['nx'][i]-1,info['nx'][i]-1,0],[0,0,info['ny'][i]-1,info['ny'][i]-1])
        vra = vcoo.ra.deg
        vdec = vcoo.dec.deg
        info['vertices_ra'][i] = vra 
        info['vertices_dec'][i] = vdec

    return Table(info)

def getgain(filename=None,head=None):
    """
    This gets the GAIN information from a FITS files 
 
    Parameters
    ----------
    filename : str, optional
       FITS filename 
    head : header, optional
      Use this header array instead of reading FITS file. 

    Returns
    -------
    gain : float
      The gain in e/ADU 
       If the gain is not found in the header then None
       is returned. 
    keyword : str
      The actual header keyword used. 
 
    Example
    -------
    
    gain,keyword = getgain(filename)
 
    By D.Nidever  May 2008 
    Translated to Python by D. Nidever, April 2022
    """

    if filename is None and head is None:
        raise ValueError('filename or head must be input')
    nfiles = dln.size(filename)
    
    # Can't use input HEAD if multiple fits files input 
    if nfiles > 1:
        head = None
     
    # More than one name input 
    if nfiles > 1:
        gain = np.zeros(nfiles,float)
        keyword = np.zeros(nfiles,(np.str,50))
        for i in range(nfiles): 
            gain1,keyword1 = getgain(filename[i])
            gain[i] = gain1
            keyword[i] = keyword1
        return gain,keyword
     
    # No header input, read from fits file 
    if head is None:
        # Check that the file exists
        if os.path.exists(filename)==False:
            raise ValueError(filename+' NOT FOUND')
        if filename[-7:]=='fits.fz':
            head = readfile(filename,1,header=True)
        else:
            head = readfile(filename,header=True)        
     
    gain = head.get('GAIN')
    egain = head.get('EGAIN')   # imacs 
     
    # Use GAIN 
    if gain is not None:
        keyword = 'GAIN' 
    # Use EGAIN 
    if gain is None and egain is not None:
        gain = egain 
        keyword = 'EGAIN' 
             
    # No GAIN 
    if gain is None and egain is None:
        print('NO GAIN FOUND')
        gain = None
        keyword = None

    return gain,keyword


def getrdnoise(filename=None,head=None):
    """
    This gets the RDNOISE information from a FITS files 
 
    Parameters
    ----------
    filename : str, optional
       FITS filename 
    head : header, optional
      Use this header array instead of reading FITS file. 

    Returns
    -------
    rdnoise : float
      The rdnoise in electrons/read
       If the rdnoise is not found in the header then None
       is returned. 
    keyword : str
      The actual header keyword used. 
 
    Example
    -------
    
    rdnoise,keyword = getrdnoise(filename)
 
    By D.Nidever  May 2008 
    Translated to Python by D. Nidever, April 2022
    """

    if filename is None and head is None:
        raise ValueError('filename or head must be input')
    nfiles = dln.size(filename)
    
    # Can't use input HEAD if multiple fits files input 
    if nfiles > 1:
        head = None
     
    # More than one name input 
    if nfiles > 1:
        rdnoise = np.zeros(nfiles,float)
        keyword = np.zeros(nfiles,(np.str,50))
        for i in range(nfiles): 
            rdnoise1,keyword1 = getrdnoise(filename[i])
            rdnoise[i] = rdnoise1
            keyword[i] = keyword1
        return rdnoise,keyword
     
    # No header input, read from fits file 
    if head is None:
        # Check that the file exists
        if os.path.exists(filename)==False:
            raise ValueError(filename+' NOT FOUND')
        if filename[-7:]=='fits.fz':
            head = readfile(filename,1,header=True)
        else:
            head = readfile(filename,header=True)        

    rdnoise = head.get('RDNOISE')                
    readnois = head.get('READNOIS')   # swope
    enoise = head.get('ENOISE')   #    imacs
    
    # Use RDNOISE
    if rdnoise is not None:
        readnoise = rdnoise
        keyword = 'RDNOISE' 
    # Use READNOIS
    if rdnoise is None and readnois is not None:
        readdnoise = readnois
        keyword = 'READNOIS'
    # Use ERDNOISE
    if rdnoise is None and readnois is None and enoise is not None:
        readnoise = enoise
        keyword = 'ENOISE' 
        
    # No RDNOISE 
    if rdnoise is None and readnois is None and enoise is None:
        print('NO READNOISE FOUND')
        readnoise = None
        keyword = None

    return readnoise,keyword


def getexptime(filename=None,head=None):
    """
    This gets the EXPTIME information from a FITS files 
 
    Parameters
    ----------
    filename : str, optional
       FITS filename 
    head : header, optional
      Use this header array instead of reading FITS file. 

    Returns
    -------
    exptime : float
      The exposure time information in seconds.
       If the exptime is not found in the header then None
       is returned. 
 
    Example
    -------
    
    exptime = getexptime(filename)
 
    By D.Nidever  May 2008 
    Translated to Python by D. Nidever, April 2022
    """

    if filename is None and head is None:
        raise ValueError('filename or head must be input')
    nfiles = dln.size(filename)
    
    # Can't use input HEAD if multiple fits files input 
    if nfiles > 1:
        head = None
     
    # More than one name input 
    if nfiles > 1:
        exptime = np.zeros(nfiles,float)
        for i in range(nfiles): 
            exptime[i] = getexptime(filename[i])
        return exptime
     
    # No header input, read from fits file 
    if head is None:
        # Check that the file exists
        if os.path.exists(filename)==False:
            raise ValueError(filename+' NOT FOUND')
        if filename[-7:]=='fits.fz':
            head = readfile(filename,1,header=True)
        else:
            head = readfile(filename,header=True)        

    exptime = head.get('EXPTIME')                
    
    # No EXPTIME 
    if exptime is None:
        print('NO EXPTIME FOUND')

    return exptime

def getpixscale(filename,head=None):
    """
    Get the pixel scale for an image. 
 
    Parameters
    ----------
    file    FITS filename 
    =head   The image header for which to determine the pixel scale. 
    /stp    Stop at the end of the program. 
 
    Returns
    -------
    scale   The pixel scale of the image in arcsec/pix. 
    =error  The error if one occurred. 
 
    Example
    -------

    scale = getpixscale('ccd1001.fits')
 
    BY D. Nidever   February 2008 
    Translated to Python by D. Nidever,  April 2022
    """

    scale = None  # bad until proven good 

     
    # No header input, read from fits file 
    fpack = 0 
    if head is None:
        # Check that the file exists
        if os.path.exists(filename)==False:
            raise ValueError(filename+' NOT FOUND')

        # Open the file
        #hdu = fits.open(filename)
        
        # Fpack or regular fits
        if filename[-7:]=='fits.fz':
            fpack = 1 
            exten = 1 
        else: 
            fpack = 0 
            exten = 0 
         
        # Read the header
        if head is None:
            head = readfile(filename,exten=exten,header=True) 
         
        # Fix NAXIS1/2 in header 
        if fpack == 1:
            head['NAXIS1'] = head['ZNAXIS1']
            head['NAXIS2'] = head['ZNAXIS2']            

    # Does the image have a SCALE parameter
    if head.get('scale') is not None:
        scale = head['scale']
    # Does the image have a PIXSCALE parameter 
    if scale is None:
        pixscale = head.get('pixscale')
        if pixscale is not None:
            scale = pixscale
    # Does the image have a PIXSCALE1 parameter 
    if scale is None: 
        pixscale1 = head.get('pixscale1')
        if pixscale1 is not None:
            scale = pixscale1
     
    # Try the WCS 
    if scale is None:
        try:
            wcs = WCS(head)
             
            # Get the coordinates for two positions 
            #  separated by 1 pixel 
            #head_xyad,head,0.0,0.0,ra1,dec1,/degree 
            #head_xyad,head,1.0,0.0,ra2,dec2,/degree 
            #dist = sphdist(ra1,dec1,ra2,dec2,/deg)*3600. 
            #scale = dist
            c1 = wcs.pixel_to_world(0,0)
            c2 = wcs.pixel_to_world(1,0)            
            dist = c1.separation(c2).arcsec
            scale = dist
            
            if scale == 0.0: 
                scale = None 
        except:
            raise ValueError('No WCS in header')
                
    # Couldn't determine the pixel scale 
    if scale == None:
        error = 'WARNING! COULD NOT DETERMINE THE PIXEL SCALE' 
        print(error)

    return scale
 
def getfilter(filename=None,setup=None,head=None,numeric=False,noupdate=False,
              silent=False,filtname=None,fold_case=False):
    """
    This gets filter information for an image 
    using the "filter" file. 
    The "short" filter names that are returned 
    are not necessarily "unique".  The filter names 
    "I c6028", "T", "T2" all have a short filter 
    name of "T". 
 
    If a filter is not found in the filter list "filters" 
    then a new short name is created for it (based on the 
    first word in the string) and added to the list. 
 
    Parameters
    ----------
    filename : str, optional
       FITS filename 
    setup : dict
       The setup information contained in the photred.setup file.
    head : header, optional
       Use this header array instead of reading FITS file. 
    numeric : boolean, optional
       Return a numeric value instead of letter.  Default is False.
    filtname : str or list, optional
       Input the filter name explicitly instead of giving the 
              FITS filename. 
    noupdate : boolean, optional
       Don't update the "filters" file if the filter is not found.
         Default is to update.
    fold_case : boolean, optional
       Ignore case. The default is for it to be case sensitive. 
    silent : boolean, optional
       Don't print anything to the screen.  Default is False.
 
    Returns
    -------
    The short filter name is output. 
 
    Example
    -------

    filter = getfilter(fitsfile,numeric=numeric,noupdate=noupdate, 
                                 filtname=filtname,fold_case=fold_case)
 
    By D.Nidever  February 2008 
    Translated to Python by D. Nidever, April 2022
    """

    # This is what the "filters" file looks like:
    # 'M Washington c6007'    M
    # 'M Washington k1007'    M
    # 'M'                     M
    # 'I c6028'               T
    # 'I Nearly-Mould k1005'  T
    # 'T'                     T
    
    nfiles = dln.size(filename)
    nfiltname = dln.size(filtname)
    # Not enough inputs
    if filename is None and filtname is None:
        raise ValueError('Need filename or filtname')
     
    # Can't use input HEAD if multiple fits files or filter names input 
    if (nfiles > 1 or nfiltname > 1): 
        head = None 
     
    # More than one FITS filename input 
    if (nfiles > 1): 
        filters = np.zeros(nfiles,(np.str,50))
        for i in range(nfile): 
            filters[i] = getfilter(filename[i],numeric=numeric,
                                   noupdate=noupdate,fold_case=fold_case) 
        return filters 
     
    # More than one filter name input 
    # ONLY if no FITS filename input 
    if (nfiles == 0 and nfiltname > 1): 
        filters = np.zeros(nfiltname,(np.str,50))
        for i in range(nfiltname): 
            filters[i] = getfilter(filtname=filtname[i],numeric=numeric,
                                   noupdate=noupdate,silent=silent,fold_case=fold_case) 
        return filters 
     
     
    # Does the "filters" file exist? 
    if os.path.exists('filters')==False:
        if setup is None:
            raise ValueError('No setup file')
        scriptsdir = setup['scriptsdir']
        if scriptsdir is None:
            raise ValueError('NO SCRIPTSDIR')
        if os.path.exists('filters'): os.remove('filters')
        shutil.copyfile(scriptsdir+'/filters','./filters')
     
    # Load the filters
    lines = dln.readlines('filters',noblank=True)
    gd, = np.where(np.char.array(lines).strip() != '') 
    if len(gd) == 0: 
        raise ValueError('NO FILTERS')
    lines = np.char.array(lines)[gd]
    longnames = [l.split("'")[1] for l in lines]
    shortnames = [l.split("'")[2].strip() for l in lines]
     
    # Get the filter information from the header 
    if (nfiles > 0): 
        # Header NOT input, read FITS files 
        if head is None:
            # Does it have the ".fits" of ".fits.fz" ending
            ext = os.path.splitext(os.path.basename(filename))[1]
            if ext != '.fits' and filename[-7:] != 'fits.fz': 
                raise ValueError(filename+' IS NOT A FITS FILE')
             
            # Make sure the file exists 
            if os.path.exists(filename)==False:
                raise ValueError(filename+' NOT FOUND')
             
            # Read the header 
            if filename[-7:] == 'fits.fz': 
                head = readfile(filename,exten=1,header=True)
            else: 
                head = readfile(filename,header=True)
         
        # Problem with header
        if head is None:
            if silent==False:
                print(filename+' - Problem loading FITS header')
            return '' 
         
        filtname = head.get('FILTER')
        if filtname is None:
            raise ValueError('NO FILTER INFORMATION IN '+filename+' HEADER')
         
    # Get the filter name from "filtname" 
    else: 
        filtname = str(filtname[0]).strip()
     
     
    # Match filter
    ind, = np.where(np.char.array(longnames).lower()==filtname.lower())    
     
    # Match found 
    if len(ind) > 0:     
        filt = shortnames[ind[0]] 
         
        # Numeric value 
        if numeric: 
            snames,ui = np.unique(shortnames,return_index=True)  # unique short names
            nui = len(ui)
            numnames = (np.arange(nui)+1).astype(str) # numbers for the unique short names 
            gg, = np.where(snames == filt)   # which short name 
            numname = numnames[gg[0]] 
            return numname 
         
        return filt
         
    # No match found 
    else:
        if silent==False:
            print('NO FILTER MATCH')
         
        # Add it to the "filters" file 
        if noupdate==False:
             
            # The IRAF task is called "ccdsubset" 
            ## CCDSUBSET -- Return the CCD subset identifier. 
            ## 
            ## 1. Get the subset string and search the subset record file for the ID string. 
            ## 2. If the subset string is not in the record file define a default ID string 
            ##    based on the first word of the subset string.  If the first word is not 
            ##    unique append a integer to the first word until it is unique. 
            ## 3. Add the new subset string and identifier to the record file. 
            ## 4. Since the ID string is used to generate image names replace all 
            ##    nonimage name characters with '_'. 
            ## 
            ## It is an error if the record file cannot be created or written when needed. 
             
            # Get first word of the ID string 
            newshortname = filtname.split()[0]
             
            # Is this already a "short" filter name
            # string comparison
            ind, = np.where(np.char.array(shortnames).lower()==newshortname.lower())
             
            # Yes, already a "short" name 
            # Append integer until unique 
            if len(ind) > 0: 
                #newshortname = newshortname+'_' 
                # Loop until we have a unique name 
                flag = 0 
                integer = 1 
                while (flag == 0): 
                    sinteger = str(integer)
                    ind, = np.where(np.char.array(shortnames).lower()==(newshortname+sinteger).lower())                    
                     
                    # Unique 
                    if len(ind) == 0: 
                        newshortname = newshortname+sinteger 
                        flag = 1 
                    # Not unique, increment integer 
                    else: 
                        integer += 1
                     
             
            # Make sure the variable is okay 
            #newshortname = IDL_VALIDNAME(newshortname,/convert_all) 
             
            # Make sure it doesn't have any weird characters 
            # such as '*', '!', '$'
            newshortname = newshortname.replace('*','_')
            newshortname = newshortname.replace('!','_')
            newshortname = newshortname.replace('$','_')
            newshortname = newshortname.replace('__','_')
            # Add new filter to the "filters" file 
            newline = "'"+filtname+"'     "+newshortname
            with open('filters','wa') as f:
                f.write(newline)
            #dln.writelines('filters',newline,append=True)
            #WRITELINE,'filters',newline,/append 
            print('Adding new filter name to "filters" list')
            print('Filter ID string:  ',filtname)
            print('Filter short name: ',newshortname)
             
             
            # Numeric value 
            if numeric:
                # Reload filters
                lines = dln.readlines('filters')
                lines = [l.strip() for l in lines]
                lines = np.char.array(lines)
                gd, = np.where(lines != '')
                ngd = len(gd)
                lines = lines[gd]
                longnames = [l.split("'")[1] for l in lines]
                shortnames = [l.split("'")[2].strip() for l in lines]                
                 
                snames,ui = np.unique(shortnames,return_index=True)
                nui = len(ui)
                numnames = (np.arange(nui)+1).astype(str)  # numbers for the unique short names 
                 
                gg, = np.where(snames == newshortname)  # which short name 
                numname = numnames[gg[0]] 
                return numname 
             
        # Don't update 
        else: 
            print('NO FILTER MATCH')
            return '' 


def readsetup(setupdir=None,fake=False,std=False):
    """

    Parameters
    ----------
    setupdir : str, optional
       The directory in which to look for the setup file. 
    fake : boolean, optional
       Read the fakered.setup file. 
    std : boolean, optional
       Read the stdred.setup file. 

    Returns
    -------
    setup : dict
      The setup information.  It is a dictionary.

    Example
    -------

    setup = readsetup()
 
    By D.Nidever  March 2008 
    Translated to Python by D. Nidever,  April 2022
    """
 
    if setupdir is None:  # default setup directory 
        setupdir = os.path.abspath('.')
     
    # Type of setup file 
    setupfile = 'photred'
    if std:
        setupfile = 'stdred'
    if fake:
        setupfile = 'fakered' 
     
    # READ THE SETUP FILE 
    #-------------------- 
    # This is a 2xN array.  First colume are the keywords 
    # and the second column are the values. 
    # Use READPAR.PRO to read it 
    setupfiles = glob(setupdir+'/'+setupfile+'.*setup')
    nsetupfiles = len(setupfiles)
    if (nsetupfiles < 1): 
        raise ValueError('NO '+str(setupfile)+' SETUP FILE')
    if (nsetupfiles > 1):
        raise ValueError('MORE THAN ONE '+str(setupfile)+' SETUP FILE FOUND')
     
    # Read the setup file
    lines = dln.readlines(setupfiles[0],comment='#')
    if lines is not None:
        nlines = len(lines)
    else:
        raise ValueError('setup file is empty')

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

def readmch(mchfile):
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

    files,trans,magoff = readmch('ccd1001.mch')
     
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

def readraw(filename):
    """
    This reads the DAOPHOT/DAOMASTER RAW file and also .mag and .makemag files.
     
    Parameters
    ----------
    filename : str
      The name of the .raw file.
     
    Returns
    -------
    phot : astropy table
      A structure with the raw data.
    head : list
      A string array with the raw header.
     
    Example
    ------
    
    phot,head = readraw('temp.raw')
     
    By D. Nidever   Feb. 2008 (basically a copy of readals.pro 
    Translated to Python by D. Nidever,  April 2022
    """
     
    if os.path.exists(filename)==False:
        raise ValueError(filename+' DOES NOT EXIST')

    base = os.path.basename(filename)
    ext = os.path.splitext(base)[1][1:]
    
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
                trial = line4[33]
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
        if ext=='mag' or ext=='makemag':
            nmag = int(  (len(instr)-(9+2*9+2*9)) / 9 / 2 )            
        ncol = nmag*2+5 
         
        # ID  X  Y  MAG1  ERR1  MAG2  ERR2 ...  CHI SHARP 
        #nmag = (ncol-5)/2 
         
        # Stars in this file
        numstar = int((dln.numlines(filename)-3 )/nstarline)

        # READING THE DATA 
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
        if ext=='mag' or ext=='makemag':
            fieldwidths = tuple([9,9,9]+(2*nmag+2)*[9])            
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
                    trial = instr1[33]
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
        
    # Not a raw file 
    else: 
        print('This is NOT a raw file')
        # try to read as a table
        phot = Table.read(filename,format='ascii')
        return phot

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


# Make meta-data dictionary for an image:
def makemeta(fluxfile=None,header=None):
    '''
    This creates a meta-data dictionary for an exposure that is used by many
    of the photometry programs.  Either the filename or the header must be input.
    Note that sometimes in multi-extension FITS (MEF) files the information needed
    is both in the primary header and the extension header.  In that case it is best
    to combine them into one header and input that to makemeta().  This can easily
    be accomplished like this:
      
       head0 = fits.getheader("image1.fits",0)
       head = fits.getheader("image1.fits",1)
       head.extend(head0,unique=True)
       meta = makemeta(header=head)

    Parameters
    ----------
    fluxfile : str, optional
             The filename of the FITS image.
    header : str, optional
           The header of the image.

    Returns
    -------
    meta : astropy header
        The meta-data dictionary which is an astropy header with additional
        keyword/value pairs added.

    Example
    -------

    Create the meta-data dictionary for `image.fits`

    .. code-block:: python

        meta = makemeta("image.fits")

    Create the meta-data dictionary from `head`.

    .. code-block:: python

        meta = makemeta(header=head)

    '''

    # You generally need BOTH the PDU and extension header
    # To get all of this information

    if (fluxfile is None) & (header is None):
        print("No fluxfile or headerinput")
        return
    # Initialize meta using the header
    if fluxfile is not None:
        header = readfile(fluxfile,header=True,exten=0)
    meta = header

    #- INSTCODE -
    if "DTINSTRU" in meta.keys():
        if meta["DTINSTRU"] == 'mosaic3':
            meta["INSTCODE"] = 'k4m'
        elif meta["DTINSTRU"] == '90prime':
            meta["INSTCODE"] = 'ksb'
        elif meta["DTINSTRU"] == 'decam':
            meta["INSTCODE"] = 'c4d'
        else:
            print("Cannot determine INSTCODE type")
            return
    else:
        print("No DTINSTRU found in header.  Cannot determine instrument type")
        return

    #- RDNOISE -
    if "RDNOISE" not in meta.keys():
        # Check DECam style rdnoise
        if "RDNOISEA" in meta.keys():
            rdnoisea = meta["RDNOISEA"]
            rdnoiseb = meta["RDNOISEB"]
            rdnoise = (rdnoisea+rdnoiseb)*0.5
            meta["RDNOISE"] = rdnoise
        # Check other names
        if meta.get('RDNOISE') is None:
            for name in ['READNOIS','ENOISE']:
                if name in meta.keys(): meta['RDNOISE']=meta[name]
        # Bok
        if meta['INSTCODE'] == 'ksb':
            meta['RDNOISE']= [6.625, 7.4, 8.2, 7.1][meta['CCDNUM']-1]
        if meta.get('RDNOISE') is None:
            print('No RDNOISE found')
            return
    #- GAIN -
    if "GAIN" not in meta.keys():
        try:
            gainmap = { 'c4d': lambda x: 0.5*(x.get('GAINA')+x.get('GAINB')),
                        'k4m': lambda x: x.get('GAIN'),
                        'ksb': lambda x: [1.3,1.5,1.4,1.4][x.get['CCDNUM']-1] }  # bok gain in HDU0, use list here
            gain = gainmap[meta["INSTCODE"]](meta)
            meta["GAIN"] = gain
        except:
            gainmap_avg = { 'c4d': 3.9845419, 'k4m': 1.8575, 'ksb': 1.4}
            gain = gainmap_avg[meta["INSTCODE"]]
            meta["GAIN"] = gain
    #- CPFWHM -
    # FWHM values are ONLY in the extension headers
    cpfwhm_map = { 'c4d': 1.5 if meta.get('FWHM') is None else meta.get('FWHM')*0.27, 
                   'k4m': 1.5 if meta.get('SEEING1') is None else meta.get('SEEING1'),
                   'ksb': 1.5 if meta.get('SEEING1') is None else meta.get('SEEING1') }
    cpfwhm = cpfwhm_map[meta["INSTCODE"]]
    meta['CPFWHM'] = cpfwhm
    #- PIXSCALE -
    if "PIXSCALE" not in meta.keys():
        pixmap = { 'c4d': 0.27, 'k4m': 0.258, 'ksb': 0.45 }
        try:
            meta["PIXSCALE"] = pixmap[meta["INSTCODE"]]
        except:
            w = WCS(meta)
            meta["PIXSCALE"] = np.max(np.abs(w.pixel_scale_matrix))

    return meta


def mkopt(base=None,meta=None,va=1,lo=7.0,th=3.5,ls=0.2,hs=1.0,lr=-1.0,hr=1.0,
          wa=-2,an=-6,ex=5,pe=0.75,pr=5.0,cr=2.5,ce=6.0,ma=50.0,red=1.0,wa2=0.0,
          fitradius_fwhm=1.0,hi=None,rd=None,ga=None,fw=None,logger=None):
    """
    Create the DAOPHOT and ALLSTAR option files (.opt and .als.opt) for an exposure.

    Parameters
    ----------
    base : str
         The base name to use for the option files.  The DAOPHOT option file will
         be called `base`.opt and the ALLSTAR option file `base`.als.opt
    meta : astropy dictionary
         The metal-data dictionary for the image.    
    va : int, default = 1
       The variable type of PSF to use.
       -1: Analytic PSF only
        0: Analytic PSF and look-up table of empirical corrections
        1: linear variations across the field
        2: quadratic variations across the field
    lo : float, default = 7.0
       Low good datum (7. works fine on most imags).
    th : float, default = 3.5
       Threshold in sigma above the background (3.5 works fine).
    ls : float, default = 0.2
       Lower sharpness cutoff.
    hs : float, default = 1.0
       High sharpness cutoff.
    lr : float, default = -1.0
       Lower roundness cutoff.
    hr : float, default = 1.0
       High roundness cutoff.
    wa : int, default = -2
       Watch progress for DAOPHOT.  Determines what output is displayed.
    an : int, default = -6
       Analytic model PSF.
        1: Gaussian (3 pararameters)
        2: Moffat function (3 parameters), beta=1.5
        3: Moffat function (3 parameters), beta=2.5
        4: Lorentz function (3 parameters)
        5: Penny function, Gauss+Lorentz (4 parameters), G+L are parallel
        6: Penny function (5 parameters), G and L can be in different directions
        A negative sign in front means to try all functions up to X and pick the best one.
    ex : int, default = 5
       Extra PSF cleaning passes.
    pe : float, default = 0.75
       Percent error due to the uncertainty in the fine-scale structure of the flat field.
    pr : float, default = 5.0
       Profile error due to the incompleteness of the PSF model.
    cr : float, default = 2.5
       Clipping range.  Used to remove outlier pixels. Parameter "a" in the formula given in
       Stetson 1987, PASP, 99, 191, section III.D.2.d "Resisting bad data".
    ce : float, default = 6.0
       Clipping exponent.  Parameter b in above clipping formula.
    ma : float, default = 50.0
       Maximum group size
    red : float, default = 1.0
        Redetermine centroid (0 = no, 1 = yes).
    wa2 : float, default = 0.0
        Watch progress for ALLSTAR.      
    fitradius_fwhm : float, default = 1.0
        The fitting radius size in units of the seeing FWHM for the area to be fit.
    hi : float, optional
       High good datum.  Normally set by `saturate` from `meta`.
    rd : float, optional
       The read noise in electrons. Normally set by `rdnoise` from `meta`.
    ga : float, optional
       The gain in electrons/ADU. Normally set by `gain` from `meta`.
    fw : float, optional
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

    if logger is None: logger=dln.basiclogger('phot')   # set up basic logger if necessary

    optfile = base+".opt"
    alsoptfile = base+".als.opt"

    if meta is None and (ga is None or rd is None or fw is None or hi is None):
        meta = makemeta(base+'.fits')

    # Get frame specific parameters from meta if necessary
    if ga is None: ga = meta['gain']
    if rd is None: rd = meta['rdnoise']
    if fw is None: fw = meta['fwhm'] / meta['pixscale']
    if hi is None: hi = meta['saturate']


    # Calculating some things
    fw = np.min([ fw , 20 ])            # daophot won't accept anything higher than 20
    re = rd/ga
    fi = np.min([ fitradius_fwhm*FW , 51 ])                  # daophot won't accept anything higher than 51
    ps = np.min([ (4.0*fw) , 51 ])       # daophot won't accept anything higher than 51
    ins = np.min([ (fi - 1.0) , 35 ])     # daophot won't accept anything higher than 35
    os = np.min([ (ps + 1.0) , 100 ])    # daophot won't accept anything higher than 100

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

    outarr = [re,ga,lo,hi,fw,th,ls,hs,lr,hr,wa,fi,ps,va,an,ex,pe,pr]
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

    outarr2 = [fi,ins,os,red,wa2,pe,pr,cr,ce,ma]
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


def readals(filename,silent=False):
    """
    This reads the ALLSTAR photometry file (.als).
 
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
    
    phot,head = readals('obj2153_1.als')
 
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
     

def readcoo(filename,silent=False):
    """
    This reads the DAOPHOT coordinates file (.coo).

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
    
    phot,head = readcoo('obj2153_1.coo')
 
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
     

def readaper(filename,silent=False):
    """
    This reads the DAOPHOT aperture photometry file (.ap).
 
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
    
    phot,head = readaper('obj2153_1.ap')
 
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
         
        # READING THE DATA 
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

    
def readtrans(transfile,silent=False,logfile=None):
    """
    Read in a photometric transformation file. 
 
    Parameters
    ----------
    transfile : str
      This gives the transformation information needed 
        to calibrate the raw photometry.  Normally there 
        is a different one of these for every night. 
 
        There need to be two lines per band. 
        FIRST line:  band name,  color name, transformation 
        equation coefficients (zero point, airmass, color 
        airmass*color, color^2) 
        SECOND line: errors in the transformation equation 
        coefficients 
 
        This is an example transfile: 
        M    M-T  -0.9990    0.1402     -0.1345    0.0000   0.0000 
                  1.094E-02  5.037E-03  2.010E-03  0.0000   0.0000 
        T    M-T  -0.0061    0.0489     0.0266     0.0000   0.0000 
                  6.782E-03  3.387E-03  1.374E-03  0.0000   0.0000 
        D    M-D  1.3251     0.1403     -0.0147    0.0000   0.0000 
                  1.001E-02  5.472E-03  2.653E-02  0.0000   0.0000 
 
        If the transfile has chip information then it should 
        look like this: 
        1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000 
                    0.0040   -0.0000    0.0001   0.0000   0.0000 
 
        2  G  G-R  -0.3617    0.1713   -0.1193   0.0000   0.0000 
                    0.0039   -0.0000    0.0001   0.0000   0.0000 
 
        3  G  G-R  -0.3457    0.1713   -0.1193   0.0000   0.0000 
                    0.0039   -0.0000    0.0001   0.0000   0.0000 
 
        If the transfile has night and chip information then it should 
        look like this: 
        55975  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000 
                           0.0040   -0.0000    0.0001   0.0000   0.0000 
 
        55975  2  G  G-R  -0.3617    0.1713   -0.1193   0.0000   0.0000 
                           0.0039   -0.0000    0.0001   0.0000   0.0000 
 
        55975  3  G  G-R  -0.3457    0.1713   -0.1193   0.0000   0.0000 
                           0.0039   -0.0000    0.0001   0.0000   0.0000 
 
        The transformation information can also be input for 
        individual file: 
        F5-00517150_43  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000 
                                 0.0040   -0.0000    0.0001   0.0000   0.0000 
 
        F5-00517150_44  G  G-R  -0.4284    0.1767   -0.1261   0.0000   0.0000 
                                 0.0030    0.0020    0.0001   0.0000   0.0000 
 
        F5-00517150_45  G  G-R  -0.5082    0.1801   -0.1215   0.0000   0.0000 
                                 0.0025    0.0023    0.0001   0.0000   0.0000 
 
        Each pair of lines can use a different format, i.e. the four format 
        types can be mixed in a file.  But this is not recommended since 
        it might be confusing what transformation to use for a given file 
        since there could be multiple "matches". 
 
    silent : boolean, optional
      Don't print anything to the screen. 
    logfile : str, optional
      The name of a logfile to write messages to. 
 
    Returns
    -------
    trans : dict
      The transformation structure.  NIGHT and CHIP will 
        always be included even if they are "blank". 
 
    Example
    -------

    trans = readtrans('n1.trans')
 
    By D.Nidever  Feb.2013 
    Added Night+chip format March 2015 
    Added filename format  July 2017 
    Translated to Python by D. Nidever,  April 2022
    """
     
    if os.path.exists(transfile) == False:
        raise ValueError(transfile+' NOT FOUND')
    
    # Logfile
    # Set up logging to screen and logfile
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    logger = logging.getLogger() 
    while logger.hasHandlers(): # some existing loggers, remove them   
        logger.removeHandler(logger.handlers[0]) 
    logger = logging.getLogger()
    if logfile is not None:
        if os.path.exists(logfile): os.remove(logfile)
        fileHandler = logging.FileHandler(logfile)
        fileHandler.setFormatter(logFormatter)
        rootLogger.addHandler(fileHandler)
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)
        logger.setLevel(logging.NOTSET)
    else:
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)
        logger.setLevel(logging.NOTSET)   
     
    ## ##################################################### 
    ## READ THE TRANSFORMATION FILE 
    ## Two lines per band. 
    ## First line:  band name,  color name, transformation equation coefficients 
    ## Second line: errors in the transformation equation coefficients 
     
    # If this has chip-specific transformations then the lines will be 
    # First line:  chip,  band name, color name, trans eqns. 
    # second line:  errors 
    #  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000 
    #              0.0040   -0.0000    0.0001   0.0000   0.0000 
    # 
    #  2  G  G-R  -0.3617    0.1713   -0.1193   0.0000   0.0000 
    #              0.0039   -0.0000    0.0001   0.0000   0.0000 
    # 
    #  3  G  G-R  -0.3457    0.1713   -0.1193   0.0000   0.0000 
    #              0.0039   -0.0000    0.0001   0.0000   0.0000 
    # 
    # If the file has night and chip information then the lines will be 
    # First line: night, chip,  band name, color name, trans eqns. 
    # second line:  errors 
    #  55975  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000 
    #                     0.0040   -0.0000    0.0001   0.0000   0.0000 
    # 
    #  55975  2  G  G-R  -0.3617    0.1713   -0.1193   0.0000   0.0000 
    #                     0.0039   -0.0000    0.0001   0.0000   0.0000 
    # 
    #  55975  3  G  G-R  -0.3457    0.1713   -0.1193   0.0000   0.0000 
    #                     0.0039   -0.0000    0.0001   0.0000   0.0000 
    # 
    # If the file has information on individual files then the lines will be 
    # First line: filename,  band name, color name, trans eqns. 
    # second line:  errors 
    #  F5-00517150_43  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000 
    #                           0.0040   -0.0000    0.0001   0.0000   0.0000 
    # 
    #  F5-00517150_44  G  G-R  -0.4284    0.1767   -0.1261   0.0000   0.0000 
    #                           0.0030    0.0020    0.0001   0.0000   0.0000 
    # 
    #  F5-00517150_45  G  G-R  -0.5082    0.1801   -0.1215   0.0000   0.0000 
    #                           0.0025    0.0023    0.0001   0.0000   0.0000 
    # 
     
    # What options? 
    # 1 - single night, single chip (no NIGHT or CHIP information) 
    # 2 - single night, multiple chips (no NIGHT but CHIP information) 
    # 3 - multiple nights, multiple chips (NIGHT and CHIP information) 
    # 4 - individual files 
    optcase = -1 
    # I DON'T THINK "optcase" IS ACTUALLY USED FOR ANYTHING ANYMORE 

    dt = [('night',int),('chip',int),('file',np.str,300),('band',np.str,10),('color',np.str,10),('colband',np.str,10),
          ('colsign',int),('zpterm',float),('amterm',float),('colterm',float),('amcolterm',float),('colsqterm',float),
          ('zptermsig',float),('amtermsig',float),('coltermsig',float),('amcoltermsig',float),('colsqtermsig',float)]
    nlines = dln.numlines(transfile)
    trans = np.zeros(nlines,dtype=np.dtype(dt))
    trans['night'] = -1
    trans['chip'] = -1

    count = 0
    with open(transfile,'r') as f:
        while True:
            line = f.readline()
            if line=='':  # we're at the end
                break
            # Blank or commented line
            if line.strip()=='' or line[0]=='#':
                continue
            line = line.replace('\n','')

            arr = line.split()
            narr = len(arr)
             
            # This is the format.  NIGHT (MJD), CHIP, BAND, COLOR, ZPTERM, AMTERM, 
            #                       COLORTERM, AMCOLTERM, AMSQTERM 
            #  55975  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000 
            #                     0.0040   -0.0000    0.0001   0.0000   0.0000 
            # NIGHT and CHIP are optional. 
            # For specific files it looks like this: 
            #  F5-00517150_43  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000 
            #                           0.0040   -0.0000    0.0001   0.0000   0.0000 
             
            # Are NIGHT and CHIP supplied? 
            # ---NIGHT and CHIP--- 
            if narr==9:         
                # Make sure NIGHT and CHIP are numbers
                if arr[0].isnumeric() and arr[1].isnumeric():
                    trans['night'][count] = arr[0]
                    trans['chip'][count] = arr[1] 
                    arr = arr[2:] 
                # Parsing error 
                else: 
                    raise ValueError('NIGHT and CHIP must be numbers.')
                optcase = np.maximum(optcase,3)
                 
            # ---CHIP or FILENAME--- 
            elif narr==8:     
                # CHIP, first value is a number
                if arr[0].isnumeric():
                    trans['chip'][count] = arr[0]
                    arr = arr[1:] 
                    optcase = np.maximum(optcase,2)
                    # FILENAME, first value is NOT a number 
                else: 
                    trans['file'][count] = arr[0].strip()
                    arr = arr[1:] 
                    optcase = np.maximum(optcase,4)
             
            # ---Neither--- 
            elif narr==7:
                # not much to do 
                optcase = np.maximum(optcase,1)
 
            # ---Not enough values--- 
            else:
                raise ValueError('Need at least 7 values in the TRANS line.')
 
            # Parse the rest of the line 
            trans['band'][count] = arr[0] 
            trans['color'][count] = arr[1] 
            trans['zpterm'][count] = arr[2] 
            trans['amterm'][count] = arr[3] 
            trans['colterm'][count] = arr[4] 
            trans['amcolterm'][count] = arr[5] 
            trans['colsqterm'][count] = arr[6] 
 
            # Reading in the error line 
            #---------------------------
            line2 = f.readline().replace('\n','')
            arr2 = line2.split()
            narr2 = len(arr2) 
 
            # Need at least 5 terms 
            if narr2 < 5: 
                raise ValueError('Need at least 5 values in the ERROR line.')
 
            # Parse the error line 
            trans['zptermsig'][count] = arr2[0] 
            trans['amtermsig'][count] = arr2[1] 
            trans['coltermsig'][count] = arr2[2] 
            trans['amcoltermsig'][count] = arr2[3] 
            trans['colsqtermsig'][count] = arr2[4] 

            count += 1
 

    # Trim extra elements
    trans = trans[0:count]
            
    # Leave in NIGHT and CHIP columns even if they are "blank".  For 
    # consistency sake 
 
    # No chip information, strip CHIP 
    #gdnight = where(trans.chip ge 0,ngdchip) 
    #gdchip = where(trans.chip ge 0,ngdchip) 
    #if ngdchip eq 0 then begin 
    #  oldtrans = trans 
    #  trans = replicate({band:'',color:'',colband:'',colsign:0,zpterm:0.0d,amterm:0.0d,colterm:0.0d,
    #            amcolterm:0.0d,colsqterm:0.0d,zptermsig:0.0d,amtermsig:0.0d,coltermsig:0.0d,
    #            amcoltermsig:0.0d,colsqtermsig:0.0d},n_elements(trans)) 
    #  STRUCT_ASSIGN,oldtrans,trans 
    #  undefine,oldtrans 
    #endif 
 
    ntrans = len(trans) 
 
    # Figure out the colband and colsign for each band/chip 
    for i in range(ntrans): 
        band = trans['band'][i].strip()
        col = trans['color'][i]
        col = "".join(col.split())  # remove any spaces
        # Splitting up the two bands
        arr = np.char.array(col.split('-'))
        ind, = np.where(arr==band)
        # colsign = 1 means band - colband 
        if (ind[0] == 0): 
            trans['colband'][i] = arr[1] 
            trans['colsign'][i] = 1 
        # colsign = 2 means colband - band 
        elif (ind[0] == 1): 
            trans['colband'][i] = arr[0] 
            trans['colsign'][i] = 2 
        elif (ind[0] == -1): 
            trans['colband'][i] = '' 
            trans['colsign'][i] = -1 
 
    # Print the transformation equations
    if silent==False:
        logger.info(' TRANSFORMATION EQUATIONS')
        logger.info('--------------------------------------------------------------------------------')
        logger.info('  NIGHT/CHIP/FILE  BAND COLOR ZERO-POINT  AIRMASS   COLOR     AIR*COL   COLOR**2 ')
        logger.info('--------------------------------------------------------------------------------')
        for i in range(ntrans):
            form1 = '%10d%6d%6s%7s%10.4f%10.4f%10.4f%10.4f%10.4f'
            form1f = '%-16s%6s%7s%10.4f%10.4f%10.4f%10.4f%10.4f'
            # FILENAME 
            if trans['file'][i] != '': 
                logger.info(form1f % (trans['file'][i],'  '+trans['band'][i],trans['color'][i],trans['zpterm'][i],
                                     trans['amterm'][i],trans['colterm'][i],trans['amcolterm'][i],trans['colsqterm'][i]))
            # NO Filename 
            else:
                logger.info(form1 % (trans['night'][i],trans['chip'][i],'  '+trans['band'][i],trans['color'][i],trans['zpterm'][i],
                                     trans['amterm'][i],trans['colterm'][i],trans['amcolterm'][i],trans['colsqterm'][i]))
            form2 = '%29s%10.4f%10.4f%10.4f%10.4f%10.4f'
            logger.info(form2 % ('',trans['zptermsig'][i],trans['amtermsig'][i],trans['coltermsig'][i],
                                 trans['amcoltermsig'][i],trans['colsqtermsig'][i]))
        logger.info('--------------------------------------------------------------------------------')
        logger.info('')
 
    return trans
 
def readtfr(tfrfile):
    """
    This reads a DAOMATCH/DAOMASTER tfr file. 
 
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

    files,tab = readtfr('ccd1001.tr')
 
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
    tab = Table(tab)

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
    #if 'extver' in meta:
    #    del meta[bd[1:]]
    #if 'checksum' in meta:
    #    del meta[bd[1:]]
    #if 'datasum' in meta:
    #    del meta[bd[1:]]
         
    # Add "COMMENT " before "BEGIN EXTENSION HEADER ---", it causes problems in daophot
    #if 'BEGIN' in meta:
    #    beg = meta['BEGIN']
    #    del meta['BEGIN']
    #    meta['COMMENT'] = beg
         
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
            versnum = np.array(plver[1:].split('.'))
            if int(versnum[0]) > 3 or (int(versnum[0]) == 3 and int(versnum[1]) >= 5): 
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


def readresource(rfile):
    """ Read resource file"""
    rlines = dln.readlines(rfile,comment='#',noblank=True)
    arr = [l.split('=') for l in rlines]
    names = [a[0].strip().lower() for a in arr]
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
        rmeta = None
        if os.path.exists(rfile): 
            # Read the resource file
            rstr = readresource(rfile)
            # There is a resource header 
            if 'header' in rstr.keys():
                if rstr['header'][0] == '/': 
                    hfile = rstr['header']
                else: 
                    hfile = fdir+'/'+rstr['header']
                if os.path.exists(hfile)==False:
                    raise ValueError('Local header file '+hfile+' NOT FOUND')
                else: 
                    rmeta = fits.Header.fromfile(hfile,sep='\n',endcard=False,padding=False)
        # Read the regular FITS file header
        meta = fits.getheader(filename)
             
        # Decide which header to use 
        if rmeta is not None:
            #print,'Using resource header ',hfile 
            meta = rmeta 
             
        # Only return the header
        if header:
            return meta 
        # Read the data in the fits file 
        im = fits.getdata(filename)
        return im,meta
         
    # If the file is empty then use the resource information 
    #======================================================= 
         
    # Create lock file
    if header==False:
        dln.touch(filename+'.lock')
         
    # Read the resource file
    rstr = readresource(rfile)
         
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
            meta = fits.Header.fromfile(hfile,sep='\n',endcard=False,padding=False)
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
    fits.PrimaryHDU(tflux,head).writeto(tfluxfile,overwrite=True)
         
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
    if 'header' in rstr.keys():
        if rstr['header'][0]=='/':
            hfile = rstr['header']
        else: 
            hfile = fdir+'/'+rstr['header']
        if os.path.exists(hfile)==False:
            raise ValueError('Local header file '+hfile+' NOT FOUND')
        else:
            meta = fits.Header.fromfile(hfile,sep='\n',endcard=False,padding=False)

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
        return readcoo(filename)
    elif ext=='ap':
        return readaper(filename)
    elif ext=='mch':
        return readmch(filename)
    elif ext=='raw' or ext=='makemag':
        return readraw(filename)
    elif ext=='tfr':
        return readtfr(filename)
    elif ext=='als' or ext=='alf':
        return readals(filename)
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
            if header: # only read header
                return dln.readlines(filename,nreadline=1)
            phot = readraw(filename)
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
        return readtrans(filename)
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
        # Read using resource file 
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
