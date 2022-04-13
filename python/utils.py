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

def loadsetup():
    pass

def make_parser(fieldwidths):
    # https://stackoverflow.com/questions/4914008/how-to-efficiently-parse-fixed-width-files
    cuts = tuple(cut for cut in accumulate(abs(fw) for fw in fieldwidths))
    pads = tuple(fw < 0 for fw in fieldwidths) # bool flags for padding fields
    flds = tuple(zip_longest(pads, (0,)+cuts, cuts))[:-1]  # ignore final one
    slcs = ', '.join('line[{}:{}]'.format(i, j) for pad, i, j in flds if not pad)
    parse = eval('lambda line: ({})\n'.format(slcs))  # Create and compile source code.
    # Optional informational function attributes.
    parse.size = sum(abs(fw) for fw in fieldwidths)
    parse.fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                                                for fw in fieldwidths)
    return parse


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
    with open(filename,'rb') as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
     
    # This is an ALS file
    arr1 = line1.split()
    if arr1[0]=='NL' and line3.strip()=='': 
        # Figure out the number of columns
        f = open(filename,'rb')
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()        
         
        # First line for the first star
        line4 = f.readline()
        instr = line4 
         
        # Check for continuation lines 
        endflag = 0 
        nstarline = 1 
        continuation = 0 
        while (endflag != 1) and ~eof(unit): 
            line4 = f.readline()
         
            # If there are too many frames/columns then these 
            # go on separate lines and lead with 27 blank spaces 
         
            # This is a continuation line 
            if line4[0:16]=='':
                trial = line4[33:34]
                if trial == ' ': 
                    nspaces = 24 
                else: 
                    nspaces = 25
                instr1 = line4[nspaces:]
                instr += instr1 
                nstarline += 1 
                continuation = 1 
            else:
                endflag = 1 
        f.close()
         
        # Now parts the long line for a single star with a formatted read 
        #fmt = '(I7,2A9,'+strtrim(2*nmag+2,2)+'F9.4)' 
        nmag = (len(instr)-(7+2*9+2*9)) / 9 / 2 
        ncol = nmag*2+5 
         
        # ID  X  Y  MAG1  ERR1  MAG2  ERR2 ...  CHI SHARP 
        #nmag = (ncol-5)/2 
         
        # Stars in this file
        numstar = (dln.numlines(filename)-3 )/nstarline 
         
        # LOADING THE DATA 
        #------------------ 
         
        # nfiles < 12 
        # Each star has data on 1 line 
        #if (nmag lt 12) then begin 
        if (continuation == 0): 
             
            fieldtypes = [3,lonarr(ncol-1)+4] 
            fieldnames = ['ID','X','Y'] 
            for i in np.arange(1,nmag+1): 
                fieldnames = [fieldnames,'MAG'+str(i,2),'ERR'+str(i,2)] 
            fieldnames = [fieldnames,'CHI','SHARP'] 

            phot = Table.read(filename)
            #phot = importascii(filename,fieldtype=fieldtypes,fieldnames=fieldnames,skip=3,/noprint) 
            head = [line1,line2] 
             
            # nfiles >= 12 
            # Each star has data on 2 lines 
            # Most of this code was copied from photcalib.pro 
        else: 
             
            # mastable is where everything is stored, id, x, y, unsolved magnitudes, chi, sharp
            mastable = np.zeros((numstar,2*nmag+5),float)
             
            # Reading in the magnitude file
            f = open(filename,'rb')
            line1 = f.readline()
            line2 = f.readline()
            line3 = f.readline()
            head = [line1,line2] 
             
            # Loop through the stars 
            for j in np.arange(numstar): 
                instr = '' 
                instr1 = ' '
                inline = np.zeros(2*nmag+5,float)
                 
                # Loop through the lines per star 
                for k in np.arange(nstarline):
                    instr1 = f.readline()
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
                 
                # WRITE (1,111) IDMAST(IMASTR), POSIT1, POSIT2, 
                # .            ((DATUM(J,I), J=1,2), I=1,NFRM), SUMCHI, SUMSHP 
                #        111       FORMAT (I7, 2A9, 12F9.4:/ (25X, 12F9.4)) 
                 
                # formatted read
                #num,x,y,mags = struct.unpack("6sx8sx9sx6sx2sx30sx6sx6sx6sx2s", instr.strip())

                from io import StringIO
                names = tuple(['id','x','y']+['mag'+str(i+1) for i in np.arange(2*nmag)]+['chi','sharp'])
                formats = tuple(['i4','f4','f4']+(2*nmag+2)*['f4'])
                out = np.loadtxt(StringIO(line), dtype={'names':names,'formats':formats})
                #num,x,y,mags = np.loadtxt(d, dtype={'names': ('a','b','c','d'),
                #                                    'formats': ('i7','s9','s9','f4')},unpack=True)
                fieldwidths = tuple([7,9,9]+(2*nmag+2)*[9])
                parser = make_parser(fieldwidths)
                out = parser(line)

                line = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789\n'
                fieldwidths = (2, -10, 24)  # negative widths represent ignored padding fields
                parse = make_parser(fieldwidths)
                fields = parse(line)
                print('format: {!r}, rec size: {} chars'.format(parse.fmtstring, parse.size))
                print('fields: {}'.format(fields))
                
                
                fmt = '(I7,2A9,'+str(2*nmag+2,2)+'F9.4)' 
                id = 0
                x=''
                y=''
                all=np.zeros(2*nmag+2,float) 
                reads,instr,id,x,y,all,format=fmt 
                inline[0] = id 
                inline[1] = x 
                inline[2] = y 
                inline[3:2*nmag+5-1] = all 
                 
                # old method, causes problems when there are no spaces between values 
                #reads,instr,inline 
                 
                mastable[j,0:2*nmag+5-1] = inline[0:2*nmag+5-1] 
             
            f.close()
             
            # Now transfer to structure 
            fieldtypes = [3,lonarr(nmag-1)+4] 
            fieldnames = ['ID','X','Y'] 
            for i in np.arange(1,nmag+1): 
                fieldnames = [fieldnames,'MAG'+str(i,2),'ERR'+str(i,2)] 
            fieldnames = [fieldnames,'CHI','SHARP'] 
             
            mastable2 = transpose(mastable) 
            if numstar > 1: 
                phot = ARR2STR(mastable2,fieldnames=fieldnames,fieldtypes=fieldtypes)
            else: 
                # arr2str needs a 2D array, make it [Ncol, 1] 
                phot = ARR2STR(reform(mastable2,len(mastable2),1),fieldnames=fieldnames,fieldtypes=fieldtypes)
                
    # Not a raw/makemag file 
    else: 
        print('This is NOT a MAKEMAG file')
        return None,None

    return phot,head,mastable


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

def file_isfits(filename):
    """
    Check if this file is a FITS file or not. 
 
    Parameters
    ----------
    filename   The name of the file to check. 
 
    Returns
    -------
    return     1 if the file is a FITS file and 0 otherwise. 
 
    Example
    ------
  
    test = file_isfits(filename) 
 
    By D. Nidever, Jan 2019 
    Based partially from is_fits.pro by Dave Bazell 
    Translated to python by D. Nidever, April 2022
    """

    # Does the file exist 
    if os.path.exists(filename)==False:
        return False
     
    # Four possible possibilities: 
    # 1) Regular FITS file (this includes fpacked FITS files) 
    # 2) Gzipped FITS file 
    # 3) ASCII file 
    # 4) Some other binary file 
     
    # Try to open the file normally
    # astropy will catch lots of problems
    #  empty file, corrupted file, not a FITS file
    try:
        hdu = fits.open(filename)
    except:
        return False

    hdu.close()
    return True
