#!/usr/bin/env python
#
# Script for removing saturated stars (with some saturated pixels) from a DAOPHOT list

import os
import sys
import numpy as np
from astropy.io import fits

def read(filename):
    with open(filename) as fp:
        contents = fp.read()
    return contents

def write(filename,lines):
    if os.path.exists(filename): os.remove(filename)
    # Write the file
    f = open(filename,'w')
    f.writelines(lines)
    f.close()

# Main command-line program
if __name__ == "__main__":
    if len(sys.argv)<2:
        print('lstfilter.py lstfilter [outfile]')
        sys.exit()
    # List filename
    filename = sys.argv[1]
    print('Input list '+filename)
    base,ext = os.path.splitext(os.path.basename(filename))
    if len(sys.argv)>2:
        outfile = sys.argv[2]
    else:
        outfile = filename
    # Load the list file
    lines = read(filename)
    lines = lines.split('\n')
    nlines = len(lines)
    if nlines<4:
        raise Exception('No stars in '+filename)
    slines = lines[3:]
    if slines[-1]=='': del slines[-1]   # remove blank last line
    if type(slines) is list:
        nstars = len(slines)
    else:
        nstars = 1
    # Load the coo file and sharp values
    if os.path.exists(base+'.coo') is False:
        raise Exception(base+'.coo NOT FOUND')
    coolines = read(base+'.coo')
    coolines = coolines.split('\n')
    if coolines[-1]=='': del coolines[-1] # remove blank last line 
    if coolines[-1]=='': del coolines[-1] # remove blank last line
    if len(coolines)<4:
        raise Exception('No stars in '+base+'.coo')
    scoolines = coolines[3:]
    sharpdict = {}
    for l in scoolines:
        cid,dum1,dum2,dum3,sharp,rnd,dum4 = l.split()
        sharpdict[cid] = np.float(sharp)
    # Load the opt file
    if os.path.exists(base+'.opt') is False:
        raise Exception(base+'.opt NOT FOUND')
    optlines = read(base+'.opt')
    optlines = optlines.split('\n')
    if optlines[-1]=='': del optlines[-1] # remove blank last line 
    if optlines[-1]=='': del optlines[-1] # remove blank last line 
    ps = 10
    for o in optlines:
        if o.find('PS')>-1: ps=np.float(o.split('=')[1])
        if o.find('HI')>-1: sat=np.float(o.split('=')[1])
    sat -= 3000  # reduce the number a bit to be safe
    rad = np.ceil(ps)
    # Load the fits file
    if os.path.exists(base+'.fits') is False:
        raise Exception(base+'.fits NOT FOUND')
    im, head = fits.getdata(base+'.fits',header=True)
    ny, nx = im.shape
    # Loop over the stars
    badind = []
    soutlines = []
    for i in range(nstars):
        slines1 = slines[i]
        arr = slines1.split()
        id, x, y, mag = arr[0:4]
        x = np.float(x)
        y = np.float(y)
        xlo = np.int(np.max([np.round(x)-rad,0]))
        xhi = np.int(np.min([np.round(x)+rad,nx]))
        ylo = np.int(np.max([np.round(y)-rad,0]))
        yhi = np.int(np.min([np.round(y)+rad,ny]))
        subim = im[ylo:yhi,xlo:xhi]
        nbad = np.sum(subim>=sat)
        # Get sharp value
        if id in sharpdict.keys():
            sharp = sharpdict[id]
        else:
            sharp = 0.0
        # Some saturated pixels or bad sharp value, remove
        if (nbad>0) | (sharp<0.2) | (sharp>1.0):
            badind.append(i)
            if nbad>0:
                print('Star '+str(id)+'  '+str(nbad)+' saturated pixels')
            else:
                print('Star '+str(id)+'  bad sharp value. '+str(sharp))
        # No saturated pixels, keep
        else:
            soutlines.append(slines1)
    # Create output lines
    outlines = lines[0:3].copy()
    if len(soutlines)>0: outlines+=soutlines
    outlines += ['']
    outline = '\n'.join(outlines)
    print(str(len(soutlines))+' out of '+str(nstars)+' stars left')
    # Write to output file
    if os.path.exists(outfile): os.remove(outfile)
    print('Writing to output list '+outfile)
    write(outfile,outline)
