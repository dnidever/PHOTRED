#!/usr/bin/env python
#
# DAOPHOT.PY - SExtractor and DAOPHOT routines
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20180823'  # yyyymmdd

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, Column
import time
import shutil
import subprocess
import logging
#from scipy.signal import convolve2d
from scipy.ndimage.filters import convolve
import astropy.stats
import struct
import tempfile
from dlnpyutils.utils import *

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run DAOPHOT PSF photometry on FITS image.")

    parser.add_argument('--file', '-f', action="store", help="The FITS file to process", default=None)
    parser.add_argument('--verbose', '-v', action="store_true", help="Print out the data as it is gathered.", default=False)
    parser.add_argument('--clobber', '-c', action="store_true", help="Overwrite the output file if it already exists.", default=False)

    args = parser.parse_args()
    file = args.file
    
    print "Running DAOPHOT PSF photometry on ", file

    #   Check that we have all of the files that we need

    #   Copy everything to a temporary directory

    #   temporarily uncompressing fits.fz file

    # 1st step : select PSF candidates and obtain initial PS

    #   run find and phot

    #   Get magnitude limit from log file

    #   run srcfilter.pro on cmn.lst file if it exists

    #   run PICKPSF

    #   increase magnitude limit if nPSF<12

    #   pick psf stars again

    #   run lstfilter on the first list of psf stars

    #   run PSF to construct first PSF

    # 2nd step : filter out bad PSF stars and construct a second version of PSF

    #   while loop running PSF and then removing stars that have saturated, defective or * or ? including running goodpsf

    # 3rd step : subtract neighbor stars and construct a final PSF

    #   subtract neighboring stars with GROUP, NSTAR and SUBSTAR

    #   run PSF again to create the final PSF

    #   while loop for rejecting stars with either '*' or '?' next to them, and run goodpsf.pro

    # 4th step : perform PSF photometry with ALLSTAR

    #   run allstar

    # Step 5 : Prepare for aperture correction

    #   run phot and allstar only for PSF stars and using image with neighbors subtracted

    #   Delete temporarily uncompressed fpack file

    #   fpack compress s.fits file

    #   copy files back to original directory

