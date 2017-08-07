#!/usr/local/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function

""" 
    Generate MCH, MAG and RAW files from IMG (.als and .head) files
"""
__author__     = "Antonio Dorta"
__copyright__  = "Copyright 2016, The Local Group in Multi-Dimensions | SIEie@IAC"
__credits__    = [""] 
__license__    = ""
__date__       = "2017-03-02"
__version__    = "0.1.3"
__maintainer__ = "Antonio Dorta"
__email__      = "adorta@iac.es"
__status__     = "Developtment"


import os
import glob
import numpy as np
import sys
import time
import copy
from subprocess import Popen, PIPE
from astropy.io.fits import getdata
from math import ceil, log10, sqrt
from numpy.linalg import inv
import random




########################################################################################

class bcolors:
    """Defining some colors for screen message
    """
    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    ERROR     = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'


########################################################################################

def print_stderr(*args, **kwargs):
    # print in stderr
    print(*args, file=sys.stderr, **kwargs)


########################################################################################

def exit_error_msg (error_msg, exit=True):
    """Print an error message and exit

     error_msg -- message to print
    print_help -- if True, help about config file will be displayed
    """

    # Print text with a special format (Error in red)
    print (bcolors.ERROR + "\n\n\tERROR!!! "+ error_msg + bcolors.ENDC + "\n")
    print (sys.exc_info())

    # Print also in STDERR
    print_stderr ("\n\n\tERROR!!! " + error_msg + "\n")
    print_stderr (sys.exc_info())

    # Check if force exit
    if exit:
        sys.exit(-1)


########################################################################################

def print_title (text):
    """Print text with a special format (Title in bold)"""
    print (bcolors.BOLD + text + bcolors.ENDC)


########################################################################################

def print_subtitle (text):
    """Print text with a special format (Subtitle in blue"""
    print (bcolors.OKBLUE + text + bcolors.ENDC)


########################################################################################

def print_warning (text):
    """Print text with a special format (Warning in yellow)"""
    print (bcolors.WARNING + text + bcolors.ENDC)


########################################################################################

def print_info (text):
    """Print text with a special format (Info in green)"""
    print (bcolors.OKGREEN + text + bcolors.ENDC)


########################################################################################

def delfile(fname):
    """Delete a file. Ignore errors (like file does not exist)
    
    fname -- filename (and path) of file to be deleted
    """

    try:
        os.remove(fname)
    except OSError:
        pass


########################################################################################

def get_radius(): 
    """Build the values of the radius needed as parameter in deomaster input
    It is a range [INIT:END:STEP] with float numbers where the last value (END) has to
    appear and right after will be repeated a number of times (REPEAT_END). 
    Values are separted by newlines and last value is a special number to quit (EXIT)

    Ex. # Critical match-up radius: 
          3 2.95 2.90 ... 1.10 1.05 1.00 1 1 ... 1 1   (1 x 15 times)
    """
    INIT         = 3.0       # Initial value of the range
    END          = 1         # Final value of the range (included)
    STEP         = -0.05     # Step
    REPEAT_END   = 15        # Number of extra times that last value will be repeated
    EXIT         = 0         # Special exit value (to specify the end of radius)
    SEP          = "\n"      # Separator between values

    # It can be also done with numpy arange and extra conversion process from float64
    # We will do using simple loops
    radius = ""
    rad_val = INIT

    # Generate main radius
    while (rad_val >= END):
         radius += str(rad_val) + SEP
         rad_val += STEP

    # Add repetitions of last value
    for i in xrange (REPEAT_END):
        radius += str(END) + SEP

    # Add final exit value (NO separator after exit!)
    radius += str(EXIT) 
    return radius


########################################################################################

def f_writeln(f,txt):
    """Simply write a file adding a newline after text.
       We use this function to make code more legible

       f   -- file descriptor (given by open(...)
       txt -- text to write to file
    """
    f.write(txt+"\n")


########################################################################################

def exec_cmd_stdin(cmd, fn_input, log=True):
    """Execute a command in shell with a stdin file. 
       If requested, two log files will be created with stdout and stderr

            cmd -- command to execute
       fn_input -- filename of input file (used as stdin)
            log -- if True, log files for stdout (fn_input.out) and stderr (fn_input.err) will be created
    """

    try:
        with open(fn_input, "r") as f:
            print ("Executing command '" + cmd + " < " + fn_input + "'")
            # Open the input file and create a subprocess to execute the command with it
            p = Popen([cmd], stdout=PIPE, stdin=f, stderr=PIPE)
            stdout_data,stderr_data = p.communicate()

        if log:
           # Logs have been requested...
            try:
                with open (fn_input+".out", "w") as f:
                    # Logging stdout in fn_input.out
                    f.write(stdout_data)
                with open (fn_input+".err", "w") as f:
                    # Logging stderr in fn_input.err
                    f.write(stderr_data)
                print ("LOGs for this command were created in files " + fn_input+".out and " + fn_input+".err")
            except:
                print_warning ("There were problems when creating logs for command '" + cmd + " < " + fn_input + "'")
                print (sys.exc_info())
                    
    except:
        exit_error_msg ("There were errors when calling  " + cmd + " with " + fn_input)


########################################################################################

def get_filename(fname_in, desc, oblig=True):
    """ Get a proper filename given a name that might contain wildcards (* or ?)
        If several files match the given name, then the first of them will be chosen
        and a message will be showed. If there are no files matching the given name
        an error will be showed (execution will be halted if oblig == True)

    fname_in -- input filename (it can include wildcards: * ?
        desc -- description (only to show in warning or error messages)
       oblig -- if True, exit if there are no files matching fname_in 
    """
    # Get MCH file of the directory:
    fname_list = glob.glob(fname_in)
    if len(fname_list) == 0:
        msg = "NO " + desc + " file has been found matching filename " + fname_in
        if oblig:
            exit_error_msg(msg) 
        else:
            print_warning(msg)
            fname_out = fname_in
    if len(fname_list) > 1:
        print_warning("Several " + desc + " files have been found! Using first of them")    
    fname_out = fname_list[0]
    #print_info ("Using " + desc + " file: " + fname_out + " (if not correct, change it in config file)")
    return fname_out;


########################################################################################

def get_mch(fn_mch_in, refine_mch, field, chip):
    """  Use information of filters to generate the MCH file 
         Use external command 'daomaster' to to get final MCH file
         'daomaster' only accept stdin arguments, so will
         generate input files each time and use stdin redirection.
         'daomaster' does NOT accept char '-' in filenames, we will 
         use soft links with char '_' that were created in preprocess.

     fn_mch_in -- MCH file with list of input files
    refine_mch -- Mode (number of coefficients). It should be 0 (no change) or 6, 12, 20
         field -- field identificator (used to generate proper filenames for new files)
          chip --  chip identificator (used to generate proper filenames for new files)

       RETURNS
        fn_mch: filename of output MCH file 
    """

    # Get output MCH filename (delete file, just in case it exists)
    fn_mch = get_mch_addstar_fn(field, chip)
    delfile(fn_mch)   
 
    if refine_mch == 0 or not refine_mch in ALLOWED_MODES:
        # If refine_mch is 0 (or is NOT a valid mode), then MCH file will NOT be changed,
        # just a symlink between old and new filenames will be created
        if refine_mch != 0:
            print_warning("Mode " + str(refine_mch) + " is NOT ALLOWED. MCH file will be NOT refined")

        try:
            os.symlink(fn_mch_in, fn_mch)
        except:
            # If link fails, most probably it already exists 
            # (or error will be found later: full disk, permissions, etc.)
            pass
        # RETURN since we will NOT run daomaster
        return fn_mch

    print_title ("\n\nUSING 'daomaster' TO GENERATE MCH FILE (from .alf files) WITH MODE " + str(refine_mch))
    print_title ("============================================================================")

    # Get filenames for current filter. Delete older files to avoid overwriting question when using daomaster
    fn_input   = field + "-daomaster_" + chip + ".input"

    #  Generate daomaster input file 
    radius = get_radius()
    print  ("Automatically creating input file for daomaster: " + fn_input)
   
    try:
        with open(fn_input, "w") as f:
            print_info ("Next files will be created in this stage: " + fn_mch)

            # daomaster input file
            f_writeln(f, fn_mch_in)          # File with list of input files (MCH)
            f_writeln(f, "2 0 99")           # Minimum number, minimum fraction, enough frames:       
            f_writeln(f, "1")                # Maximum sigma
            f_writeln(f, str(refine_mch))    # Mode: See ALLOWED_MODES
            f_writeln(f, radius)             # Critical match-up radius
            f_writeln(f, "n")                # NO:  Assign new star IDs? 
            f_writeln(f, "n")                # NO:  A file with mean magnitudes and scatter? 
            f_writeln(f, "n")                # NO:  A file with corrected magnitudes and errors?
            f_writeln(f, "n")                # NO: A file with raw magnitudes and errors? 
            f_writeln(f, "y")                # YES: A file with the new transformations?  
            f_writeln(f, fn_mch)             # Transfromations filename (MCH)
            f_writeln(f, "n")                # NO:  A file with the transfer table? 
            f_writeln(f, "n")                # NO:  Individual .COO files? O
            f_writeln(f, "n")                # NO:  Simply transfer star IDs? 

        # Execute command with generated input!!
        exec_cmd_stdin("daomaster", fn_input)

    except:
        exit_error_msg ("There were errors when running daomaster < "+ fn_input)

    return fn_mch


########################################################################################

def get_input_stars(fn_stars, cols_order, stars_shuffle):
    """ Read input stars from file and store it in a dictionary using columns order
    
    fn_stars -- file containing stars magnitudes
    cols_order -- array with the order of the filters
    stars_shuffle -- if True, then stars from input file will be shuffled before using them

    RETURNS:
      * stars: Dictionary with stars ID (@) and magnitudes per filter 
    """

    print_title ("\n\nGETTING STARS FROM " + fn_stars)
    print_title ("============================================================================")

    try:
        # Get output filename and read RAW input file (skip first lines of file header)
        if not stars_shuffle:
            stars_raw = np.loadtxt(fn_stars, unpack=True, ndmin=2) # No shuffling, direct read (Unpack)
        else:
            stars_raw = np.loadtxt(fn_stars, ndmin=2) # We will shuffle stars. np.shuffle() only works in the 
            print("Shuffling " + str(len(stars_raw)) + " stars.")   # first dimension, so we will NOT unpack 
            np.random.shuffle(stars_raw)              # Shuffle rows (each col is a filter)
            stars_raw = np.transpose(stars_raw)       # Transpose (equivalent to Unpack, now each row is a filter)
            print("Shuffling done!")
     
        pos = 0
        stars = {}
        for col in cols_order:
            if col != IGNORE_COL:  # Ignore marked columns
                stars[col] = stars_raw[pos]
                if col == STAR_ID:
                	  print ("Column " + str(pos+1) + ": stars ID. Total stars: " + str(len(stars[col])))
                else:
                	  print ("Column " + str(pos+1) + ": filter " + col + ". Total stars: " + str(len(stars[col])))
            else:
                print ("Column " + str(pos+1) + ": IGNORED")
            pos += 1

    except:
        exit_error_msg ("There were errors when processing data from "+ fn_stars + " to generate stars input file")
	
    return stars
    

########################################################################################

def in_rectangle(corners, point):
    """Check if a point is contained in a rectangle given by two opposite corners

     corners -- array with 4 values that cointains the two opposite corners of
                the rectangle that defines the CCD
                corners = [x0, y0, x1, y1]  ->   x0 <= x1, y0 <= y1
       point -- array with the two coordinates (X,Y) 
                point = [x, y]

                           +--------------o (x1, y1)
                           |              |
                           | o (x, y)     |
                           |              |
                  (x0, y0) o--------------+

    RETURNS:
     * True if point is located inside (or in the borders) of the rectangle. False if not
    """
    #return ((corners[0] <= point[0]) and (point[0] <= corners[2]) and
    #        (corners[1] <= point[1]) and (point[1] <= corners[3]))

    return (corners[0] <= point[0] <= corners[2] and
            corners[1] <= point[1] <= corners[3])


########################################################################################

def get_mode(ncols, fname, main_mode, ALLOWED_MODES):
    """Returns the "mode" of the MCH file (number of coefficients Cn) used in DAOMASTER.
       MCH file format:
      'filename 1' C(1,1) C(1,2) C(1,3) C(1,4) C(1,5) C(1,6) DMAG(1) SMAG(1) C(1,7) C(1,8) ... C(1,20)
      ...
      'filename n' C(n,1) C(n,2) C(n,3) C(n,4) C(n,5) C(n,6) DMAG(n) SMAG(n) C(n,7) C(n,8) ... C(n,20)
      MODE (number of Cn) is given by the number of cols in file MCH, EXCEPT filename, DMAG and SMAG -> ncols - 3

              ncols -- number of columns in MCH file
              fname -- name of the MCH file (to display it in case of error)
          main_mode -- mode of the main MCH file, in order to compare it is the same (default None)
      ALLOWED_MODES -- allowed models in DAOMASTER

    RETURNS:
      * "mode" of the MCH file
    """
    NUM_IGNORE_FIELDS = 3     # Number of fields that are no coefficents in each road: filename, dmag, smag = 3
   
    # Get mode 
    mode = ncols - NUM_IGNORE_FIELDS

    # Check values
    if main_mode and mode != main_mode:
        exit_error_msg ("Mode of MAIN MCH FILE (%s) and MCH FILE %s (%s) are NOT the same" % (fname, main_mode, mode),  False)
    if ALLOWED_MODES and mode not in ALLOWED_MODES:
        exit_error_msg("Mode %s NOT allowed. Allowed modes: %s   Check MCH file (%s)!" % (mode, ALLOWED_MODES, fname), False)

    return mode
 
########################################################################################

def crowdingmultipro(max_iters, field, chip, mode, mch_fnames, mch_data, caja, mtrans, 
                     chip_info, radcent, dimfield, distance, corners, MAX_COEF, numcaj, 
                     numstar, offset, shift):
    """
    Apply the transformations and generate the .add files with results
    Make .add files (one per line in file.mch) to be used with add task in DAOPHOT.
    Developed by adorta starting from a MATLAB code Version 2.7 05-Jul-2010  
     - Code changed to perform several transformation in cascade
     - Code changed to transform also magnitudes
     - Code changed to write also magnitude file (.mag) per mock
     - Code changed to allow the limit of Mocks (iters) that will be generated

       
   max_iters -- if > 0, it will limit the number of iters (mocks generated)
   
       field -- field we are processing (usually FX: F1, F2, etc.)
 
        chip -- chip we are processing (usually YY: 01, 02 ... 62)

        mode -- mode of the main MCH file (used in DAOMASTER. See ALLOWE_MODES)

  mch_fnames -- List of filenames that appears in the MCH file 
                provided by DAOMASTER task used to obtain the photometry.

    mch_data -- Data (coefficients) obtained from MCH file provided by DAOMASTER 
                task used to obtain the photometry.

        caja -- dictionary containing the magnitudes of the stars to be injected:
                Magnitude must be red filter
		            shift in first mag
	              shift to second mag

      mtrans -- Transformation matrix of main MCH file
  
   chip_info -- info about images from chip file needed by magnitude transformations

     radcent -- Distance in pixel between the centroids of the stars
                The centroids of two adjacent stars will be separated this quantity

    dimfield -- Field size covered by all the images in file.mch and offset in X and Y
                to be applied: [xmin ymin xmax ymax]
    
    distance -- Distance Modulus to calculate calibrated magnitude

     corners -- Coordinates of the corners of the images
                [x0, y0, x1, y1]  x0 <= x1, y0 <= y1

                         +--------------o (x1, y1)
                         |              |
                         |              |
                (x0, y0) o--------------+

    MAX_COEF -- Maximum number of coefficients in MCH file

      numcaj -- When the total number of synthetic stars cannot be fitted in one image, 
                this algorithm splits it in several files. 
                Numcaj is the number for the first file. By default numcaj=1
                Ex: numcaj=1 gives file1.add, file2.add, ...
                    numcaj=8 gives file8.add, file9.add, ...
                This can be used when running crowdingmultipro several times, to avoid
	              overwritting the old files and star from the (last+1).

     numstar -- Number to add to the ID star. 0 by default
                Ex: If numstar=100, file1.add contains:
                   101 X1 Y1 mag1
                   102 X2 Y2 mag2
                   ...
	             Similar to NUMCAJ, it must be used to assign correct ID to stars when
               more than one run of crowdingmultipro are required.

     offset -- random values of offset [offx, offy]

      shift -- shift in magnitude to move the stars in the CMD after every iteration

    """

    print_title ("\n\nPROCESSING TRANSFORMATIONS AND CREATING ADD FILES...")
    print_title ("============================================================================")

    # -----------------------------------------
    # SET DEFAULT VALUES (if empty)
    # -----------------------------------------

    if numcaj is None:
        numcaj = 1
    if numstar is None:
        numstar = 0
    if offset is None:
        lim = radcent/2
        #offset= [offx                    , offy                    ]
        offset = [random.uniform(-lim,lim), random.uniform(-lim,lim)]
    if shift is None:
        lim = 0.07
        shift = random.uniform(-lim,lim)

    # Get dimensions
    radcent /= 2
    xmin = dimfield[0]
    ymin = dimfield[1]
    xmax = dimfield[2]
    ymax = dimfield[3]

    offx = offset[0]
    offy = offset[1]
    ncol = xmax
    nrow = ymax

    print ("Random parameters: Mag shift=%s. offx=%s, offy=%s" % (shift, offx, offy))

#-----------------------------------------------------------------------------

# PROCESS MCH FILE

    """
    #DISABLED!!!
    mode = 0
    mch_fnames = []
    mch_data = []
    try:
        with open(mch) as f:
            for raw_line in f:

                line = raw_line.replace("'", "").split()
                # If mode is not set, get mode
                if mode == 0:
                    mode = get_mode (len(line), mch, main_mode, None)

                # First col: filename
                mch_fnames.append(mch.split(".")[0]+"-"+line[0].split(".")[0]+"_")
                # Other cols: vals (convert from string to float
                mch_data.append([float(x) for x in line[1:]])
    except:
        print_warning ("\nWARNING!!! MCH file '%s' does not exist or it has an invalid format. Processing next mch file...\n" % mch)
        return None
    """

    # Get number of images (number of lines in MCH files)
    nimages = len(mch_data)
    xsize   = int(((xmax - xmin + 1) / (radcent * 2)))
    ysize   = int(((ymax - ymin + 1) / (radcent * 2)))
    arin    = ceil(xsize * sqrt(2)) * ysize;
    print("Mode: %s, nimages: %s, xsize: %s, ysize: %s" % (mode, nimages, xsize, ysize))

    # SEPARATION between lines (distribute artificial stars in triangles, NOT squares)
    XSEP = sqrt(0.5)
 
    # Get number of stars
    cn = len(caja[STAR_ID]);
    print ('Total number of artificial stars: ' + str(cn))
    print ('Maximun number of artificial stars in each frame: ' + str(int(arin)))
    cnmax = cn
    if cnmax < arin:
        print_warning ('Do not waste my calculation time, you can still add ' + str(int(arin-cnmax)) + " stars!")
        ########### sys.exit(-1)


    # Get max number of stars per output file
    addmax = int(ceil(float(cnmax)/arin))
    if cnmax > arin:
        print ('Too much artificial stars for a single file')
        print ('It will be used ' + str(addmax) + ' files to complete the total number of artificial stars')
    if max_iters > 0 and max_iters < addmax:
        print_warning("Limiting number of Mocks to " + str(max_iters))
        addmax = max_iters
    max_files = int(round(999/nimages)-1)
    if addmax > max_files:
       exit_error_msg('Too many open files\nMaximun number is: ' + str(max_files) + "\nDecrease the total number of artificial stars", False)


    #-----------------------------------------------------------------------------

    # Now the matrix will be built, depending on MODE
    # Build also the filename array
    mchdat = np.empty((0,MAX_COEF), float)
    filenames = [] 

    # Create "add" files (we will follow MATLAB script and use same variables and names):
    for iadd in xrange(addmax):
        for i in xrange(nimages):
            # FILL with 0.0 till reach MAX_COEF  (current number of elements: mode + 2)
            mchdat = np.append(mchdat, np.array([mch_data[i]]), axis=0)
            fname = get_iter_filename(mch_fnames[i], field, chip, iadd+numcaj, FILE_ITER_EXT)
            filenames.append(fname)
            

    ###########################################
    ###### T R A N S F O R M A T I O N S ######
    ###########################################
    """
    we assign coeficients in order to transform the coordinates:
    'image 1' C(1,1) C(1,2) C(1,3) C(1,4) C(1,5) C(1,6) DMAG(1) SMAG(1) C(1,7) C(1,8) ... C(1,20)
    ...
    'image n' C(n,1) C(n,2) C(n,3) C(n,4) C(n,5) C(n,6) DMAG(n) SMAG(n) C(n,7) C(n,8) ... C(n,20)
   
    to solve the system:
   
    Mode 2
    X =  C(n,1) + C(n,3)*Xn + C(n,5)*Yn
    Y =  C(n,2) + C(n,4)*Xn + C(n,6)*Yn
                                 
    where C(n,3) = C(n,6) = 1 and  C(n,5) = C(n,4) = 0
   
    Mode 4:
    X = C(n,1) + C(n,3)*Xn - Sign(cros)*C(n,4)*Yn
    Y = C(n,2) + C(n,4)*Xn + Sign(cros)*C(n,3)*Yn
                                 
    Mode 6:
    X =  C(n,1) + C(n,3)*Xn + C(n,5)*Yn
    Y =  C(n,2) + C(n,4)*Xn + C(n,6)*Yn
   
    Mode 12
    X =  C(n,1) + C(n,3)*Xn + C(n,5)*Yn + C(n,7)*XDOS + C(n, 9)*XY  + C(n,11)*YDOS
    Y =  C(n,2) + C(n,4)*Xn + C(n,6)*Yn + C(n,8)*XDOS + C(n,10)*XY  + C(n,12)*YDOS
   
    Mode 20
    X = C(n,1) + C(n,3)*Xn + C(n,5)*Yn + C(n,7)*XDOS + C(n, 9)*XY + C(n,11)*YDOS + C(n,13)*XS*XDOS + C(n,15)*YS*XDOS + C(n,17)*XS*YDOS + C(n,19)*YS*YDOS
    Y = C(n,2) + C(n,4)*Xn + C(n,6)*Yn + C(n,8)*XDOS + C(n,10)*XY + C(n,12)*YDOS + C(n,14)*XS*XDOS + C(n,16)*YS*XDOS + C(n,18)*XS*YDOS + C(n,20)*YS*YDOS
   
    where:
    Sign(cros) -> sign of (C(n,3)*C(n,6) - C(n,4)*C(n,5))
    XS = 2*(Xn-1)/(NCOL-1) - 1
    YS = 2*(Yn-1)/(NROW-1) - 1
    XY = XS*YS
    XDOS = 1.5*XS^2 - 0.5
    YDOS = 1.5*YS^2 - 0.5
    DMAG = magnitud media pesada imagen-a-imagen respecto de la primera
    SMAG = varianza de DMAG imagen-a-imagen
    """
    
    # FIRST TRANSFORMATION, data from MAIN MCH file (NOT NEEDED HERE!)
    ###mtrans_all = [main_mtrans, None]
    #mtrans_all = [mtrans]
 
    # SECOND TRANSFORMATION, data from SECONDARY MCH files
    mtrans = np.delete(mchdat,np.s_[6,7],1)
    dmag   = mchdat[:,[6]]


    # Set init values
    numcajorg = numcaj
    nfopen  = 0
    fidimag = []

    # Open output files for writing
    for k in xrange(nimages):
        fidimag.append(open(filenames[k], "w"))

    # Create the MAG file with stars to be added
    fn_mag = get_iter_filename(mch_fnames[0], field, chip, numcaj, ".mag").replace("_", "-add_")
    f_mag = open(fn_mag, "w+")
    # Copy HEADER from ref. image .als (3 first lines)
    mag_header = ""
    with open(mch_fnames[0]+".als", "r") as f_als:
        n_line = 0
        for als_line in f_als:
            mag_header += als_line
            n_line += 1
            if n_line >= 3:
                break
    f_mag.write(mag_header)


    # Init arrays of stars inside/outside rectangle (just to print some stats)
    stars_in  = [0] * nimages
    stars_out = [0] * nimages

    # Apply transformations for ALL STARS
    # ACCESS to the CCD will be performed SEQUENTIALLY
    # (No random accesses like in Matlab code, to speed up
    # the execution since we avoid to compute the visited positions
    xaux = 0 
    yaux = 0 

 
    rowpos = 0
    print ('Writing set of files %s...' % numcaj)
    for star_pos,star_id in enumerate(data_in[STAR_ID]):

        # Get initial positions
        xpos_init = xmin + (2 * (xaux + 1) - 1) * radcent + offx
        ypos_init = ymin + (2 * (yaux + 1) - 1) * radcent + offy

        # Write line in file .mag for the current star. Format:
        # ID  XREF YREF   F1 F1ERR  F2 F2ERR ... FN FNERR  CHI  SHARP  FLAG  PROB
        # ALL ERRORS ARE 0.0000, CHI: 1.0000 SHARP=0.0000 FLAG:0 PROB:1.00
        # Number of decimals and spaces are STRICT!!!
        # Each line can only have 12 cols (first line could have 3 extra cols)
        # Second and consecutive lines has 27 spaces before first data
        MAX_MAG_COLS = 12
        num_cols = 0
        mag_line = "%9d %8.3f %8.3f" % (star_id, xpos_init, ypos_init)
       
        # -------------------------------------------------------------------------

        #############################################
        # TRANSFORM ALL IMAGES
        #############################################
        for nimg in range(nimages):

            # Second transformation, data from current image
            #mtrans_all[0] = mtrans[nimg,:]
            mtrans_all = [mtrans[nimg,:]]

            # Set initial positions
            xpos = xpos_init
            ypos = ypos_init

           
            # BEGIN TRANSFORMATIONS
            # This code has been prepared in order to work with several consecutive transformations,
            # that's why we use lists of transformations instead of the variable
            for mtransn in mtrans_all:

                # TRANSFORMATION
                matref = np.array([[1, mtransn[0], mtransn[1]],  [0, mtransn[2], mtransn[3]],  [0, mtransn[4], mtransn[5]]])
                invmatref = inv(matref)

                xs = 2 * (xpos - 1) / (ncol - 1) - 1
                ys = 2 * (ypos - 1) / (nrow - 1) - 1
                xy = xs * ys

                x2 = 1.5 * (xs ** 2) - 0.5
                y2 = 1.5 * (ys ** 2) - 0.5
            
                Tx = np.asanyarray(- mtransn[6::2].dot([x2, xy, y2, xs*x2, ys*x2, xs*y2, ys*x2])).sum()
                Ty = np.asanyarray(- mtransn[7::2].dot([x2, xy, y2, xs*x2, ys*x2, xs*y2 ,ys*x2])).sum()
                #posstar = [[1, xpos, ypos] + [0, Tx, Ty]] * invmatref
                posstar = np.dot([1, Tx+xpos, Ty+ypos], invmatref)

                # Use same xpos and ypos to concatenate next transformation (if any)
                xpos = posstar[1]
                ypos = posstar[2]

            # TRANSFORMATIONS DONE

            
            # After transformations, SAVE TO FILE (xpos, ypos) ONLY if it is inside rectangle given by corners
            if in_rectangle (corners, [xpos, ypos]):
                # This star is INSIDE the limits, transform the MAGNITUDE and write it to ADD
                try:
                    # Get calibrated color: colsign * (band - colband)
                    chip_img = chip_info[nimg]
                    color = chip_img['COLSIGN'] * (data_in[chip_img['BAND']][star_pos] - data_in[chip_img['COLBAND']][star_pos]) 
                    # PERFORM MAGNITUDE TRANSFORMATION
                    calmag = get_mag(chip_img, data_in[chip_img['FILTER']][star_pos], color, chip_img['FILTER'], distance, magext)
                except:
                    exit_error_msg("There was a problem when processing magnitude for image " + mch_fnames[nimg])

                # Write it to ADD file 
                fidimag[nimg + (numcaj - numcajorg) * nimages].write('%6i %8.3f %8.3f %8.3f\n' % 
                                (star_id, xpos, ypos, calmag))
                mag_line += "%9.4f   0.0000" % calmag
                stars_in[nimg]  += 1
            else:
                # This point is outside corners, just ignore it
                mag_line += "%9.4f   9.9999" % 99.9999
                stars_out[nimg] += 1

            # Check if we have to go to a new line in MAG file (every 12 cols we need to add new line)
            num_cols += 2
            if num_cols % MAX_MAG_COLS == 0: mag_line += "\n%27s" % ''
 

#-----------------------------------------------------------------------------
        # END MAG file with fixed fields: CHI SHARP FLAG PROB

        num_cols += 2
        mag_line += "   1.0000   0.0000"
        if num_cols % MAX_MAG_COLS == 0: mag_line += "\n%23s" % ''
        mag_line += "    0   1.00\n"
        #mag_line += "   1.0000   0.0000    0   1.00\n"
        f_mag.write(mag_line)


#-----------------------------------------------------------------------------

        # UPDATE COUNTERS: Get ready for next iteration over the stars loop 
        # (Go to next star and next col of CCD)
        cn -= 1
        yaux += 1
    		# If we are in the last col, go to first position of next row 
        if (yaux >= ysize):
            xaux += XSEP
            rowpos += 1
            # If we are in even row, y position is half in order to distribute stars in triangles
            if (rowpos % 2) == 0:
                yaux = 0.0
            else:
                yaux = 0.5


        # Check if CCD is FULL (xaux >= xsize) and there are still some stars to process (cn > 0)
        if (xaux >= xsize) and (cn > 0):

            # Check if we have reached the limit of iters (Mocks). If so, quit loop!
            if numcaj >= addmax: break        

            # CCD IS FULL, write in new files
            xaux = 0
            yaux = 0
            numcaj += 1
            nfopen += 1

            # Close all open files
            for f in fidimag:
                if not f.closed:
                    f.close()
            # Open new files
            for k in xrange(nimages):
                fidimag.append(open(filenames[nfopen * nimages + k], 'w'))

            # Close previous MAG file and create new one for current iteration
            f_mag.close()
            fn_mag = get_iter_filename(mch_fnames[0], field, chip, numcaj, ".mag").replace("_", "-add_")
            f_mag = open(fn_mag, "w+")
            f_mag.write(mag_header)
 
            print ('Writing set of files %s...' % numcaj)


#-----------------------------------------------------------------------------

    # ALL OPERATIONS DONE
    # Close all open files
    for f in fidimag:
        if not f.closed:
            f.close()
    f_mag.close()


#-----------------------------------------------------------------------------

    # It could happen that due to the discards one (or more) of the mocks 
    # has 0 stars. Then .add files will be empty and the creation of .fits
    # files will fail. If that happens (any of the mocks is empty), then
    # we will totally REMOVE that mock.

    last_add_list = glob.glob(field+FILE_ITER_SEP+str(addmax)+"-*_"+chip+"."+FILE_ITER_EXT)
    for fn_add in last_add_list:
        f_info = os.stat(fn_add)
        #print(fn_add + " -> " + str(f_info.st_size))
        if f_info.st_size == 0:
            # One of the add files is empty... Delete all files related to this mock!!!
            print_warning(fn_add + " file is empty!! Removing all files of iter #" + str(addmax))
            all_last_mock = glob.glob(field+FILE_ITER_SEP+str(addmax)+"-*")
            for fn_rm in all_last_mock:
                #print("Removing file: " + fn_rm)
                os.remove(fn_rm)

            # Update the max number of mocks
            addmax -= 1
            break
      

    
#-----------------------------------------------------------------------------

    # Show final stats about number of included and discarded stars
    print ("All files written:")
    for i in xrange(nimages):
        print ("  Image %2d (%s). Total stars written (inside rectangle): %s of %s (discarded: %s)" % 
                  (i+1, filenames[i], stars_in[i], stars_out[i] + stars_in[i], stars_out[i]))

    return addmax

#-----------------------------------------------------------------------------


########################################################################################

def print_config_file_help():
    """Print help about the config file syntax. Returns None
    """

    print_title("You must specify a CONFIG FILE with the next format:")
    print ("\tLine 1: main MCH file (wildcards allowed: * ?)")
    print ("\tLine 2: File with data from chips (wildcards allowed: * ?)")
    print ("\tLine 3: Input file with magnitudes (wildcards allowed: * ?)")
    print ("\tLine 4: Order of input file columns (separated by space, it MUST contain stars ID marked with " 
                + STAR_ID + ". Ignore cols using: " + IGNORE_COL + ")")
    print ("\tLine 5: Max CCD size [format]: xmax ymax")
    print ("\tLine 6: radcent (set a float value OR use * to calculate radcent using .opt files)")
    print ("\tLine 7: Field size covered by all the images (dimfield). You can specify a MCH file to calculate values OR set them: xmin xmax ymin ymax")
    print ("\tLine 8: Distance Modulus used to calculate calibrated magnitude (it will be added to absorption)")
    print ("\t# Empty lines or those beginning with '#' will be ignored\n")
    print_title("\nEXAMPLE OF CONFIG FILE:\n")
    print ("\t# Main MCH file\n\tF*_*_*.alf.mch\n")
    print ("\t# Chips data\n\tField??_chips.fits\n")
    print ("\t# Input file\n\tinput_stars*.txt\n")
    print ("\t# order of input file cols (" + STAR_ID + " for ID and " + IGNORE_COL +" to ignore)\n\t" + STAR_ID + " u g r i z " +  IGNORE_COL + "\n")
    print ("\t# Max CCD size: xmax ymax\n\t#          +------o (xmax,ymax)\n\t#          |      |\n\t#    (0,0) o------+\n\t2048 4096\n")
    print ("\t# radcent: use * to autocalculate OR set the value\n\t*\n\t# 13.5\n")
    print ("\t# dimfield: file OR xmin ymin xmax ymax\n\tFMchN.mch\n\t# -8 -4 1040 516\n")
    print ("\t# distance: 0\n")


########################################################################################

def clean_comments(f):
    """Read config file and get data from it, cleaning comments (lines beginning with "#") 
       and empty lines, after stripping.
       Returns line content of next non-empty line (removing comments). 
       Raise an exception if file has no more data

    f -- file descriptor of config file
    """

    for line in f:
        # Strip: remove all blanks 
        data = line.strip()
        if data == "" or data.startswith("#"):
           # Skip empty lines or those beginning with '#' for comments
           continue
        return data

    # No data, launch exception!
    raise

########################################################################################

def get_radcent_from_psf():
    """ THIS FUNCTION IS NOT USED ANYMORE SINCE radcent IS CALCULATED 
        FROM PS AND FI INFO IN OPT FILES!!
        Returns value of radcent calculated from PSF files in the CWD
        Data from PSF files are obtained always from the first two
        columns of the second line of each PSF file. There is ONE blank
        space before first col and then each col has a fixed size of 13 chars.
        
      

        PSF Format:
           Line 1:...Header...
           Line 2:_11111111111112222222222222...nnnnnnnnnnnnn
           Line 3:...
    """
    FILE_EXT     = "*.psf"  # Extension of files to process
    PROCESS_LINE = 2       # Number of the line where data is located (first line is 1)
    IGNORE_CHARS = 1       # Number of chars to ignore from the beginning of line
    FIELD_LENGTH = 13      # Number of chars of each value
    NUM_FIELDS   = 2       # Number of consecutive values to process
    FIRST_FIELD  = 0       # Position of the first value (first position is 0)

    sum_values = 0.0
    num_values = 0
    num_files  = 0

    print_title ("\n\nPROCESSING PSF FILES TO CALCULATE VALUE OF RADCENT:")
    print_title ("============================================================================")

    fname_list = glob.glob(FILE_EXT)
    for fname in fname_list:
         with open(fname) as f:
            num_line = 0
            num_files += 1
            file_ok = False
            for line in f:                         # We read the current file
                num_line += 1
                if (num_line == PROCESS_LINE):     # We are in the correct line to get data
                    if len(line) >= (NUM_FIELDS + FIRST_FIELD) * FIELD_LENGTH:   # Check if it is long enough!
                        for i in xrange (FIRST_FIELD, NUM_FIELDS+FIRST_FIELD):
                            try:
                                # Get all values, convert them to FLOAT and add it to previous values
                                sum_values += float(line[IGNORE_CHARS+FIELD_LENGTH*i : IGNORE_CHARS+FIELD_LENGTH*(i+1)])
                                num_values += 1
                                file_ok = True
                            except:
                                # There was an error reading the current value, go to next file 
                                file_ok = False
                                break
                    # We have processed the data line, break loop and go to next file
                    break
              
            if file_ok:
               print ("File %s: Processed..." % fname)
            else:
               print_warning("File %s: THERE WERE ERRORS AND SOME DATA WAS IGNORED!!! Check file format and content!!" % fname)

    if num_values == 0:
        exit_error_msg("There are no PSF files in current directory or they have a wrong format") 

    radcent = sum_values/num_values * 5 * 2.35  # 5 * ( 2.35 * sigma) = 5 * FWHM   [assuming a gaussian]
    print_info ("%s files and %s values were processed.\nAverage: %s, radcent = %s" % (num_files, num_values, sum_values/num_values, radcent))


    return radcent

########################################################################################

def get_radcent_from_opt(field, chip):
    """ Returns value of radcent calculated from OPT data stored in .dat files in the CWD
     radcent is the distance in pixel between the centroids of the stars
     The centroids of two adjacent stars will be separated this quantity

        PSF Format:
           Line 1:...Header...
           Line 2:_11111111111112222222222222...nnnnnnnnnnnnn
           Line 3:...

       field -- field we are processing (usually FX: F1, F2, etc.)
        chip -- chip we are processing (usually YY: 01, 02 ... 62)

    """

    print_title ("\n\nCALCULATE VALUE OF RADCENT USING FI AND PS INFO OF OPT FILES:")
    print_title ("============================================================================")


    try:
        FI_file = field + "-chipFI_" + chip + ".dat"
        PS_file = field + "-chipPS_" + chip + ".dat"
        max_FI = max(np.loadtxt(FI_file, unpack=True))
        max_PS = max(np.loadtxt(PS_file, unpack=True))
        radcent = max_PS + 2 * max_FI + 1
        print_info("MAX PS: " + str(max_PS) + ". MAX FI: " + str(max_FI) + ". Radcent (PS+2FI+1): " + str(radcent))
        return radcent


    except:
        exit_error_msg("There are no files " + FI_file + " and/or " + PS_file + " for current chip or they have a wrong format") 

########################################################################################

def get_dimfield_from_file(fname, corners, radcent):
    """Return dimfield ([xmin,ymin, xmax,ymax]) obtained from a MCH file. dimfield is calculated 
       adding the min/max values of X and Y in MCH File to the respective values in corners.
       We will use ceil(min) and floor(max) in order to convert float to int
    
         fname -- filename of MCH file
       corners -- original dims [x0,y0, x1,y1]
       radcent -- distance in pixel between the centroids of the stars. The centroids of two adjacent stars will be separated this quantity

       RETURNS:
         dimfield
    """ 

    COL_X = 1     # Position of column X in MCH file (first col is 0)
    COL_Y = 2     # Position of column Y in MCH file (first col is 0)
    xmin = xmax = ymin = ymax = None

    try:
        with open(fname) as f:
            for raw_line in f:

                # Process line, remove filename delimiter(') and split in columns
                line = raw_line.replace("'", "").split()

                if xmin is None or float(line[COL_X]) < xmin:
                    xmin = float(line[COL_X])
                if xmax is None or float(line[COL_X]) > xmax:
                    xmax = float(line[COL_X])
                if ymin is None or float(line[COL_Y]) < ymin:
                    ymin = float(line[COL_Y])
                if ymax is None or float(line[COL_Y]) > ymax:
                    ymax = float(line[COL_Y])
       
        # Reduce size using radcent/2 and convert to int with floor(max-radcent/2) and ceil(min+radcent/2)  
        xmin = int(ceil(xmin + radcent/2))
        ymin = int(ceil(ymin + radcent/2))
        xmax = int(xmax - radcent/2)
        ymax = int(ymax - radcent/2)

        # After processing the whole file, we add the floor(max)/ceil(min) values,
        # but ONLY if those values are lower for min or bigger for max than those in corners
        dimfield = [min(corners[0], corners[0] + xmin), 
                    min(corners[1], corners[1] + ymin),
                    max(corners[2], corners[2] + xmax),
                    max(corners[3], corners[3] + ymax)]

        # Display values
        print_title ("\n\nTO CALCULATE DIMFIELD USING FILE %s" % fname)
        print_title ("============================================================================")
        print ("Original corners: %s" % corners) 
        print ("Values obtained from %s: %s (after reducing all limits with radcent/2: %.3f) " % (fname, [xmin, ymin, xmax, ymax], radcent/2))
        print_info ("Final dimfield: %s" % dimfield)

        return dimfield

    except:
        exit_error_msg ("File '%s' used to calculate dimfield does not exist or it has an invalid format...\n" % fname)
        return None






########################################################################################

def process_mch(mch_file, input_data, cols_order, star_ini, number_stars, ALLOWED_MODES, MAX_COEF):
    """
    Read Main MCH file and the input file, in order to establish the correspondences between them, since we
    need to know which column of input file corresponds to each secondary MCH file of those that appear in the
    main MCH file.

         mch_file -- filename of the main MCH file (with .alf images)
       input_data -- data in input file 
       cols_order -- information about columns of input file given by user in the config file
         star_ini -- first star to process (count begins in 0)
     number_stars -- number of starts to process (defined by user)
    ALLOWED_MODES -- List of allowed models of MCH files 
         MAX_COEF -- Maximum number of coeffiecients in MCH files (with NO filename field)
    
    RETURNS:
                mode: Mode of the CH file
          mch_fnames: array of filenames of all files listed in the MCH file
            mch_data: matrix with all coefficients of MCH file 
          input_data: matrix with all data from columns of input file 
              nstars: number of stars to process (after checking that number_star is not over the limit)
    """
    #--------------------------------------
    # LOAD AND PROCESS MAIN MCH FILE
    #--------------------------------------
    mode = 0
    mch_data = np.empty((0,MAX_COEF), float)
    mch_fnames = []

    try:
        with open(mch_file) as f:
            for raw_line in f:
                line = raw_line.replace("'", "").split()
                # If mode is not set
                if (mode == 0): 
                   mode = get_mode(len(line), mch_file, None, ALLOWED_MODES)

                # First col: filename
                mch_fnames.append(line[0].replace(".alf", ""))
   
                # Other cols: Cn values (convert from string to float) and add 0.0 till reach MAX_COEF
                row = [float(x) for x in line[1:]] + ([0.0] * (MAX_COEF - (mode + 2)))
                mch_data = np.append(mch_data, np.array([row]), axis=0)
    except:
        exit_error_msg("Main MCH file '%s' does not exist or it has an invalid format" % mch_file, False)
    
    #-----------------------------------------------------------------
    # LOAD INPUT FILE AND FIND ASSOCIATIONS BETWEEN COLS AND MCH FILES
    #-----------------------------------------------------------------

    if number_stars == 0:
        nstars = len(input_data) 
    elif len(input_data) >= number_stars:
        nstars = number_stars
    else:
        nstars = len(input_data)
        print_warning ("WARNING!!! Parameters start_ini=%s and number_stars=%s are out of range. Correcting to number_stars=%s -> stars: [%s, %s]" % (star_ini+1, number_stars, nstars, star_ini+1,star_ini+nstars))
      

    return [mode, mch_fnames, mch_data, input_data, nstars]

########################################################################################

def process_daophot(field, chip, numiters):
    """ 
       FITS files will be created from .add files using daophot.

       field -- field we are processing (usually FX: F1, F2, etc.)
        chip -- chip we are processing (usually YY: 01, 02 ... 62)
    numiters -- Number of iterations (mocks)
    """
    IMG_EXT  = ".fits"
    PSF_EXT  = ".psf"
    DAO_EXT  = ".opt"
    DAO_GAIN = "GA"
    DAO_SEP  = "="

    print_title ("\n\nUSING 'daophot' TO GENERATE FITS FILES (from .add files)")
    print_title ("============================================================================")

    # Get all images names (.fits files):
    images = glob.glob(field+"-*_"+chip+IMG_EXT)
    if len(images) == 0:
        exit_error_msg ("There are no " + IMG_EXT + " images!!")
    print_info("Processing " + str(len(images)) + " " + IMG_EXT + " images")

    for img in images:
        # Get basename removing extension
        fn_img = img.replace(IMG_EXT, "")

        # Get GAIN value from .opt
        gain = ""
        try:
            with open(fn_img+DAO_EXT, "r") as f:
                for line in f:
                    if line.startswith(DAO_GAIN):
                        gain = line.split(DAO_SEP)[1].strip()
                        break

        except:
            exit_error_msg ("There were errors when getting GAIN value from file " + fn_img + DAO_EXT)

        if gain == "":
            exit_error_msg ("There were errors when getting GAIN value from file: " + fn_img + DAO_EXT)

        # Get filenames for current filter. Delete older files to avoid overwriting question when using daophot
        fn_input = fn_img + "-daophot.input"

        # Get base image filename

        #  Generate daophot input file 
        print ("Automatically creating input file for daophot: " + fn_input)
   
        try:
            with open(fn_input, "w") as f:
                # daophot input file
                f_writeln(f, "1")                   # Fake value   
                f_writeln(f, "1")                   # Fake value
                f_writeln(f, "OPTIONS")             #    
                f_writeln(f, fn_img+DAO_EXT)        #       
                f_writeln(f, "")                    #    
                f_writeln(f, "ATTACH " + fn_img)    #    
                f_writeln(f, "ADDSTAR")             #    
                f_writeln(f, fn_img+PSF_EXT)        #        
                f_writeln(f, "1")                   #    
                f_writeln(f, gain)                  #    

                for itr in xrange(1,numiters+1):
                    f_writeln(f, get_iter_filename(fn_img, field, chip, itr, "."+FILE_ITER_EXT))    #        
                    f_writeln(f, "")                       #    

                f_writeln(f, "EXIT")                   #    

            # Execute command with generated input!!
            exec_cmd_stdin("daophot", fn_input, True)

        except:
            exit_error_msg ("There were errors when running daophot < "+ fn_input)

    return None

########################################################################################

def process_argv(args, mch_ext):
    """
    Check and process command-line arguments

      args -- array of arguments (sys.argv)
   mch_ext -- Extension of MCH files

        RETURNS: 
        * main_mch_file: filename of the main MCH file
        * chips_file:    filename of the chips data file
        * input_file:    filename of the input file 
        * cols_order:    array with associations between Main MCH files and Input columns 
        * corners:       array with opposite corners of CCD: [x0,y0, x1,y1]
        * radcent:       value of radcent
        * dimfield:      array: [xmin, xmax, ymin, ymax]
        * distance:      value of Distance Modulus used to calculate the calibrated magnitudes
    """
    argc           = 1
    #ARG_FIELD      = argc; argc+=1
    #ARG_CHIP       = argc; argc+=1
    ARG_MCHFILE    = argc; argc+=1
    ARG_CHIPSFILE  = argc; argc+=1
    ARG_MAXMOCKS   = argc; argc+=1
    ARG_STARSCOLS  = argc; argc+=1
    ARG_STARSSHUF  = argc; argc+=1
    ARG_MAGEXT     = argc; argc+=1
    ARG_MAXCCDSIZE = argc; argc+=1
    ARG_DIMFIELD   = argc; argc+=1
    ARG_RADCENT    = argc; argc+=1
    ARG_DISTANCE   = argc; argc+=1
    ARG_REFINEMCH  = argc; argc+=1
    NUM_ARGS       = argc

    error_msg = ""

    # ---------------------------
    # CHECK INPUT VALUES
    # ---------------------------

    if len(args) != NUM_ARGS:
        # Wrong number of parameters
        exit_error_msg("Syntax: %s <1:mch_file> <2:chips_file> <3:maxmocks> <4:starscols> <5:starsshufle> " 
                       "<6:magext> <7:maxccdsize> <8:dimfield> <9:radcent> <10:distance> <11:refine_mch>" % args[0], True)


    try:
        # MCH filename has format FX-NNNNNN_YY.mch) -> Get FIELD (FX) and CHIP (YY) from it!
        mch_fn = os.path.basename(args[ARG_MCHFILE])         # Remove path, get only FILENAME
        field  = mch_fn.split("-")[0].strip()                # Field (FX) goes BEFORE "-"
        chip   = mch_fn.split("_")[1].split(".")[0].strip()  # Chip  (YY) goes BETWEEN "_" and "." 
        print_info ("Processing Field " + field + " and Chip " + chip)
    except:
        error_msg = " * Field and chip data\n"


    # Process config file and read each field
    # ---------------------------------------

    try:
        # Main MCH file: string
        # Filename format is FX_NNNNNN_YY.alf.mch -> replace FX-NNNNNN_YY.mch to get new format
        mch_args   = mch_fn.replace("-", "_").replace(".mch", ".alf.mch")
        main_mch_file = get_filename(mch_args, "Main MCH", True)
    except:
        error_msg += " * Main MCH file\n"
 
    try:
        # chips file: string
        chips_file = get_filename(args[ARG_CHIPSFILE], "Chips data", True)
    except:
        #print ("Unexpected error:", str(sys.exc_info()))
        error_msg += " * Chips data file\n"

    try:
        # max_iters: int (limit in the number of mocks)
        max_iters = int(args[ARG_MAXMOCKS])
    except:
        error_msg += " * Max Number of Mocks\n"

    try:
        # starsShuffle: bool
        stars_shuffle = (args[ARG_STARSSHUF] != '0') and (args[ARG_STARSSHUF] != '')
    except:
        error_msg += " * Stars Shuffle\n"
 
    
    try:
        # Input file: string
        stars_file = field + "-stars_" + chip + ".txt"
        input_file    = get_filename(stars_file, "Input stars", True)
    except:
        error_msg += " * Input file\n"

    try:
        # Column order: array of string
        cols_order    = args[ARG_STARSCOLS].replace(" ","").split(',')
        if not STAR_ID in cols_order:
            error_msg += " * Column order (it MUST contain " + STAR_ID + " to specify stars ID)\n"
    except:
        error_msg += " * Column order\n"

    try:
        magext = {}
        # MAGEXT values could be directly given or a file could be specified

        ##########################################################
        # MAGEXT -> VALUES
        ##########################################################
        # if magext begins with "values:", it should be a common float value or array of floats (one per filter)
        if args[ARG_MAGEXT].lower().startswith("values:"):
            num_filters = 0
            mag_values = args[ARG_MAGEXT].split(":")[1]
            # Try to get the values splitting the list
            magext_aux  = [float(x) for x in mag_values.replace(" ","").split(',')]

            # Go for each filter:
            for flt in cols_order:

                # Check whether we have already processed all filters
                if len(magext_aux) != 1 and num_filters >= len(magext_aux):
                    num_filters += 1
                    break
 
                # Skip STAR_ID or columns to ignore 
                if flt != STAR_ID and flt != IGNORE_COL:
                    if len(magext_aux) == 1:
                        magext[flt] = magext_aux[0]             # One common value
                    else:
                        magext[flt] = magext_aux[num_filters]   # List of values
                    num_filters += 1

            if len(magext) != num_filters or (len(magext_aux) != 1 and len(magext_aux) != num_filters):
                error_msg += " * MagExt should be one common value or has a value for each filter\n"
        else:
            ##########################################################
            # MAGEXT -> FILE
            ##########################################################
            # Magext should be a FILE
            MAGEXT_FILE_DEFAULT = "extinction"
            magext_file = ""
            # Check if specified file exists
            if args[ARG_MAGEXT] and os.path.isfile(args[ARG_MAGEXT]):
                magext_file = args[ARG_MAGEXT]
            elif os.path.isfile(MAGEXT_FILE_DEFAULT):
                magext_file = MAGEXT_FILE_DEFAULT
            else: 
                # There is no file, use 0 instead 
                for flt in cols_order:
                    if flt != STAR_ID and flt != IGNORE_COL:
                        magext[flt] = 0             # One common value
                #error_msg += " * MagExt is not a valid file or an array of values\n"
         
            # Process file 
            if magext_file != "": 
                with open(magext_file) as f:
                    for line in f:
                        # Ignore comments 
                        if not line.strip().startswith("#"):
                            # Split line: left filter, right value
                            line_data = line.strip().split()
                            # Check if filter is in list and assign value! 
                            if len(line_data) == 2 and line_data[0] in cols_order:
                                magext[line_data[0]] = float(line_data[1])

                # Check we have data for all filters
                for flt in cols_order:
                    if flt != STAR_ID and flt != IGNORE_COL and not flt in magext:
                        error_msg += " * File " + args[ARG_MAGEXT] + " does NOT cointain Mag.Ext. value for filter '" + flt + "'\n"
                        break

    except:
        error_msg += " * MagExt\n"

    try:
        # max_ccd_size: array of 2 floats: xmax ymax
        max_ccd_size  = [float(x) for x in args[ARG_MAXCCDSIZE].replace(" ","").split(',')]
        if len(max_ccd_size) == 2:
            corners = [0.0, 0.0, max_ccd_size[0], max_ccd_size[1]]
        else:
            error_msg += " * Max CCD size must have 2 values: xmax,ymax\n"
    except:
        error_msg += " * Max CCD size\n"

    try:
        # Radcent: if user specifies "*", this value is calculated. 
        # If not, get float value in config file
        if args[ARG_RADCENT].strip() == "*":
            # If read data is "*", calculate radcent value from PSF files
            radcent = get_radcent_from_opt(field, chip)
        else:
            # Get float value from config file
            radcent = float(args[ARG_RADCENT])
            if radcent <= 0:
                error_msg += " * Radcent should be > 0 (" + str(radcent) + " was set)\n"
    except:
        error_msg += " * Radcent\n"

    try:
        # dimfield: filename or array of 4 floats
        value = args[ARG_DIMFIELD].split(',')
        if len(value) == 1 and value[0] == '*':
            # Only one value, it must be the filename to calculate dimfield
            dimfield = get_dimfield_from_file(main_mch_file, corners, radcent)
        elif len(value) == 4:
            dimfield = [float(x) for x in value]
        else:
            error_msg += " * dimfield must have 4 values: xmin,ymin,xmax,ymax\n"
    except:
        error_msg += " * dimfield\n"

    try:
        # distance: float value
        distance = float(args[ARG_DISTANCE])
    except:
        error_msg += " * Distance modulus\n"

    try:
        # refine_mch: int value
        refine_mch = int(args[ARG_REFINEMCH])
    except:
        error_msg += " * Refine MCH\n"



    # Chech whether there were some errors
    if error_msg:
        # Error, the config file is not valid
        # (some data is missing or has no proper format)
        exit_error_msg("Format of config file is not correct\nThe next fields are missing or have incorrect format:\n" + error_msg, True)
        return None

    else:
        # No errors: print and return values 
        print_title ("\n\nPROCESSING INPUT DATA:")   
        print_title ("============================================================================")
        print_subtitle("    Main MCH file: " + main_mch_file)       
        print_subtitle("  Chips data file: " + chips_file)       
        print_subtitle("       Input file: " + input_file)
        print_subtitle("     Column order: " + str(cols_order))
        print_subtitle("          Corners: " + str(corners))
        print_subtitle("          Radcent: " + str(radcent))
        print_subtitle("         Dimfield: " + str(dimfield))  
        print_subtitle("            Field: " + field)  
        print_subtitle("             Chip: " + chip) 
        print_subtitle("  Mag. Extinction: " + str(magext))
        print_subtitle(" Distance Modulus: " + str(distance))
        print_subtitle("       Refine MCH: " + str(refine_mch))
        if max_iters > 0: print_subtitle("      Limit Mocks: " + str(max_iters)) 
        if stars_shuffle: print_subtitle("    Shuffle stars: YES") 
        else: print_subtitle("    Shuffle stars: NO") 


        return [main_mch_file, chips_file, input_file, stars_shuffle, cols_order, magext, corners, radcent, 
                dimfield, distance, field, chip, max_iters, refine_mch]

########################################################################################

def get_mag(data, calmag, colsub, star_filt, distance, magext):
    """ Perform magnitude transformation following next equation:
        color = colsub  * 3.07 * EBV
        instmag = calmag + zpterm + amterm*X + colterm*COLOR + apcorr - 2.5*alog10(exptime)  [X = airmass]

       data   -- chip values for the current image
       calmag -- magnitude from star star[chip[FILTER] (it should be = star[chip[BAND]]
       colsub -- color given by: chip[COLSIGN] * (star[chip[BAND]] - star[chip[COLBAND]])
    star_filt -- filter used in this image
     distance -- Distance Modulus used to calculate the calibrated magnitude
       magext -- magnitude extinctions

    RETURNS:
     * calculated magnitude according to the given equation
    """

    absorption = magext[star_filt] * data['EBV']
    calmag += absorption + distance
    instmag = calmag + data['ZPTERM'] + data['AMTERM'] * data['AIRMASS'] + data['COLTERM'] * colsub + data['APCOR'] - 2.5 * log10(data['EXPTIME'])
    return instmag
 

########################################################################################

def get_chip_info(fn_chips, fnames_img):
    """ Read chip file to extract the info for current images

    fn_chips -- filename of the chip info file
  fnames_img -- list of current images filenames

    RETURNS:
          data: matrix with data for current images
    """

    # Get Data from FITS
    data = []
    try:
        data_raw = getdata(fn_chips, 1)
    except:
        exit_error_msg("File " + fn_chips + " not found or it has no proper format")
        
    # Save time creating a structure where chip info for each image could be directly accessed
    # using image index (nimg), instead of searching each image every time we process a star
    # Profiling shows that doing that in this way could considerably reduce the execution time
    for img in fnames_img:
        # Filter using base name (filename without extension)
        row = data_raw[data_raw['BASE'] == img]

        # Check for errors (no row of data found) or warnings (several rows were found)
        if len(row) == 0:
            exit_error_msg("File " + fn_chips + " does not have info for image " + img)
        elif len(row) > 1:
            print_warning("There are several different data for image " + img + " in chips info file " + fn_chips + ". Using the first of them")

        # Append first search result
        data.append(row[0])
            
    return data



########################################################################################

def get_iter_filename(orig_fn, field, chip, iteration, new_ext=None, sep=None):
    """ 
    Get the filename INCLUDING the Mock number (iteration) using orig_fn as basename
    For instance, FX-001234_YY.eee will be FXM1-001234_YY.eee for the first Mock (iter=1)
    It will also change extension if a new extension (new_ext) is given

     orig_fn -- Original filename
       field -- field we are processing (usually FX: F1, F2, etc.)
        chip -- chip we are processing (usually YY: 01, 02 ... 62)
   iteration -- current iteration (mock number)
     nex_ext -- new extension (if present, old extension will be replaced by the new one)
        sep  -- separator between field and mock number (iteration)

     RETURNS:
        new_fn: New filename after applying the modifications
    """

    if sep is None:
        sep = FILE_ITER_SEP

    new_fn = orig_fn.replace(field+"-", field + sep + str(iteration) + "-", 1)
 

    if not new_ext is None:
        # Check if extension should be also changed
        if new_ext != "" and not "." in new_ext:
            new_ext = "." + new_ext
        new_fn = os.path.splitext(new_fn)[0] + new_ext
 
    return new_fn

########################################################################################

def get_mch_addstar_fn(field, chip):
    """ 
    Get the filename of the MCH file according to field and chip

       field -- field we are processing (usually FX: F1, F2, etc.) 
        chip -- chip we are processing (usually YY: 01, 02 ... 62)

    RETURNS:
       New filename of MCH file according to current field and chip
    """

    return field + "-" + MCH_ADDSTAR + "_" + chip + ".mch"

########################################################################################

def create_mch(main_mch, field, chip, numiters, dest_mch):
    """
    Create the new MCH files (replacing .alf with .als) and adding the Mock 
    number (iteration). Since we are introducing new filenames, we need some 
    related files (.psf, .mag, etc.)
    Those files are the same as the original ones, we only need to copy them, 
    but it's better to just make a symbolic link (symlink). 
    For example, when creating FXMN-001234_YY.alf, we also need the FXMN-001234_YY.psf, 
    but this file is the same for all iterations, so we need to do: 
    XMi-001234_YY.psf -> FX-001234_YY.psf for each iteration from 1..N
    We need to that for ALL extension in SYMLINKS_EXT.

    main_mch -- Filename of the main MCH file
       field -- field we are processing (usually FX: F1, F2, etc.)
        chip -- chip we are processing (usually YY: 01, 02 ... 62)
    numiters -- Number of iterations (number of Mocks)
    dest_mch -- Original name of the MCH file that will be trasnformed 
    """

    OLD_EXT = ".alf"
    NEW_EXT = ".als"
    SYMLINKS_EXT = [".als", ".raw", ".opt", ".psf", ".als.opt", ".ap" , ".mag", ".log",
                    ".weights", ".scale", ".zero", "_comb.psf", "_comb.opt", "_comb.als.opt", "_shift.mch"]

    # Correct DESTINATION MCH FILENAME in case it is wrong: FX_YYYYYY_NN.alf.mch -> FX-YYYYYY_NN.mch
    dest_mch = dest_mch.replace(OLD_EXT, "").replace(field+"_", field+"-")


    print_title ("\n\nCreating new MCH files from file " + main_mch + "...")
    print_title ("============================================================================")
    # Get all images names:

    try:
        mch_f = [None] * numiters
        for i in xrange(numiters):
            # Build all new MCH files (one per Mock)
            dest_fn = get_iter_filename(dest_mch, field, chip, i+1)
            mch_f[i] = open(dest_fn, "w")
            print("Creating MCH file: " + dest_fn)
        
        # Process old MCH file (original one, we will use it like a template)
        with open(get_mch_addstar_fn(field, chip), "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    for i in xrange(numiters):
                        # Get OLD filename (Original one from MCH. 
                        # Then remove everything to get ONLY the basename (remove extension, etc.)
                        old_fn = line.split(" ")[0].replace("'","").replace(OLD_EXT,"")
                        # Build the new names that includes different extension and MOCK number
                        new_fn = get_iter_filename(old_fn, field, chip, i+1) 
                        # Replace old names with new ones
                        xline = line.replace(old_fn+OLD_EXT, new_fn+NEW_EXT)
                        f_writeln(mch_f[i], xline)
                          
                        # Create SYMLINKS
                        # We need symlinks from new files with Mock number to the old original files
                        # for the given extension (SYMLINKS_EXT)
                        for ext in SYMLINKS_EXT:
                            if os.path.isfile(old_fn+ext):
                                try:
                                    os.symlink(old_fn+ext, new_fn+ext)
                                except:
                                    continue
        # Close all open files
        for i in xrange(numiters):
            mch_f[i].close() 
 
    except:
        exit_error_msg ("Error when creating MCH files!")

########################################################################################

def update_apcor(field, chip, mch_fnames, numiters):
    """
    After creating the new filenames for each mock, we have to update the Aperture Correction (ApCor)
    Info about ApCor is located in file with name apcor.lst that has been transfer to the main
    directory and also to each field/chip if they were separated
    We will create a new partial ApCor file just with the info of the current images, we need to
    replicate the info of each image:
    Original file:
    FX-......-NNa.del    0.0000000
    New file:
    FX-......-NNa.del    0.0000000
    FXMY-......-NNa.del    0.0000000

       field -- field we are processing (usually FX: F1, F2, etc.) 
        chip -- chip we are processing (usually YY: 01, 02 ... 62)
  mch_fnames -- List of filenames of all MCH files
    numiters -- Number of iterations (number of Mocks)
    """

    # Filenames (original one and partial one)
    APCOR_FILENAME = "apcor.lst"
    FNAME2RM = "a.del"       # In ApCor filenames have extra text to be removed         

    print_title ("\n\nUPDATING ApCor INFORMATION FOR THIS CHIP...")
    print_title ("============================================================================")


    partial_filename = field+"-partial_apcor_"+chip+".lst"
    # Use the original names of the images (FXX-......-NN)
    fnames_apcor = copy.copy(mch_fnames)
    try:
        with open(APCOR_FILENAME, "r") as f_old:
            with open(partial_filename, "w+") as f_new:
                for line in f_old:
                    if not line.strip().startswith('#'):   # Ignore comments
                        line_data = line.strip().split()   # Separate filename (first col) and values (second col)
                        fname = line_data[0].replace(FNAME2RM, "")  # Remove extra chars

                        # Check if current line has 2 columns and this filename is in our list
                        if len(line_data) == 2 and fname in fnames_apcor:  
                            # Valid image!! Write the original data
                            f_new.write(line)
                            # Add mocks with same value (1..N)
                            for iadd in xrange(1,numiters+1):
                                fname_new = get_iter_filename(fname, field, chip, iadd, FNAME2RM)
                                f_new.write(fname_new + "   " + line_data[1]+"\n")
                            # Remove current image from the pending list
                            fnames_apcor.remove(fname)

                        if len(fnames_apcor) == 0:
                            break   # If there are no more images to be processed, break!!

        # List of pending images should be empty... If not, print warning!!
        if len(fnames_apcor) != 0:
            print_warning("Not all images were added to apcor.lst!! Missing images: " + str(fnames_apcor))
        print_info("Update done! Partial ApCor data for field " + field + " and chip " + chip + " is stored in file " + partial_filename)

    except:
        print_warning("There was an error when processing apcor.lst file. " 
                     +"Aperture Correction info may not be properly updated, CALIB stage could fail!")
        print (sys.exc_info())


 
###################################################################################

def create_trans_eq_calib(field, chip, numiters, chip_info):
    """
    Create the Transforamtion Equations that will be used in CALIB stage
    All data are extracted from Chips info File. Output has 3 lines per image and mock:

    1) First line: filename,  band name,      color name, zero-point, airmass, color,   airmass*color, color^2
       Fields:     BASE,      FILTER or BAND  note1       ZPTERM,     AMTERM   COLTERM, AMCOLTERM,     COLSQTERM
    2) Second line:  errors (note2)
    3) Third line: empty

    note1: COLSIGN indicates how to construct the color (COLSIGN=1 means color=BAND-COLBAND. COLSIGN=-1 means color=COLBAND-BAND)
    note2: XXXXSIG are the error terms that go on the second line (example: ZPTERMSIG, AMTERMSING, etc.)
    Example:
    F5-00517150_43  g  g-r  -0.4089    0.1713   -0.1193   0.0000   0.0000
                             0.0040   -0.0000    0.0001   0.0000   0.0000
       field -- field we are processing (usually FX: F1, F2, etc.) 
        chip -- chip we are processing (usually YY: 01, 02 ... 62)
    numiters -- Number of iterations (number of Mocks)
   chip_info -- Info about each image (filters, coefficients, etc.)


    """

    print_title ("\n\nCREATING TRANSFORMATION EQUATIONS FOR THIS CHIP:")   
    print_title ("============================================================================")

    try:
        line_format = "%18s  %s  %s  %8.4f %8.4f %8.4f %8.4f %8.4f\n"
        partial_filename = field+"-partial_calib_"+chip+".trans"
        with open(partial_filename, "w+") as f:
            for info in chip_info:
                for iadd in xrange(1, numiters+1):
                    # Add Mock number to base name
                    fname = info['BASE'].replace("-", "M"+str(iadd)+"-")
                   
                    # Get color name based on COLSIGN (see note1)
                    if info['COLSIGN'] == 1: color = info['BAND'] + "-" + info['COLBAND']
                    else:                    color = info['COLBAND'] + "-" + info['BAND']

                    line  = line_format % (fname, info['FILTER'], color,
                            info['ZPTERM'],    info['AMTERM'],    info['COLTERM'],    info['AMCOLTERM'],    info['COLSQTERM'])

                    line += line_format % ('', ' ', '   ',
                            info['ZPTERMSIG'], info['AMTERMSIG'], info['COLTERMSIG'], info['AMCOLTERMSIG'], info['COLSQTERMSIG'])
                    
                    line += "\n"
                    f.write(line)

        print_info("Done! Partial Transf. Eq. for field " + field + " and chip " + chip + " are stored in file " + partial_filename)

    except:
        print_warning("There was an error when generating Transformation Equations for this chip. " 
                     +"Info may not be properly created and CALIB stage could fail!")
        print (sys.exc_info())


    

###################################################################################
###################################################################################
############################                        ###############################
############################    M    A    I    N    ###############################
############################                        ###############################
###################################################################################
###################################################################################

###################################################################################
# GLOBAL CONFIG
###################################################################################

DEF_CCD_SIZE  = [0, 0, 4096, 2048]      # Default values of CCD size
#ALLOWED_MODES = [2, 4, 6, 12, 20]  # Allowed modes in DAOMASTER
ALLOWED_MODES = [6, 12, 20]  # Allowed modes in DAOMASTER (Transformations of mode 2 and 4 NOT implemented!!)
MAX_COEF      = 22                      # Max. number of coefficients in MCH files (with NO filename)
STAR_ID       = "@"
IGNORE_COL    = "*"
FILE_ITER_SEP = 'M'			           # Separator used when creating the "add" files
FILE_ITER_EXT = 'add'              # Add files extension
MCH_FILE_EXT  = '.alf.mch'         # Extension of MCH files
MCH_ADDSTAR   = 'addstar'


# python3 support:
if sys.version_info.major >= 3:
    xrange = range

###################################################################################
# BEGIN PROCESS (Read and proccess input params)
###################################################################################


init_time = time.time()
print ("INIT TIME " + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

# Check command-line arguments and get PARAMS
[orig_mch_file,   \
 chips_file,      \
 stars_file,      \
 stars_shuffle,   \
 cols_order,      \
 magext,          \
 corners,         \
 radcent,         \
 dimfield,        \
 distance,        \
 field,           \
 chip,            \
 max_iters,       \
 refine_mch]      \
 = process_argv(sys.argv, MCH_FILE_EXT)


# Get and process mch_file of MCH 
mch_file = get_mch(orig_mch_file, refine_mch, field, chip)

# Get and store stars magnitude input file
stars = get_input_stars(stars_file, cols_order, stars_shuffle)

star_init = number_stars = 0
# Get data from MAIN MCH and INPUT files
[mode,            \
 mch_fnames,      \
 mch_data,        \
 data_in,         \
 nstars]          \
 = process_mch(mch_file, stars, cols_order, star_init, number_stars, ALLOWED_MODES, MAX_COEF)

# Get info from chips 
chip_info = get_chip_info(chips_file, mch_fnames)


################################################
## Get ready to apply transformations
################################################

# Remove cols #6 (DMAG) and #7 (SMAG) in MAIN TRANSFORM MATRIX
mmtrans  = np.delete(mch_data,np.s_[6,7],1)
#MAX_OFFSET = 10
MAX_OFFSET = (radcent/2)-1
#offset = [offx,                                   offy]
offset  = [random.uniform(-MAX_OFFSET,MAX_OFFSET), random.uniform(-MAX_OFFSET,MAX_OFFSET)]
shift   = 0

###########################################################
# PROCESS TRANSFORMATIONS
###########################################################
numiters = crowdingmultipro(max_iters, field, chip, mode, mch_fnames, mch_data, data_in, mmtrans, 
                            chip_info, radcent, dimfield, distance, corners, MAX_COEF, 
                            None, star_init, offset, shift)

###########################################################
# DAOPHOT: Process each Iteration
###########################################################

process_daophot(field, chip, numiters)
create_mch(mch_file, field, chip, numiters, orig_mch_file)

##########################################################
# AUXILIAR OPERATIONS: Update apcor (APerture CORrection) 
# and create Transformation Equation for CALIB stage
# (partial data only for this chip, later all partial
# data will be gathered and joined in one single file)
##########################################################
update_apcor(field, chip, mch_fnames, numiters)
create_trans_eq_calib(field, chip, numiters, chip_info)

##########################################################
print ("\nEND TIME " + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
print ("Final elapsed time (in seconds): " + str(time.time() - init_time))


