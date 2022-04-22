#!/usr/bin/env python

import os
import time
import numpy as np
import subprocess
import tempfile
from astropy.io import fits
from dlnpyutils import utils as dln
from glob import glob
 
def check(irafdir=None,silent=False):
    """
    This program checks that you can run IRAF from IDL. 
 
    Parameters
    ----------
    irafdir : str, optional
       The path to the IRAF directory 
         where login.cl should be located.  If 
         this is not entered then it assumes that 
         ~/iraf/ is the IRAF directory. 
    silent : boolean, optional
       Don't print anything to the screen. 
 
    Returns
    -------
    test : boolean
       True if the test was successful, False if it was not, and -1 if there was an error. 
    out : list
       The IRAF output from the test 
 
    Example
    -------
    test = check_iraf(irafdir=irafdir)
 
    By D.Nidever  Nov.2007 
    Translated to Python by D. Nidever,  April 2022
    """

    test = False  # bad to start with 
     
    # Input IRAF directory 
    diriraf = None
    if irafdir is not None:
        diriraf = os.path.expanduser(irafdir)
        if len(glob(diriraf))==0:
            if silent==False:
                print('DIRECTORY '+str(irafdir)+' DOES NOT EXIST.  TRYING ~/iraf/ INSTEAD.')
            diriraf = None
            
    # No IRAF directory yet 
    if diriraf is None:
        diriraf = os.path.expanduser('~/iraf/')
        if len(glob(diriraf))==0:
            raise ValueError('NO IRAF DIRECTORY. RETURNING')
        diriraf = dln.first_el(diriraf)
        dirtest = os.path.exists(diriraf) 
        if dirtest == False:
            raise ValueError('NO IRAF DIRECTORY. RETURNING')
         
    curdir = os.getcwd()
                           
    # Write a test IRAF script
    cmd = []
    cmd += ['print("")']   # adorta: FIRST LINE WILL BE IGNORED!! 
    cmd += ['print("hello world")']
    cmd += ['logout']
    tid,cmdfile = tempfile.mkstemp(suffix='.cl',prefix="temp")
    dln.writelines(cmdfile,cmd)
     
    # Goto the IRAF directory 
    os.chdir(diriraf)
                      
    # Running IRAF 
    out = subprocess.check_output('cl < '+cmdfile,shell=True)
    if type(out) is bytes:
        out = out.decode()
    out = out.split('\n')  # split into lines

    # Return to original directory 
    os.chdir(curdir)
                      
    # Erasing the temporary files
    if os.path.exists(cmdfile): os.remove(cmdfile)
     
    # See if it printed the message
    gd = dln.grep(out,'hello world',index=True)
    if len(gd) > 0:  # YES, it worked 
        test = True 
     
    # Explain how to fix the login.cl file 
    if test == False and silent==False:
        print('Running IRAF from Python failed!')
        print('Please edit your login.cl file so that it:es not')
        print('print anything to the screen on IRAF startup.')
        print('The most likely cause are the 9 lines after')
        print('"Set the terminal type".')
        print('On PLEIONE the 4 lines after') 
        print('"# Delete any old MTIO lock (magtape position) files."')
        print('can be problematic.')


    return test
              

def run(scriptname,irafdir=None,silent=False):
    """
    This runs an IRAF script. 
    
    Parameters
    ----------
    scriptname : str
       The absolute filename of the IRAF script 
         It must end with "logout". 
    irafdir : str, optional
       The IRAF home directory. 
    silent : boolean, optional
       Don't print anything to the screen 
 
    Returns
    -------
    The IRAF script will be executed and the output printed 
    the screen (unless /silent is set). 
    out : list
        The IRAF output 
 
    Example
    -------

    out = run(scriptname,irafdir)
 
    By D.Nidever  August 2007 
    Similar to a perl script written by Armin Rest 
    Translated to Python by D. Nidever,  April 2022
    """

    # Important directories 
    if irafdir is None:
        irafdir = os.path.expanduser('~/iraf/')
    irafdir = glob(irafdir)
    if len(irafdir) == 0:
        raise ValueError('NO IRAF DIRECTORY')
    irafdir = dln.first_el(irafdir)
    curdir = os.getcwd()
              
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if check(irafdir)==False:
        raise ValueError('IRAF TEST FAILED')
     
    # Go to IRAF home directory
    os.chdir(irafdir)
     
    # Message to the screen 
    if silent==False:
        print('Running IRAF...')
        print('with IRAF home directory: ',irafdir)
        print('Executing script: ',scriptname)
     
    # Execute the script 
    out = subprocess.check_output('cl < '+scriptname,shell=True)
    if type(out) is bytes:
        out = out.decode()
    out = out.split('\n')  # split into lines

    # Print the output to the screen
    if silent==False:
        for l in out: print(l)
     
    # Go back to original directory
    os.chdir(curdir)
                      
 
def imalign(input,reference,coords,output,shifts='',boxsize=7,
            bigbox=11,negative='no',background='INDEF',
            lower='INDEF',upper='INDEF',niterate=3,tolerance=0,
            maxshift='INDEF',shiftimages='yes',interp_type='spline3',
            boundary_type='constant',constant=0.0,trimimages='yes',
            verbose='yes',irafdir=None):
    """
    This runs IRAF's IMALIGN that aligns and shifts images. 
 
    Parameters
    ----------
    input 
        The  input  images  to  be shifted and trimmed.  The input image 
        list should contain the reference image so that its borders  are 
        used in the computation of the overlap region. 
 
    reference 
        The reference image to which the input images will be aligned. 
 
    coords 
        A  text  file  containing the reference image coordinates of the 
        registration objects to be centered in each  image,  one  object 
        per  line  with  the  x and y coordinates in columns one and two 
        respectively. 
 
    output 
        The output images. 
 
    shifts = "" 
        A text file containing the initial estimate for  each  image  of 
        the  shift  in each axis relative to the reference image.  These 
        estimates are used to modify the coordinates of the registration 
        objects  prior  to  centering.   The  format  of the file is one 
        image per line with the x and y shifts in columns  one  and  two 
        respectively.    The   sense   of   the  shifts  is  such  that: 
        Xshift=Xref-Xin and  Yshift=Yref-Yin.   If  shifts  is  null,  a 
        coarse  centering  pass will be made to attempt to determine the 
        initial shifts. 
 
    boxsize = 7 
        The size in pixels of the box to use for  the  final  centering, 
        during  which  all  the sources in coords are recentered in each 
        image using the initial estimate of the relative shift for  each 
        image.   Care should be taken to choose an appropriate value for 
        this parameter, since it is highly data dependent. 
 
    bigbox = 11 
        The size in pixels of the box to use for coarse centering.   The 
        coarse  pass  through  the  centering algorithm is made with the 
        box centered at the nominal position of the first source in  the 
        coordinate  list.   Coarse  centering  is  performed only if the 
        shifts file is undefined.  Care should be  taken  to  choose  an 
        appropriate  value  for  this parameter, since it is highly data 
        dependent.  Large values  should  be  suspect  until  the  final 
        results  are  checked to see that the centering did not converge 
        on the wrong coordinates,  although  the  usual  result  for  an 
        inappropriate  bigbox  size  is  that  the  algorithm  fails  to 
        converge and the task aborts. 
 
    negative = no 
        Are the features negative ? 
 
    background = INDEF 
        The  absolute  reference  level  for   the   marginal   centroid 
        calculation.   If  background  is INDEF, this is set to the mean 
        value (between the thresholds) of the individual sources. 
 
    lower = INDEF 
        The lower threshold for the data.  Individual pixels  less  than 
        this value will be given zero weight in the centroids. 
 
    upper = INDEF 
        The  upper  threshold  for  the data.  Individual pixels greater 
        than this value will be given zero weight in the centroids. 
 
    niterate = 3 
        The maximum number of  centering  iterations  to  perform.   The 
        centering  will  halt  when  this  limit  is reached or when the 
        desired Itolerance is achieved. 
 
    tolerance = 0 
        The tolerance for convergence of the centering algorithm.   This 
        is  the  integral  shift of the centering box from one iteration 
        to the next. 
 
    maxshift = INDEFR 
        The maximum permitted difference  between  the  predicted  shift 
        and  the the computed shift for each object. Objects with shifts 
        greater than maxshift are ignored. If maxshift is  undefined  no 
        shift checking is done. 
 
    shiftimages = yes 
        If  shiftimages  is  yes, the IMSHIFT task will be used to align 
        the images.  If shiftimages  is  no,  the  images  will  not  be 
        aligned, but the coordinates will still be centered. 
 
    interp_type = "spline3" 
        The interpolation function used by the IMSHIFT task. 
 
    boundary_type = "constant" 
        The boundary extension type used by the IMSHIFT task. 
 
    constant = 0. 
        The  constant  used  by  the  IMSHIFT  task  if boundary_type is 
        "constant". 
 
    trimimages = yes 
        If trimimages is yes, the  output  images  will  be  trimmed  to 
        include  only  the region over which they all overlap.  The trim 
        section that is actually used  may  differ  slightly  from  that 
        reported   by   IMCENTROID,  due  to  a  correction  applied  to 
        compensate for the boundary extension "contamination"  near  the 
        edges of the images. 
 
    verbose = yes 
        Print the centers, shifts, and trim section? 
 
    irafdir 
        The IRAF home diretory. 
 
 
    Returns
    -------
    Shifted images with the names in the "output" file. 
    xoff : numpy array
       The offset in X between the original and shifted images 
         xorig = xshift + xoff 
    yoff : numpy array
       The offset in Y between the original ahd shifted images 
         yorig = yshift + yoff 
    trans : numpy array
       The X/Y shifts in the format [nfiles,X/Y]. 
    trimsection : numpy array
       The trim section in the format [Xstart,Xend,Ystart,Yend] 
 
    Example
    -------
  
    out = imalign('@inlist','ref.fits','coordfile','@outlist')
 
 
 
                                    I R A F 
                     Image Reduction and Analysis Facility 
    PACKAGE = immatch 
       TASK = imalign 
    
    input   =         @fits_8.list  Input images 
    referenc=        obj034_8.fits  Reference image 
    coords  =          chip8.coord  Reference coordinates file 
    output  =      @outfits_8.list  Output images 
    (shifts =          chip8.shift) Initial shifts file 
    (boxsize=                    7) Size of the small centering box 
    (bigbox =                   11) Size of the big centering box 
    (negativ=                   no) Are the features negative ? 
    (backgro=                INDEF) Reference background level 
    (lower  =                INDEF) Lower threshold for data 
    (upper  =                INDEF) Upper threshold for data 
    (niterat=                    5) Maximum number of iterations 
    (toleran=                    0) Tolerance for convergence 
    (maxshif=                INDEF) Maximum acceptable pixel shift 
    (shiftim=                  yes) Shift the images ? 
    (interp_=               linear) Interpolant 
    (boundar=              nearest) Boundary type 
    (constan=                   0.) Constant for constant boundary extension 
    (trimima=                  yes) Trim the shifted images ? 
    (verbose=                  yes) Print the centers, shifts, and trim section ? 
    (list   =                     ) 
    (mode   =                   ql) 
 
 
    By D. Nidever    February 2008 
    Translated to Python by D. Nidever, April 2022
    """

     
    # Important directories 
    if irafdir is None:
        irafdir = os.path.expanduser('~/iraf/')
    irafdir = glob(irafdir)
    if len(irafdir)==0:
        raise ValueError('NO IRAF DIRECTORY')
    irafdir = dln.first_el(irafdir)
    curdir = os.getcwd()
     
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #---------------------------------------------------------------
    if check(irafdir=irafdir)==False:
        raise ValueError('IRAF TEST FAILED')
     
    # Input strings 
    sinput = str(input).strip()
    sreference = str(reference).strip()
    scoords = str(coords).strip()
    soutput = str(output).strip()
    sshifts = str(shifts).strip()
    sboxsize = str(boxsize).strip()
    sbigbox = str(bigbox).strip()
    snegative = str(negative).strip()
    sbackground = str(background).strip()
    slower = str(lower).strip()
    supper = str(upper).strip()
    sniterate = str(niterate).strip()
    stolerance = str(tolerance).strip()
    smaxshift = str(maxshift).strip()
    sshiftimages = str(shiftimages).strip()
    sinterp_type = str(interp_type).strip()
    sboundary_type = str(boundary_type).strip()
    sconstant = str(constant).strip()
    strimimages = str(trimimages).strip()
    sverbose = str(verbose).strip()
     
    # Write IRAF script
    cmd = []
    cmd += ['print("")']  # first line will be ignored 
    cmd += ['cd '+curdir ]
    cmd += ['immatch']
    cmd += ['imalign.boxsize = '+sboxsize]
    cmd += ['imalign.bigbox = '+sbigbox]
    cmd += ['imalign.negative = '+snegative]
    cmd += ['imalign.background = '+sbackground]
    cmd += ['imalign.lower = '+slower]
    cmd += ['imalign.upper = '+supper]
    cmd += ['imalign.niterate = '+sniterate]
    cmd += ['imalign.tolerance = '+stolerance] 
    cmd += ['imalign.maxshift = '+smaxshift]
    cmd += ['imalign.shiftimages = '+sshiftimages]
    cmd += ['imalign.interp_type = "'+sinterp_type+'"']
    cmd += ['imalign.boundary_type = "'+sboundary_type+'"']
    cmd += ['imalign.constant = '+sconstant]
    cmd += ['imalign.trimimages = '+strimimages]
    cmd += ['imalign.verbose = '+sverbose]
    cmd += ['imalign("'+sinput+'","'+sreference+'","'+scoords+'","'+soutput+'",shifts="'+sshifts+'")']
    cmd += ['logout']
    tid,cmdfile = tempfile.mkstemp(prefix="temp")  # absolute filename
    dln.writelines(cmdfile,cmd)
     
    # Goto the IRAF directory
    curdir = os.getcwd()
    os.chdir(irafdir)
     
    # Running IRAF 
    out = subprocess.check_output('cl < '+cmdfile,shell=True)
    if type(out) is bytes:
        out = out.decode()
    out = out.split('\n')   # split into lines

    # IMALIGN can fail if it doesn't find any matching sources in 
    # the images. 
     
    # The output
    lo = dln.grep(out,'Shifts',index=True)  # where to start output
    if len(lo)==0:
        lo = 0
    else:
        lo = lo[0]
    for l in out[lo:]: print(l)
     
    # Return to original directory 
    os.chdir(curdir)
     
    # Erasing the temporary files
    if os.path.exists(cmdfile): os.remove(cmdfile)
     
    # What were the shifts 
    # this is not very robust yet.
    infiles = dln.loadinput(sinput)
    ninfiles = len(infiles)
    trans = np.zeros((ninfiles,2),float)
    out_trans = out[lo+1:lo+1+ninfiles] 
    transarr = out_trans.split()
    names = transarr[0,:]
    xshift = np.array(transarr[1,:]).astype(float)
    yshift = np.array(transarr[3,:]).astype(float)
    trans[:,0] = xshift 
    trans[:,1] = yshift 
     
     
    # Getting trim section 
    # [Xstart,Xend,Ystart,Yend]
    gdtrim = dln.grep(out,'#Trim_Section',index=True)
    #gdtrim , = np.where(stregex(out,'#Trim_Section',/fold_case,/boolean) eq 1,ngtrim) 
    tarr = out[gdtrim[0]].split()
    tsection = tarr[2]  # i.e. [4:2046,1:2041] 
    trimsection = [4] 
    tlen = len(tsection) 
    tsection = tsection[0:tlen]
    tsection = tsection[1:tlen-1] 
    tsecarr = tsection.split(',')
    trimsection = tsecarr.split(':')
    trimsection = np.array(trimsection).astype(float)
    #trimsection = (trimsection)(*) 
     
     
    # Because of the trimming this can cause shifts between the 
    # original and shifted/combined images 
    head1 = fits.getheader(sreference) 
    crpix1a = head1['CRPIX1']
    crpix2a = head1['CRPIX2']
    # Image has WCS 
    if str(crpix1a,2) != '0': 
        base = os.path.splitext(os.path.basename(sreference))[0]
        head2 = fits.getheader(base+'.shft.fits') 
        crpix1b = head2['CRPIX1']
        crpix2b = head2['CRPIX2']
        xoff = crpix1a-crpix1b 
        yoff = crpix2a-crpix2b 
        # No WCS, Use LTV1/LTV2 
    else: 
        base = os.path.splitext(os.path.basename(sreference))[0]
        head = fits.getheader(base+'.shft.fits') 
        # Take negative 
        ltv1 = head['LTV1']
        ltv2 = head['LTV2']
        # Take negative 
        # xorig = xshift + xoff 
        # yorig = yshift + yoff 
        xoff = float(-ltv1) 
        yoff = float(-ltv2) 
         
     
    print('Xoff='+str(xoff)) 
    print('Yoff='+str(yoff))

    return xoff,yoff,trans,trimimages

 
def imcombine(input,output,headers='',bpmasks='',rejmask='',
              nrejmasks='',expmasks='',sigma='',logfile='STDOUT',
              combine='average',reject='none',project='no',outtype='real',
              outlimits='',offsets='none',masktype='none',
              maskvalue=0.0,blank=0.0,scale='non',zero='none',weight='none',
              statsec='',expname='',lthreshold='INDEF',hthreshold='INDEF',
              nlow=1,nhigh=1,nkeep=1,mclip='yes',lsigma=3.0,hsigma=3.0,
              rdnoise=0.0,gain=1.0,snoise=0.0,sigscale=0.1,pclip=-0.5,
              grow=0.0,irafdir=None):
    """
    This runs IRAF's IMCOMBINE that combines images 
 
    Parameters
    ----------
    input 
        List of input images to combine.  If the  project  parameter  is 
        "no"  then  all  input  images must have the same dimensionality 
        though they may be of different  sizes.   Otherwise  each  input 
        image   is  handled  separately  and  they  may  have  different 
        dimensionalities. 
 
 
    When the project parameter is "no" all the input images are combined 
    into  a  single  output file.  In this case the following parameters 
    specify only a single file name.   Otherwise  each  input  image  is 
    combined  by  projecting (combining across) the highest dimension to 
    produce a lower dimensional  image.   For  this  type  of  combining 
    there  is  one output for each input and so the following parameters 
    specify matching lists. 
 
 
    output 
        Output combined image(s).  If there are  fewer  than  100  input 
        images  the  names  of  the  input images are recorded in header 
        keywords IMCMBnnn. 
 
    headers = "" (optional) 
        Optional output multiextension  FITS  file(s).   The  extensions 
        are dataless headers from each input image. 
 
    bpmasks = "" (optional) 
        Optional  output bad pixel mask(s) with good values of 0 and bad 
        values of 1.  Output pixels are marked  as  bad  when  no  input 
        pixels  contributed  to the output pixel.  The file name is also 
        added to the output image header under the keyword BPM. 
 
    rejmask = "" (optional) 
        Optional output mask file(s) identifying  rejected  or  excluded 
        pixels.   The  pixel  mask  is  the size of the output image but 
        there is one extra dimension with length equal to the number  of 
        input  images.   Each element of the highest dimension is a mask 
        corresponding to an input image with values of  1  for  rejected 
        or  excluded  pixels and values of 0 for pixels which were used. 
        The order of the masks is the order  of  the  input  images  and 
        image  header  keywords,  indexed by the pixel coordinate of the 
        highest dimension identify the  input  images.   Note  that  the 
        pixel positions are in the output pixel coordinate system. 
 
    nrejmasks = "" (optional) 
        Optional  output pixel mask(s) giving the number of input pixels 
        rejected or excluded from the input images. 
 
    expmasks = "" (optional) 
        Optional output exposure mask(s) giving the sum of the  exposure 
        values   of   the   input  images  with  non-zero  weights  that 
        contributed  to  that  pixel.   Since  masks  are  integer,  the 
        exposure  values  may  be  scaled  to preserve dynamic range and 
        fractional significance.  The scaling values are  given  in  the 
        header  under  the  keywords  MASKSCAL  and  MASKZERO.  Exposure 
        values are computed from the mask values  by  scale  *  value  + 
        zero  where  scale is the value of the MASKSCAL keyword and zero 
        is the value of the MASKZERO keyword. 
 
    sigma = "" (optional) 
        Optional output sigma  image(s).   The  sigma  is  the  standard 
        deviation,  corrected  for  a  finite  population,  of the input 
        pixel  values  (excluding  rejected  pixels)  about  the  output 
        combined pixel values. 
 
 
    logfile = "STDOUT" (optional) 
        Optional  output  log file.  If no file is specified then no log 
        information is produced.  The special filename  "STDOUT"  prints 
        log information to the terminal. 
 
 
    combine = "average" (average|median|sum) 
        Type  of  combining  operation  performed  on  the  final set of 
        pixels   (after   offsetting,   masking,    thresholding,    and 
        rejection).   The  choices  are  "average",  "median", or "sum". 
        The median uses the average of the two central values  when  the 
        number  of  pixels  is even.  For the average and sum, the pixel 
        values are multiplied by the weights (1 if no weighting is used) 
        and  summed.   The average is computed by dividing by the sum of 
        the weights.  If the  sum  of  the  weights  is  zero  then  the 
        unweighted average is used. 
 
    reject = "none" (none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip) 
        Type  of  rejection  operation performed on the pixels remaining 
        after offsetting, masking and thresholding.  The algorithms  are 
        described  in  the  DESCRIPTION  section.  The rejection choices 
        are: 
 
              none - No rejection 
            minmax - Reject the nlow and nhigh pixels 
           ccdclip - Reject pixels using CCD noise parameters 
          crreject - Reject only positive pixels using CCD noise parameters 
           sigclip - Reject pixels using a sigma clipping algorithm 
         avsigclip - Reject pixels using an averaged sigma clipping algorithm 
             pclip - Reject pixels using sigma based on percentiles 
 
    project = no 
        Project (combine) across the  highest  dimension  of  the  input 
        images?   If  "no"  then all  the input images are combined to a 
        single output  image.   If  "yes"  then  the  highest  dimension 
        elements  of  each  input  image are combined to an output image 
        and optional pixel list and sigma images.  Each element  of  the 
        highest dimension may have a separate offset. 
 
    outtype = "real" (none|short|ushort|integer|long|real|double) 
        Output  image pixel datatype.  The pixel datatypes are "double", 
        "real", "long", "integer", unsigned short "ushort", and  "short" 
        with  highest precedence first.  If "none" is specified then the 
        highest precedence datatype of the input images is  used.   When 
        there  is  a  mixture  of  short  and  unsigned short images the 
        highest  precedence  become  integer.   The  datatypes  may   be 
        abbreviated to a single character. 
 
    outlimits = "" 
        Output  region limits specified as pairs of whitespace separated 
        values.  The first two numbers are the limits  along  the  first 
        output  image  dimension,  the  next  two numbers are the limits 
        along the second dimension, and so on.  If the higher  dimension 
        limits  are  not  specified  they  default  to  the  full range. 
        Therefore, if no limits are specified then the  full  output  is 
        created.   Note  that  the  output size is computed from all the 
        input images including offsets if specified and the  coordinates 
        are relative to that size. 
 
 
    offsets = "none" (none|wcs|world|physical|grid|<filename>) 
        Integer offsets to add to each image axes.  The options are: 
 
        "none" 
            No offsets are applied. 
 
        "wcs" or "world" 
            The  world  coordinate  system (wcs) in the image is used to 
            derive  the  offsets.   The  nearest  integer  offset   that 
            matches  the  world  coordinate  at  the center of the first 
            input image is used. 
 
        "physical" 
            The physical coordinate system defined by the  IRAF  LTM/LTV 
            keywords define the offsets. 
 
        "grid" 
            A  uniform  grid  of offsets is specified by a string of the 
            form 
 
                    grid [n1] [s1] [n2] [s2] ... 
 
            where ni is the number of images in dimension i  and  si  is 
            the  step  in  dimension  i.  For example "grid 5 100 5 100" 
            specifies a 5x5 grid with origins offset by 100 pixels. 
 
        <filename> 
            The offsets are given  in  the  specified  file.   The  file 
            consists  of  one  line  per  image with the offsets in each 
            dimension forming the columns. 
 
    masktype = "none" (none|goodvalue|badvalue|goodbits|badbits|!<keyword>) 
        Type of pixel masking to use.  If "none" or  ""  then  no  pixel 
        masking  is  done  even  if  an  image  has an associated  pixel 
        mask.  The other choices are to select the value  in  the  pixel 
        mask  to be treated as good (goodvalue) or bad (badvalue) or the 
        bits (specified as a value) to be treated as good (goodbits)  or 
        bad  (badbits).   In  these cases the pixel mask file name comes 
        from the image header keyword BPM.  If the  parameter  is  given 
        as  "!<keyword>"  where  <keyword> is a header keyword, the mask 
        file comes from the value of that keyword  and  the  mask  value 
        interpretation  is the same as "goodvalue".  Note, if the number 
        of input images becomes too large (currently about 4090 .imh  or 
        2045  projection.   This means all the images must have the same 
        size and dimensionality. 
 
    maskvalue = 0 
        Mask value used with the masktype parameter.  If the  mask  type 
        selects  good  or bad bits the value may be specified using IRAF 
        notation for decimal, octal, or hexidecimal; i.e  12,  14b,  0cx 
        to select bits 3 and 4. 
 
    blank = 0. 
        Output value to be used when there are no pixels. 
 
    scale = "none" (none|mode|median|mean|exposure|@<file>|!<keyword>) 
        Multiplicative  image  scaling  to  be applied.  The choices are 
        none, multiply by the reciprocal of the mode,  median,  or  mean 
        of  the specified statistics section, multiply by the reciprocal 
        of the exposure time  in  the  image  header,  multiply  by  the 
        values  in  a  specified  file, or multiply by a specified image 
        header keyword.  When specified in a file  the  scales  must  be 
        one per line in the order of the input images. 
 
    zero = "none" (none|mode|median|mean|@<file>|!<keyword>) 
        Additive  zero  level  image  shifts to be applied.  The choices 
        are none, add the negative of the mode, median, or mean  of  the 
        specified  statistics  section,  add the values given in a file, 
        or add the values  given  by  an  image  header  keyword.   When 
        specified  in a file the zero values must be one per line in the 
        order of the input images.  File or keyword zero  offset  values 
        do not allow a correction to the weights. 
 
    weight = "none" (none|mode|median|mean|exposure|@<file>|!<keyword>) 
        Weights  to  be applied during the final averaging.  The choices 
        are none, the mode, median, or mean of the specified  statistics 
        section,  the  exposure  time, values given in a file, or values 
        given by an image header keyword.  When specified in a file  the 
        weights  must  be  one per line in the order of the input images 
        and the only adjustment made by the task is for  the  number  of 
        images  previously  combined.    In this case the weights should 
        be those appropriate for the scaled images which would  normally 
        be the inverse of the variance in the scaled image. 
 
    statsec = "" 
        Section  of  images  to  use  in  computing image statistics for 
        scaling and weighting.  If no section is given then  the  entire 
        region  of  the  input is sampled (for efficiency the images are 
        sampled if they are big enough).  When  the  images  are  offset 
        relative  to  each  other one can precede the image section with 
        one of the modifiers "input", "output",  "overlap".   The  first 
        interprets  the  section  relative  to the input image (which is 
        equivalent to not specifying a modifier), the second  interprets 
        the  section  relative to the output image, and the last selects 
        the common overlap and any following section is ignored. 
 
    expname = "" 
        Image header keyword to be used with the  exposure  scaling  and 
        weighting  options.   Also  if  an exposure keyword is specified 
        that keyword will be added to the output image using a  weighted 
        average of the input exposure values. 
 
                            Algorithm Parameters 
 
    lthreshold = INDEF, hthreshold = INDEF 
        Low  and  high  thresholds  to  be  applied to the input pixels. 
        This is done before any scaling, rejection, and  combining.   If 
        INDEF the thresholds are not used. 
 
    nlow = 1,  nhigh = 1 (minmax) 
        The  number  of  low  and  high  pixels  to  be  rejected by the 
        "minmax" algorithm.  These numbers are  converted  to  fractions 
        of  the  total  number  of input images so that if no rejections 
        have taken place the specified number  of  pixels  are  rejected 
        while  if pixels have been rejected by masking, thresholding, or 
        nonoverlap, then the fraction of the remaining pixels, truncated 
        to an integer, is used. 
 
    nkeep = 1 
        The  minimum number of pixels to retain or the maximum number to 
        reject when using the clipping  algorithms  (ccdclip,  crreject, 
        sigclip,  avsigclip,  or pclip).  When given as a positive value 
        this is the minimum number to keep.  When given  as  a  negative 
        value  the  absolute value is the maximum number to reject.  The 
        latter is in addition to pixels missing due  to  non-overlapping 
        offsets, bad pixel masks, or thresholds. 
 
    mclip = yes (ccdclip, crreject, sigclip, avsigcliip) 
        Use  the  median  as  the estimate for the true intensity rather 
        than the average with  high  and  low  values  excluded  in  the 
        "ccdclip",  "crreject",  "sigclip",  and "avsigclip" algorithms? 
        The median is a better estimator in the presence of  data  which 
        one  wants  to  reject than the average.  However, computing the 
        median is slower than the average. 
 
    lsigma = 3., hsigma = 3. (ccdclip, crreject, sigclip, avsigclip, pclip) 
        Low  and  high  sigma  clipping  factors  for   the   "ccdclip", 
        "crreject",  "sigclip",  "avsigclip",  and  "pclip"  algorithms. 
        They multiply a "sigma" factor  produced  by  the  algorithm  to 
        select  a  point below and above the average or median value for 
        rejecting  pixels.   The  lower  sigma  is   ignored   for   the 
        "crreject" algorithm. 
 
    rdnoise = "0.", gain = "1.", snoise = "0." (ccdclip, crreject) 
        CCD  readout  noise  in  electrons,  gain  in  electrons/DN, and 
        sensitivity noise as a  fraction.   These  parameters  are  used 
        with  the  "ccdclip"  and "crreject" algorithms.  The values may 
        be either numeric or an image header keyword which contains  the 
        value.  The noise model for a pixel is: 
 
            variance in DN = (rdnoise/gain)^2 + DN/gain + (snoise*DN)^2 
            variance in e- = (rdnoise)^2 + (gain*DN) + (snoise*(gain*DN))^2 
                           = rdnoise^2 + Ne + (snoise * Ne)^2 
 
        where  DN  is the data number and Ne is the number of electrons. 
        Sensitivity noise typically comes from noise  introduced  during 
        flat fielding. 
 
    sigscale = 0.1 (ccdclip, crreject, sigclip, avsigclip) 
        This  parameter  determines when poisson corrections are made to 
        the computation of a  sigma  for  images  with  different  scale 
        factors.   If all relative scales are within this value of unity 
        and all relative zero level offsets are within this fraction  of 
        the  mean  then  no correction is made.  The idea is that if the 
        images are all similarly  though  not  identically  scaled,  the 
        extra  computations  involved  in making poisson corrections for 
        variations in the sigmas can be skipped.  A value of  zero  will 
        apply  the  corrections except in the case of equal images and a 
        large value can be used if the sigmas of pixels  in  the  images 
        are independent of scale and zero level. 
 
    pclip = -0.5 (pclip) 
        Percentile  clipping  algorithm  parameter.  If greater than one 
        in absolute value then it specifies a number of pixels above  or 
        below  the  median  to use for computing the clipping sigma.  If 
        less than one in absolute value then it specifies  the  fraction 
        of  the  pixels  above  or  below the median to use.  A positive 
        value selects a point above the  median  and  a  negative  value 
        selects  a  point below the median.  The default of -0.5 selects 
        approximately the quartile point.  See the  DESCRIPTION  section 
        for further details. 
 
    grow = 0. 
        Radius  in  pixels  for  additional  pixel  to be rejected in an 
        image  with  a  rejected  pixel  from  one  of   the   rejection 
        algorithms.   This applies only to pixels rejected by one of the 
        rejection algorithms and not the masked  or  threshold  rejected 
        pixels. 
 
                           Environment Variables 
 
 
    imcombine_option (default = 1) 
        This   environment   variable   is   used   to   select  certain 
        experimental or diagnostic options.  If this  variable  has  the 
        value  1,  the default when the variable is undefined, then when 
        the number of images exceeds the number of  files  that  can  be 
        kept  open  under  IRAF  (currently  this  means  more than 4000 
        images) the images are closed and opened as needed.  This is  in 
        contrast  to  the  previous  method,  when  the variable has the 
        value 0, which first builds a single stacked image of  a  higher 
        dimension  from  the  input  images.   This  method requires the 
        images all have the same size and also that there be  sufficient 
        disk  space  for  the  stacked image and that the image  be less 
        than 2Gb in size. 
 
    irafdir 
        The IRAF home directory. 
 
    Returns
    -------
    Combined images 
 
    Example
    -------
    
    imcombine(input,output)
 
 
 
                                   I R A F 
                     Image Reduction and Analysis Facility 
    PACKAGE = immatch 
       TASK = imcombine 
 
    input   =      @outfits_8.list  List of images to combine 
    output  =          allf_8.fits  List of output images 
    (headers=                     ) List of header files (optional) 
    (bpmasks=                     ) List of bad pixel masks (optional) 
    (rejmask=                     ) List of rejection masks (optional) 
    (nrejmas=                     ) List of number rejected masks (optional) 
    (expmask=                     ) List of exposure masks (optional) 
    (sigmas =                     ) List of sigma images (optional) 
    (logfile=               STDOUT) Log file 
    
    (combine=              average) Type of combine operation 
    (reject =            avsigclip) Type of rejection 
    (project=                   no) Project highest dimension of input images? 
    (outtype=                 real) Output image pixel datatype 
    (outlimi=                     ) Output limits (x1 x2 y1 y2 ...) 
    (offsets=                 none) Input image offsets 
    (masktyp=                 none) Mask type 
    (maskval=                   0.) Mask value 
    (blank  =                   0.) Value if there are no pixels 
    
    (scale  =                 none) Image scaling 
    (zero   =                 none) Image zero point offset 
    (weight =             @weights) Image weights 
    (statsec=                 none) Image section for computing statistics 
    (expname=                 none) Image header exposure time keyword 
 
    (lthresh=                INDEF) Lower threshold 
    (hthresh=                INDEF) Upper threshold 
    (nlow   =                    1) minmax: Number of low pixels to reject 
    (nhigh  =                    1) minmax: Number of high pixels to reject 
    (nkeep  =                    1) Minimum to keep (pos) or maximum to reject (neg) 
    (mclip  =                  yes) Use median in sigma clipping algorithms? 
    (lsigma =                   3.) Lower sigma clipping factor 
    (hsigma =                   3.) Upper sigma clipping factor 
    (rdnoise=             !rdnoise) ccdclip: CCD readout noise (electrons) 
    (gain   =                !gain) ccdclip: CCD gain (electrons/DN) 
    (snoise =                   0.) ccdclip: Sensitivity noise (fraction) 
    (sigscal=                  0.1) Tolerance for sigma clipping scaling corrections 
    (pclip  =                 -0.5) pclip: Percentile clipping parameter 
    (grow   =                   0.) Radius (pixels) for neighbor rejection 
    (mode   =                   ql) 
 
 
    By D. Nidever    February 2008 
    Translated to Python by D. Nidever, April 2022
    """                     
     
    # Important directories 
    if irafdir is None:
        irafdir = os.path.expanduser('~/iraf/')
    irafdir = glob(irafdir)
    if len(irafdir) == 0:
        raise ValueError('NO IRAF DIRECTORY')
    irafdir = dln.first_el(irafdir)
    curdir = os.getcwd()
              
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if check(irafdir)==False:
        raise ValueError('IRAF TEST FAILED')
     
    # Input strings 
    sinput = str(input).strip()
    soutput = str(output).strip()
    sheaders = str(headers).strip()
    sbpmasks = str(bpmasks).strip()
    srejmask = str(rejmask).strip()
    snrejmasks = str(nrejmasks).strip()
    sexpmasks = str(expmasks).strip()
    ssigma = str(sigma).strip()
    slogfile = str(logfile).strip()
    scombine = str(combine).strip()
    sreject = str(reject).strip()
    sproject = str(project).strip()
    souttype = str(outtype).strip()
    soutlimits = str(outlimits).strip()
    soffsets = str(offsets).strip()
    smasktype = str(masktype).strip()
    smaskvalue = str(maskvalue).strip()
    sblank = str(blank).strip()
    sscale = str(scale).strip()
    szero = str(zero).strip()
    sweight = str(weight).strip()
    sstatsec = str(statsec).strip()
    sexpname = str(expname).strip()
    slthreshold = str(lthreshold).strip()
    shthreshold = str(hthreshold).strip()
    snlow = str(nlow).strip()
    snhigh = str(nhigh).strip()
    snkeep = str(nkeep).strip()
    smclip = str(mclip).strip()
    slsigma = str(lsigma).strip()
    shsigma = str(hsigma).strip()
    srdnoise = str(rdnoise).strip()
    sgain = str(gain).strip()
    ssnoise = str(snoise).strip()
    ssigscale = str(sigscale).strip()
    spclip = str(pclip).strip()
    sgrow = str(grow).strip()
     
     
    # Write IRAF script
    cmd = []
    cmd += ['print("")']  # first line will be ignored 
    cmd += ['cd '+curdir]
    cmd += ['reset min_lenuserarea = 200000']
    cmd += ['flpr']
    cmd += ['immatch']
    cmd += ['imcombine.headers = "'+sheaders+'"']
    cmd += ['imcombine.bpmasks = "'+sbpmasks+'"']
    cmd += ['imcombine.rejmask = "'+srejmask+'"']
    cmd += ['imcombine.nrejmasks = "'+snrejmasks+'"']
    cmd += ['imcombine.expmasks = "'+sexpmasks+'"']
    cmd += ['imcombine.sigma = "'+ssigma+'"']
    cmd += ['imcombine.logfile = "'+slogfile+'"']
    cmd += ['imcombine.combine = "'+scombine+'"']
    cmd += ['imcombine.reject = "'+sreject+'"']
    cmd += ['imcombine.project = '+sproject]
    cmd += ['imcombine.outtype = "'+souttype+'"']
    cmd += ['imcombine.outlimits = "'+soutlimits+'"']
    cmd += ['imcombine.offsets = "'+soffsets+'"']
    cmd += ['imcombine.masktype = "'+smasktype+'"'] 
    cmd += ['imcombine.maskvalue = '+smaskvalue]
    cmd += ['imcombine.blank = '+sblank]
    cmd += ['imcombine.scale = "'+sscale+'"']
    cmd += ['imcombine.zero = "'+szero+'"']
    cmd += ['imcombine.weight = "'+sweight+'"']
    cmd += ['imcombine.statsec = "'+sstatsec+'"'] 
    cmd += ['imcombine.expname = "'+sexpname+'"']
    cmd += ['imcombine.lthreshold = '+slthreshold] 
    cmd += ['imcombine.hthreshold = '+shthreshold]
    cmd += ['imcombine.nlow = '+snlow]
    cmd += ['imcombine.nhigh = '+snhigh] 
    cmd += ['imcombine.nkeep = '+snkeep]
    cmd += ['imcombine.mclip = '+smclip]
    cmd += ['imcombine.lsigma = '+slsigma] 
    cmd += ['imcombine.hsigma = '+shsigma]
    cmd += ['imcombine.rdnoise = "'+srdnoise+'"']
    cmd += ['imcombine.gain = "'+sgain+'"']
    cmd += ['imcombine.snoise = "'+ssnoise+'"']
    cmd += ['imcombine.sigscale = '+ssigscale]
    cmd += ['imcombine.pclip = '+spclip]
    cmd += ['imcombine("'+sinput+'","'+soutput+'")']
    cmd += ['logout']
    tid,cmdfile = tempfile.mkstemp(prefix="temp")  # absolute filename
    dln.writelines(cmdfile,cmd)
                      
    # Goto the IRAF directory
    curdir = os.getcwd()
    os.chdir(irafdir)
                      
    # Running IRAF 
    out = subprocess.check_output('cl < '+cmdfile,shell=True)
    if type(out) is bytes:
        out = out.decode()
    out = out.split('\n')   # split into lines
     
    #print,out 
     
    # The output
    lo = dln.grep(out,'IMCOMBINE',index=True)  # where to start output
    if len(lo)==0:
        lo = 0
    else:
        lo = lo[0]
    for l in out[lo:]: print(l)
     
    # Return to original directory 
    os.chdir(curdir)
     
    # Erasing the temporary files
    if os.path.exists(cmdfile): os.remove(cmdfile)

 
def imshift(input,output,xshift,yshift,shifts_file=None,interp_type='spline3',
            boundary_type='constant',constant=0.0,verbose=True,
            irafdir=None):
    """
    This runs IRAF's IMSHIFT that shifts images 
 
    Parameters
    ----------
    input 
        List of images to be transformed. 
 
    output 
        List of output images. 
 
    xshift, yshift 
        Fractional  pixel shift in x and y such that xout = xin + xshift 
        and yout = yin + yshift. 
 
    shifts_file = "" 
        The name of the text file containing the shifts for  each  input 
        image.  If  no file name is supplied each input image is shifted 
        by xshift and yshift. Shifts are listed in the text file, 1  set 
        of  shifts  per image, with the x and y shift in columns 1 and 2 
        respectively. The number of shifts in the file  must  equal  the 
        number of input images. 
 
    interp_type = "linear" 
        The  interpolant  type use to computed the output shifted image. 
        The choices are the following: 
 
        nearest 
            nearest neighbor. 
 
        linear 
            bilinear interpolation in x and y. 
 
        poly3 
            third order interior polynomial in x and y. 
 
        poly5 
            fifth order interior polynomial in x and y. 
 
        spline3 
            bicubic spline. 
 
        sinc 
            2D  sinc  interpolation.  Users   can   specify   the   sinc 
            interpolant   width  by  appending  a  width  value  to  the 
            interpolant string, e.g. sinc51 specifies a 51 by  51  pixel 
            wide  sinc  interpolant.  The  sinc  width input by the user 
            will be rounded up to the nearest odd  number.  The  default 
            sinc width is 31 by 31. 
 
        drizzle 
            2D  drizzle  resampling. Users can specify the drizzle pixel 
            fractions in x and y by appending  values  between  0.0  and 
            1.0  in  square  brackets  to  the  interpolant string, e.g. 
            drizzle[0.5]. The default value is 1.0.  The  value  0.0  is 
            increased to 0.001. Drizzle resampling with a pixel fraction 
            of 1.0 in x and y is identical to bilinear interpolation. 
 
    boundary_type = "nearest" 
        The choices are: 
 
        nearest 
            Use the value of the nearest boundary pixel. 
 
        constant 
            Use a constant value. 
 
        reflect 
            Generate value by reflecting about the boundary. 
 
        wrap 
            Generate a value by wrapping around to the opposite side  of 
            the image. 
 
    constant    The constant value to use for boundary_type="constant". 
    irafdir 
        The IRAF home diretory. 
 
 
    Returns
    -------
    Shifted images with the names in the "output" file. 

    Examle
    ------

    imshift('@inlist','@outlist',shifts_file='file.shift')
 
 
                                    I R A F 
                     Image Reduction and Analysis Facility 
    PACKAGE = immatch 
       TASK = imshift 
    
    input   =                       Input images to be fit 
    output  =                       Output images 
    xshift  =                       Fractional pixel shift in x 
    yshift  =                       Fractional pixel shift in y 
    (shifts_=                     ) Text file containing shifts for each image 
    (interp_=               linear) Interpolant (nearest,linear,poly3,poly5,spline3,sinc,drizzle) 
    (boundar=              nearest) Boundary (constant,nearest,reflect,wrap) 
    (constan=                   0.) Constant for boundary extension 
    (mode   =                   ql) 
 
    By D. Nidever    February 2013 
    Translated to Python by D. Nidever,  April 2022
    """
                      
    # xshift/yshift must be scalar
    nxshift = np.array(xshift).size
    nyshift = np.array(Yshift).size    
    if nxshift > 1 or nyshift > 1: 
        raise ValueError('XSHIFT/YSHIFT must be scalar')

    # Important directories 
    if irafdir is None:
        irafdir = os.path.expanduser('~/iraf/')
    irafdir = glob(irafdir)
    if len(irafdir) == 0:
        raise ValueError('NO IRAF DIRECTORY')
    irafdir = dln.first_el(irafdir)
    curdir = os.getcwd()
              
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if check(irafdir)==False:
        raise ValueError('IRAF TEST FAILED')

     
    # Multiple files input 
    if len(shifts_file) > 0: 
        if input[0]=='@':
            inpnames = dln.readlines(input[1:])
        if output[0]=='@':
            outnames = dln.readlines(output[1:])
        xyshifts = dln.readlines(shifts_file)
        xshifts = []
        yshifts = []
        for l in xyshifts:
            arr = l.split()
            xshifts.append(arr[0])
            ysfhits.append(arr[1])
        xshifts = np.array(xshifts).astype(float)
        yshifts = np.array(yshifts).astype(float)              
    # One file 
    else: 
        inpnames = input 
        outnames = output 
    nfiles = len(inpnames) 
     
    # Input strings 
    sinterp_type = str(interp_type).strip()
    sboundary_type = str(boundary_type).strip()
    sconstant = str(constant).strip()
     
    # Write IRAF script
    cmd = []
    cmd += ['print("")'] # first line will be ignored 
    cmd += ['cd '+curdir]
    cmd += ['immatch']
    cmd += ['imshift.input=""']
    cmd += ['imshift.output=""'] 
    cmd += ['imshift.shifts_file = ""']  # don't use shifts_file use x/yshifts 
    cmd += ['imshift.interp_type = "'+sinterp_type+'"']
    cmd += ['imshift.boundary_type = "'+sboundary_type+'"']
    cmd += ['imshift.constant = '+sconstant]
     
    # Loop through files 
    for i in range(nfiles): 
        sinpname1 = str(inpnames[i]).strip()
        soutname1 = str(outnames[i]).strip()
        # if shift is very small round to zero, otherwise can cause problems 
        if abs(xshift[i]) < 1e-4: 
            xshift1 = 0.0 
        else: 
            xshift1 = xshift[i] 
        if abs(yshift[i]) < 1e-4: 
            yshift1 = 0.0 
        else: 
            yshift1 = yshift[i] 
        sxshift1 = str(xshift1) 
        syshift1 = str(yshift1) 
        cmd += ['imshift("'+sinpname1+'","'+soutname1+'",xshift='+sxshift1+',yshift='+syshift1+')']
        cmd += ['flprcache']
     
    cmd += ['logout']
    tid,cmdfile = tempfile.mkstemp(prefix="temp") # absolute filename
    dln.writelines(cmdfile,cmd)
     
     
    # Goto the IRAF directory
    curdir = os.getcwd()
    os.chdir(irafdir)
     
    # Running IRAF 
    out = subprocess.check_outut('cl < '+cmdfile,shell=True)

    if type(out) is bytes:
        out = out.decode()
    out = out.split('\n')   # split into lines
     
    if verbose:
        for l in out: print(l)
     
    # Return to original directory
    os.chdir(curdir)
     
    # Erasing the temporary files
    if os.path.exists(cmdfile): os.remove(cmdfile)


