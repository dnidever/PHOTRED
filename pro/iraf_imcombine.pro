pro iraf_imcombine,input,output,headers=headers,bpmasks=bpmasks,rejmask=rejmask,$
              nrejmasks=nrejmasks,expmasks=expmasks,sigma=sigma,logfile=logfile,$
              combine=combine,reject=reject,project=project,outtype=outtype,$
              outlimits=outlimits,offsets=offsets,masktype=masktype,$
              maskvalue=maskvalue,blank=blank,scale=scale,zero=zero,weight=weight,$
              statsec=statsec,expname=expname,lthreshold=lthreshold,$
              hthreshold=hthreshold,nlow=nlow,nhigh=nhigh,nkeep=nkeep,mclip=mclip,$
              lsigma=lsigma,hsigma=hsigma,rdnoise=rdnoise,gain=gain,snoise=snoise,$
              sigscale=sigscale,pclip=pclip,grow=grow,irafdir=irafdir,error=error

;+
;
; IRAF_IMCOMBINE
; 
; PURPOSE:
;  This runs IRAF's IMCOMBINE that combines images
;
; INPUTS:
;
;    input
;        List of input images to combine.  If the  project  parameter  is
;        "no"  then  all  input  images must have the same dimensionality
;        though they may be of different  sizes.   Otherwise  each  input
;        image   is  handled  separately  and  they  may  have  different 
;        dimensionalities.
;    
;    
;    When the project parameter is "no" all the input images are combined
;    into  a  single  output file.  In this case the following parameters
;    specify only a single file name.   Otherwise  each  input  image  is
;    combined  by  projecting (combining across) the highest dimension to
;    produce a lower dimensional  image.   For  this  type  of  combining
;    there  is  one output for each input and so the following parameters
;    specify matching lists.
;    
;    
;    output
;        Output combined image(s).  If there are  fewer  than  100  input
;        images  the  names  of  the  input images are recorded in header
;        keywords IMCMBnnn.
;    
;    headers = "" (optional)
;        Optional output multiextension  FITS  file(s).   The  extensions
;        are dataless headers from each input image.
;    
;    bpmasks = "" (optional)
;        Optional  output bad pixel mask(s) with good values of 0 and bad
;        values of 1.  Output pixels are marked  as  bad  when  no  input
;        pixels  contributed  to the output pixel.  The file name is also
;        added to the output image header under the keyword BPM.
;    
;    rejmask = "" (optional)
;        Optional output mask file(s) identifying  rejected  or  excluded
;        pixels.   The  pixel  mask  is  the size of the output image but
;        there is one extra dimension with length equal to the number  of
;        input  images.   Each element of the highest dimension is a mask
;        corresponding to an input image with values of  1  for  rejected
;        or  excluded  pixels and values of 0 for pixels which were used.
;        The order of the masks is the order  of  the  input  images  and
;        image  header  keywords,  indexed by the pixel coordinate of the
;        highest dimension identify the  input  images.   Note  that  the
;        pixel positions are in the output pixel coordinate system.
;    
;    nrejmasks = "" (optional)
;        Optional  output pixel mask(s) giving the number of input pixels
;        rejected or excluded from the input images.
;    
;    expmasks = "" (optional)
;        Optional output exposure mask(s) giving the sum of the  exposure
;        values   of   the   input  images  with  non-zero  weights  that 
;        contributed  to  that  pixel.   Since  masks  are  integer,  the 
;        exposure  values  may  be  scaled  to preserve dynamic range and
;        fractional significance.  The scaling values are  given  in  the
;        header  under  the  keywords  MASKSCAL  and  MASKZERO.  Exposure
;        values are computed from the mask values  by  scale  *  value  +
;        zero  where  scale is the value of the MASKSCAL keyword and zero
;        is the value of the MASKZERO keyword.
;    
;    sigma = "" (optional)
;        Optional output sigma  image(s).   The  sigma  is  the  standard
;        deviation,  corrected  for  a  finite  population,  of the input
;        pixel  values  (excluding  rejected  pixels)  about  the  output 
;        combined pixel values.
;    
;    
;    logfile = "STDOUT" (optional)
;        Optional  output  log file.  If no file is specified then no log
;        information is produced.  The special filename  "STDOUT"  prints
;        log information to the terminal.
;    
;    
;    combine = "average" (average|median|sum)
;        Type  of  combining  operation  performed  on  the  final set of
;        pixels   (after   offsetting,   masking,    thresholding,    and  
;        rejection).   The  choices  are  "average",  "median", or "sum".
;        The median uses the average of the two central values  when  the
;        number  of  pixels  is even.  For the average and sum, the pixel
;        values are multiplied by the weights (1 if no weighting is used)
;        and  summed.   The average is computed by dividing by the sum of
;        the weights.  If the  sum  of  the  weights  is  zero  then  the
;        unweighted average is used.
;    
;    reject = "none" (none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip)
;        Type  of  rejection  operation performed on the pixels remaining
;        after offsetting, masking and thresholding.  The algorithms  are
;        described  in  the  DESCRIPTION  section.  The rejection choices
;        are:
;                
;              none - No rejection
;            minmax - Reject the nlow and nhigh pixels
;           ccdclip - Reject pixels using CCD noise parameters
;          crreject - Reject only positive pixels using CCD noise parameters
;           sigclip - Reject pixels using a sigma clipping algorithm
;         avsigclip - Reject pixels using an averaged sigma clipping algorithm
;             pclip - Reject pixels using sigma based on percentiles
;        
;    project = no
;        Project (combine) across the  highest  dimension  of  the  input
;        images?   If  "no"  then all  the input images are combined to a
;        single output  image.   If  "yes"  then  the  highest  dimension
;        elements  of  each  input  image are combined to an output image
;        and optional pixel list and sigma images.  Each element  of  the
;        highest dimension may have a separate offset.
;    
;    outtype = "real" (none|short|ushort|integer|long|real|double)
;        Output  image pixel datatype.  The pixel datatypes are "double",
;        "real", "long", "integer", unsigned short "ushort", and  "short"
;        with  highest precedence first.  If "none" is specified then the
;        highest precedence datatype of the input images is  used.   When
;        there  is  a  mixture  of  short  and  unsigned short images the
;        highest  precedence  become  integer.   The  datatypes  may   be 
;        abbreviated to a single character.
;    
;    outlimits = ""
;        Output  region limits specified as pairs of whitespace separated
;        values.  The first two numbers are the limits  along  the  first
;        output  image  dimension,  the  next  two numbers are the limits
;        along the second dimension, and so on.  If the higher  dimension
;        limits  are  not  specified  they  default  to  the  full range.
;        Therefore, if no limits are specified then the  full  output  is
;        created.   Note  that  the  output size is computed from all the
;        input images including offsets if specified and the  coordinates
;        are relative to that size.
;    
;    
;    offsets = "none" (none|wcs|world|physical|grid|<filename>)
;        Integer offsets to add to each image axes.  The options are:
;        
;        "none"
;            No offsets are applied.
;        
;        "wcs" or "world"
;            The  world  coordinate  system (wcs) in the image is used to
;            derive  the  offsets.   The  nearest  integer  offset   that 
;            matches  the  world  coordinate  at  the center of the first
;            input image is used.
;        
;        "physical"
;            The physical coordinate system defined by the  IRAF  LTM/LTV
;            keywords define the offsets.
;        
;        "grid"
;            A  uniform  grid  of offsets is specified by a string of the
;            form
;            
;                    grid [n1] [s1] [n2] [s2] ...
;            
;            where ni is the number of images in dimension i  and  si  is
;            the  step  in  dimension  i.  For example "grid 5 100 5 100"
;            specifies a 5x5 grid with origins offset by 100 pixels.
;        
;        <filename>
;            The offsets are given  in  the  specified  file.   The  file
;            consists  of  one  line  per  image with the offsets in each
;            dimension forming the columns.
;    
;    masktype = "none" (none|goodvalue|badvalue|goodbits|badbits|!<keyword>)
;        Type of pixel masking to use.  If "none" or  ""  then  no  pixel
;        masking  is  done  even  if  an  image  has an associated  pixel
;        mask.  The other choices are to select the value  in  the  pixel
;        mask  to be treated as good (goodvalue) or bad (badvalue) or the
;        bits (specified as a value) to be treated as good (goodbits)  or
;        bad  (badbits).   In  these cases the pixel mask file name comes
;        from the image header keyword BPM.  If the  parameter  is  given
;        as  "!<keyword>"  where  <keyword> is a header keyword, the mask
;        file comes from the value of that keyword  and  the  mask  value
;        interpretation  is the same as "goodvalue".  Note, if the number
;        of input images becomes too large (currently about 4090 .imh  or
;        2045  projection.   This means all the images must have the same
;        size and dimensionality.
;    
;    maskvalue = 0
;        Mask value used with the masktype parameter.  If the  mask  type
;        selects  good  or bad bits the value may be specified using IRAF
;        notation for decimal, octal, or hexidecimal; i.e  12,  14b,  0cx
;        to select bits 3 and 4.
;    
;    blank = 0.
;        Output value to be used when there are no pixels.
;    
;    scale = "none" (none|mode|median|mean|exposure|@<file>|!<keyword>)
;        Multiplicative  image  scaling  to  be applied.  The choices are
;        none, multiply by the reciprocal of the mode,  median,  or  mean
;        of  the specified statistics section, multiply by the reciprocal
;        of the exposure time  in  the  image  header,  multiply  by  the
;        values  in  a  specified  file, or multiply by a specified image
;        header keyword.  When specified in a file  the  scales  must  be
;        one per line in the order of the input images.
;    
;    zero = "none" (none|mode|median|mean|@<file>|!<keyword>)
;        Additive  zero  level  image  shifts to be applied.  The choices
;        are none, add the negative of the mode, median, or mean  of  the
;        specified  statistics  section,  add the values given in a file,
;        or add the values  given  by  an  image  header  keyword.   When
;        specified  in a file the zero values must be one per line in the
;        order of the input images.  File or keyword zero  offset  values
;        do not allow a correction to the weights.
;    
;    weight = "none" (none|mode|median|mean|exposure|@<file>|!<keyword>)
;        Weights  to  be applied during the final averaging.  The choices
;        are none, the mode, median, or mean of the specified  statistics
;        section,  the  exposure  time, values given in a file, or values
;        given by an image header keyword.  When specified in a file  the
;        weights  must  be  one per line in the order of the input images
;        and the only adjustment made by the task is for  the  number  of
;        images  previously  combined.    In this case the weights should
;        be those appropriate for the scaled images which would  normally
;        be the inverse of the variance in the scaled image.
;    
;    statsec = ""
;        Section  of  images  to  use  in  computing image statistics for
;        scaling and weighting.  If no section is given then  the  entire
;        region  of  the  input is sampled (for efficiency the images are
;        sampled if they are big enough).  When  the  images  are  offset
;        relative  to  each  other one can precede the image section with
;        one of the modifiers "input", "output",  "overlap".   The  first
;        interprets  the  section  relative  to the input image (which is
;        equivalent to not specifying a modifier), the second  interprets
;        the  section  relative to the output image, and the last selects
;        the common overlap and any following section is ignored.
;    
;    expname = ""
;        Image header keyword to be used with the  exposure  scaling  and
;        weighting  options.   Also  if  an exposure keyword is specified
;        that keyword will be added to the output image using a  weighted
;        average of the input exposure values.
;    
;                            Algorithm Parameters
;    
;    lthreshold = INDEF, hthreshold = INDEF
;        Low  and  high  thresholds  to  be  applied to the input pixels.
;        This is done before any scaling, rejection, and  combining.   If
;        INDEF the thresholds are not used.
;    
;    nlow = 1,  nhigh = 1 (minmax)
;        The  number  of  low  and  high  pixels  to  be  rejected by the
;        "minmax" algorithm.  These numbers are  converted  to  fractions
;        of  the  total  number  of input images so that if no rejections
;        have taken place the specified number  of  pixels  are  rejected
;        while  if pixels have been rejected by masking, thresholding, or
;        nonoverlap, then the fraction of the remaining pixels, truncated
;        to an integer, is used.
;    
;    nkeep = 1
;        The  minimum number of pixels to retain or the maximum number to
;        reject when using the clipping  algorithms  (ccdclip,  crreject,
;        sigclip,  avsigclip,  or pclip).  When given as a positive value
;        this is the minimum number to keep.  When given  as  a  negative
;        value  the  absolute value is the maximum number to reject.  The
;        latter is in addition to pixels missing due  to  non-overlapping
;        offsets, bad pixel masks, or thresholds.
;    
;    mclip = yes (ccdclip, crreject, sigclip, avsigcliip)
;        Use  the  median  as  the estimate for the true intensity rather
;        than the average with  high  and  low  values  excluded  in  the
;        "ccdclip",  "crreject",  "sigclip",  and "avsigclip" algorithms?
;        The median is a better estimator in the presence of  data  which
;        one  wants  to  reject than the average.  However, computing the
;        median is slower than the average.
;    
;    lsigma = 3., hsigma = 3. (ccdclip, crreject, sigclip, avsigclip, pclip)
;        Low  and  high  sigma  clipping  factors  for   the   "ccdclip", 
;        "crreject",  "sigclip",  "avsigclip",  and  "pclip"  algorithms. 
;        They multiply a "sigma" factor  produced  by  the  algorithm  to
;        select  a  point below and above the average or median value for
;        rejecting  pixels.   The  lower  sigma  is   ignored   for   the 
;        "crreject" algorithm.
;    
;    rdnoise = "0.", gain = "1.", snoise = "0." (ccdclip, crreject)
;        CCD  readout  noise  in  electrons,  gain  in  electrons/DN, and
;        sensitivity noise as a  fraction.   These  parameters  are  used
;        with  the  "ccdclip"  and "crreject" algorithms.  The values may
;        be either numeric or an image header keyword which contains  the
;        value.  The noise model for a pixel is:
;        
;            variance in DN = (rdnoise/gain)^2 + DN/gain + (snoise*DN)^2
;            variance in e- = (rdnoise)^2 + (gain*DN) + (snoise*(gain*DN))^2
;                           = rdnoise^2 + Ne + (snoise * Ne)^2
;        
;        where  DN  is the data number and Ne is the number of electrons.
;        Sensitivity noise typically comes from noise  introduced  during
;        flat fielding.
;    
;    sigscale = 0.1 (ccdclip, crreject, sigclip, avsigclip)
;        This  parameter  determines when poisson corrections are made to
;        the computation of a  sigma  for  images  with  different  scale
;        factors.   If all relative scales are within this value of unity
;        and all relative zero level offsets are within this fraction  of
;        the  mean  then  no correction is made.  The idea is that if the
;        images are all similarly  though  not  identically  scaled,  the
;        extra  computations  involved  in making poisson corrections for
;        variations in the sigmas can be skipped.  A value of  zero  will
;        apply  the  corrections except in the case of equal images and a
;        large value can be used if the sigmas of pixels  in  the  images
;        are independent of scale and zero level.
;    
;    pclip = -0.5 (pclip)
;        Percentile  clipping  algorithm  parameter.  If greater than one
;        in absolute value then it specifies a number of pixels above  or
;        below  the  median  to use for computing the clipping sigma.  If
;        less than one in absolute value then it specifies  the  fraction
;        of  the  pixels  above  or  below the median to use.  A positive
;        value selects a point above the  median  and  a  negative  value
;        selects  a  point below the median.  The default of -0.5 selects
;        approximately the quartile point.  See the  DESCRIPTION  section
;        for further details.
;    
;    grow = 0.
;        Radius  in  pixels  for  additional  pixel  to be rejected in an
;        image  with  a  rejected  pixel  from  one  of   the   rejection 
;        algorithms.   This applies only to pixels rejected by one of the
;        rejection algorithms and not the masked  or  threshold  rejected
;        pixels.
;    
;                           Environment Variables
;    
;    
;    imcombine_option (default = 1)
;        This   environment   variable   is   used   to   select  certain 
;        experimental or diagnostic options.  If this  variable  has  the
;        value  1,  the default when the variable is undefined, then when
;        the number of images exceeds the number of  files  that  can  be
;        kept  open  under  IRAF  (currently  this  means  more than 4000
;        images) the images are closed and opened as needed.  This is  in
;        contrast  to  the  previous  method,  when  the variable has the
;        value 0, which first builds a single stacked image of  a  higher
;        dimension  from  the  input  images.   This  method requires the
;        images all have the same size and also that there be  sufficient
;        disk  space  for  the  stacked image and that the image  be less
;        than 2Gb in size.
;
;  =irafdir
;        The IRAF home directory.
;
; OUTPUTS:
;  Combined images
;  =error   The error message if there was one.
;
; USAGE:
;  IDL>Iraf_imcombine,input,output
;
;
;
;                                   I R A F  
;                     Image Reduction and Analysis Facility
; PACKAGE = immatch
;    TASK = imcombine
;
; input   =      @outfits_8.list  List of images to combine
; output  =          allf_8.fits  List of output images
; (headers=                     ) List of header files (optional)
; (bpmasks=                     ) List of bad pixel masks (optional)
; (rejmask=                     ) List of rejection masks (optional)
; (nrejmas=                     ) List of number rejected masks (optional)
; (expmask=                     ) List of exposure masks (optional)
; (sigmas =                     ) List of sigma images (optional)
; (logfile=               STDOUT) Log file
;
; (combine=              average) Type of combine operation
; (reject =            avsigclip) Type of rejection
; (project=                   no) Project highest dimension of input images?
; (outtype=                 real) Output image pixel datatype
; (outlimi=                     ) Output limits (x1 x2 y1 y2 ...)
; (offsets=                 none) Input image offsets
; (masktyp=                 none) Mask type
; (maskval=                   0.) Mask value
; (blank  =                   0.) Value if there are no pixels
;
; (scale  =                 none) Image scaling
; (zero   =                 none) Image zero point offset
; (weight =             @weights) Image weights
; (statsec=                 none) Image section for computing statistics
; (expname=                 none) Image header exposure time keyword
;
; (lthresh=                INDEF) Lower threshold
; (hthresh=                INDEF) Upper threshold
; (nlow   =                    1) minmax: Number of low pixels to reject
; (nhigh  =                    1) minmax: Number of high pixels to reject
; (nkeep  =                    1) Minimum to keep (pos) or maximum to reject (neg)
; (mclip  =                  yes) Use median in sigma clipping algorithms?
; (lsigma =                   3.) Lower sigma clipping factor
; (hsigma =                   3.) Upper sigma clipping factor
; (rdnoise=             !rdnoise) ccdclip: CCD readout noise (electrons)
; (gain   =                !gain) ccdclip: CCD gain (electrons/DN)
; (snoise =                   0.) ccdclip: Sensitivity noise (fraction)
; (sigscal=                  0.1) Tolerance for sigma clipping scaling corrections
; (pclip  =                 -0.5) pclip: Percentile clipping parameter
; (grow   =                   0.) Radius (pixels) for neighbor rejection
; (mode   =                   ql)
;
;
;
; By D. Nidever    February 2008 
;-

undefine,error

; Not enough inputs
ninput = n_elements(input)
noutput = n_elements(output)
if ninput eq 0 or noutput eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - iraf_imcombine,input,output,headers=headers,bpmasks=bpmasks,rejmask=rejmask,'
  print,'              nrejmasks=nrejmasks,expmasks=expmasks,sigma=sigma,logfile=logfile,'
  print,'              combine=combine,reject=reject,project=project,outtype=outtype,'
  print,'              outlimits=outlimits,offsets=offsets,masktype=masktype,'
  print,'              maskvalue=maskvalue,blank=blank,scale=scale,zero=zero,weight=weight,'
  print,'              statsec=statsec,expname=expname,lthreshold=lthreshold,'
  print,'              hthreshold=hthreshold,nlow=nlow,nhigh=nhigh,nkeep=nkeep,mclip=mclip,'
  print,'              lsigma=lsigma,hsigma=hsigma,rdnoise=rdnoise,gain=gain,snoise=snoise,'
  print,'              sigscale=sigscale,pclip=pclip,grow=grow,irafdir=irafdir,error=error'
  return
endif


; Important directories
;irafdir = '/net/home/dln5q/iraf/'
if n_elements(irafdir) eq 0 then irafdir='~/iraf/'
irafdir = FILE_SEARCH(irafdir,/fully_qualify,count=nirafdir)
if nirafdir eq 0 then begin
  error = 'NO IRAF DIRECTORY'
  print,'NO IRAF DIRECTORY'
  return 
endif    
CD,current=curdir

; Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL
;---------------------------------------------------------------
CHECK_IRAF,iraftest,irafdir=irafdir
if iraftest eq 0 then begin
  error = 'IRAF TEST FAILED'
  print,'IRAF TEST FAILED.  EXITING'
  return
endif


; Default values
; These are the IRAF defaults
if n_elements(headers) eq 0 then headers=''
if n_elements(bpmasks) eq 0 then bpmasks=''
if n_elements(rejmask) eq 0 then rejmask=''
if n_elements(nrejmasks) eq 0 then nrejmasks=''
if n_elements(expmasks) eq 0 then expmasks=''
if n_elements(sigma) eq 0 then sigma=''
if n_elements(logfile) eq 0 then logfile='STDOUT'
if n_elements(combine) eq 0 then combine='average'
if n_elements(reject) eq 0 then reject='none'
if n_elements(project) eq 0 then project='no'
if n_elements(outtype) eq 0 then outtype='real'
if n_elements(outlimits) eq 0 then outlimits=''
if n_elements(offsets) eq 0 then offsets='none'
if n_elements(masktype) eq 0 then masktype='none'
if n_elements(maskvalue) eq 0 then maskvalue=0.
if n_elements(blank) eq 0 then blank=0.
if n_elements(scale) eq 0 then scale='none'
if n_elements(zero) eq 0 then zero='none'
if n_elements(weight) eq 0 then weight='none'
if n_elements(statsec) eq 0 then statsec=''
if n_elements(expname) eq 0 then expname=''
if n_elements(lthreshold) eq 0 then lthreshold='INDEF'
if n_elements(hthreshold) eq 0 then hthreshold='INDEF'
if n_elements(nlow) eq 0 then nlow=1
if n_elements(nhigh) eq 0 then nhigh=1
if n_elements(nkeep) eq 0 then nkeep=1
if n_elements(mclip) eq 0 then mclip='yes'
if n_elements(lsigma) eq 0 then lsigma=3.
if n_elements(hsigma) eq 0 then hsigma=3.
if n_elements(rdnoise) eq 0 then rdnoise=0.
if n_elements(gain) eq 0 then gain=1.
if n_elements(snoise) eq 0 then snoise=0.
if n_elements(sigscale) eq 0 then sigscale=0.1
if n_elements(pclip) eq 0 then pclip=-0.5
if n_elements(grow) eq 0 then grow=0.




; Input strings
sinput = strtrim(input,2)
soutput = strtrim(output,2)
sheaders = strtrim(headers,2)
sbpmasks = strtrim(bpmasks,2)
srejmask = strtrim(rejmask,2)
snrejmasks = strtrim(nrejmasks,2)
sexpmasks = strtrim(expmasks,2)
ssigma = strtrim(sigma,2)
slogfile = strtrim(logfile,2)
scombine = strtrim(combine,2)
sreject = strtrim(reject,2)
sproject = strtrim(project,2)
souttype = strtrim(outtype,2)
soutlimits = strtrim(outlimits,2)
soffsets = strtrim(offsets,2)
smasktype = strtrim(masktype,2)
smaskvalue = strtrim(maskvalue,2)
sblank = strtrim(blank,2)
sscale = strtrim(scale,2)
szero = strtrim(zero,2)
sweight = strtrim(weight,2)
sstatsec = strtrim(statsec,2)
sexpname = strtrim(expname,2)
slthreshold = strtrim(lthreshold,2)
shthreshold = strtrim(hthreshold,2)
snlow = strtrim(nlow,2)
snhigh = strtrim(nhigh,2)
snkeep = strtrim(nkeep,2)
smclip = strtrim(mclip,2)
slsigma = strtrim(lsigma,2)
shsigma = strtrim(hsigma,2)
srdnoise = strtrim(rdnoise,2)
sgain = strtrim(gain,2)
ssnoise = strtrim(snoise,2)
ssigscale = strtrim(sigscale,2)
spclip = strtrim(pclip,2)
sgrow = strtrim(grow,2)


; Write IRAF script
push,cmd,'cd '+curdir
push,cmd,'immatch'
push,cmd,'imcombine.headers = "'+sheaders+'"'
push,cmd,'imcombine.bpmasks = "'+sbpmasks+'"'
push,cmd,'imcombine.rejmask = "'+srejmask+'"'
push,cmd,'imcombine.nrejmasks = "'+snrejmasks+'"'
push,cmd,'imcombine.expmasks = "'+sexpmasks+'"'
push,cmd,'imcombine.sigma = "'+ssigma+'"'
push,cmd,'imcombine.logfile = "'+slogfile+'"'
push,cmd,'imcombine.combine = "'+scombine+'"'
push,cmd,'imcombine.reject = "'+sreject+'"'
push,cmd,'imcombine.project = '+sproject
push,cmd,'imcombine.outtype = "'+souttype+'"'
push,cmd,'imcombine.outlimits = "'+soutlimits+'"'
push,cmd,'imcombine.offsets = "'+soffsets+'"'
push,cmd,'imcombine.masktype = "'+smasktype+'"'
push,cmd,'imcombine.maskvalue = '+smaskvalue
push,cmd,'imcombine.blank = '+sblank
push,cmd,'imcombine.scale = "'+sscale+'"'
push,cmd,'imcombine.zero = "'+szero+'"'
push,cmd,'imcombine.weight = "'+sweight+'"'
push,cmd,'imcombine.statsec = "'+sstatsec+'"'
push,cmd,'imcombine.expname = "'+sexpname+'"'
push,cmd,'imcombine.lthreshold = '+slthreshold
push,cmd,'imcombine.hthreshold = '+shthreshold
push,cmd,'imcombine.nlow = '+snlow
push,cmd,'imcombine.nhigh = '+snhigh
push,cmd,'imcombine.nkeep = '+snkeep
push,cmd,'imcombine.mclip = '+smclip
push,cmd,'imcombine.lsigma = '+slsigma
push,cmd,'imcombine.hsigma = '+shsigma
push,cmd,'imcombine.rdnoise = "'+srdnoise+'"'
push,cmd,'imcombine.gain = "'+sgain+'"'
push,cmd,'imcombine.snoise = "'+ssnoise+'"'
push,cmd,'imcombine.sigscale = '+ssigscale
push,cmd,'imcombine.pclip = '+spclip
push,cmd,'imcombine("'+sinput+'","'+soutput+'")'
push,cmd,'logout'
;cmdfile = maketemp('temp','.cl')
cmdfile = MKTEMP('temp')    ; absolute filename
WRITELINE,cmdfile,cmd

; Goto the IRAF directory
CD,current=curdir
CD,irafdir

; Running IRAF
undefine,out
;SPAWN,'cl < '+curdir+'/'+cmdfile,out,errout
SPAWN,'cl < '+cmdfile,out,errout

;print,out

; The output
lo = first_el(where(stregex(out,'IMCOMBINE',/boolean) eq 1))  ; where to start output
lo = lo > 0
printline,out[lo:*]

; Return to original directory
CD,curdir


; Erasing the temporary files
FILE_DELETE,cmdfile,/allow,/quiet


if keyword_set(stp) then stop

end
