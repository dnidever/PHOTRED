;+
;
; IRAF_IMALIGN
; 
; PURPOSE:
;  This runs IRAF's IMALIGN that aligns and shifts images.
;
; INPUTS:
;    input
;        The  input  images  to  be shifted and trimmed.  The input image
;        list should contain the reference image so that its borders  are
;        used in the computation of the overlap region.
;    
;    reference
;        The reference image to which the input images will be aligned.
;    
;    coords
;        A  text  file  containing the reference image coordinates of the
;        registration objects to be centered in each  image,  one  object
;        per  line  with  the  x and y coordinates in columns one and two
;        respectively.
;    
;    output
;        The output images.
;    
;    shifts = ""
;        A text file containing the initial estimate for  each  image  of
;        the  shift  in each axis relative to the reference image.  These
;        estimates are used to modify the coordinates of the registration
;        objects  prior  to  centering.   The  format  of the file is one
;        image per line with the x and y shifts in columns  one  and  two
;        respectively.    The   sense   of   the  shifts  is  such  that: 
;        Xshift=Xref-Xin and  Yshift=Yref-Yin.   If  shifts  is  null,  a
;        coarse  centering  pass will be made to attempt to determine the
;        initial shifts.
;    
;    boxsize = 7
;        The size in pixels of the box to use for  the  final  centering,
;        during  which  all  the sources in coords are recentered in each
;        image using the initial estimate of the relative shift for  each
;        image.   Care should be taken to choose an appropriate value for
;        this parameter, since it is highly data dependent.
;
;    bigbox = 11
;        The size in pixels of the box to use for coarse centering.   The
;        coarse  pass  through  the  centering algorithm is made with the
;        box centered at the nominal position of the first source in  the
;        coordinate  list.   Coarse  centering  is  performed only if the
;        shifts file is undefined.  Care should be  taken  to  choose  an
;        appropriate  value  for  this parameter, since it is highly data
;        dependent.  Large values  should  be  suspect  until  the  final
;        results  are  checked to see that the centering did not converge
;        on the wrong coordinates,  although  the  usual  result  for  an
;        inappropriate  bigbox  size  is  that  the  algorithm  fails  to 
;        converge and the task aborts.
;    
;    negative = no
;        Are the features negative ?
;    
;    background = INDEF
;        The  absolute  reference  level  for   the   marginal   centroid 
;        calculation.   If  background  is INDEF, this is set to the mean
;        value (between the thresholds) of the individual sources.
;    
;    lower = INDEF
;        The lower threshold for the data.  Individual pixels  less  than
;        this value will be given zero weight in the centroids.
;    
;    upper = INDEF
;        The  upper  threshold  for  the data.  Individual pixels greater
;        than this value will be given zero weight in the centroids.
;    
;    niterate = 3
;        The maximum number of  centering  iterations  to  perform.   The
;        centering  will  halt  when  this  limit  is reached or when the
;        desired Itolerance is achieved.
;    
;    tolerance = 0
;        The tolerance for convergence of the centering algorithm.   This
;        is  the  integral  shift of the centering box from one iteration
;        to the next.
;    
;    maxshift = INDEFR
;        The maximum permitted difference  between  the  predicted  shift
;        and  the the computed shift for each object. Objects with shifts
;        greater than maxshift are ignored. If maxshift is  undefined  no
;        shift checking is done.
;
;    shiftimages = yes
;        If  shiftimages  is  yes, the IMSHIFT task will be used to align
;        the images.  If shiftimages  is  no,  the  images  will  not  be
;        aligned, but the coordinates will still be centered.
;    
;    interp_type = "spline3"
;        The interpolation function used by the IMSHIFT task.
;    
;    boundary_type = "constant"
;        The boundary extension type used by the IMSHIFT task.
;    
;    constant = 0.
;        The  constant  used  by  the  IMSHIFT  task  if boundary_type is
;        "constant".
;    
;    trimimages = yes
;        If trimimages is yes, the  output  images  will  be  trimmed  to
;        include  only  the region over which they all overlap.  The trim
;        section that is actually used  may  differ  slightly  from  that
;        reported   by   IMCENTROID,  due  to  a  correction  applied  to 
;        compensate for the boundary extension "contamination"  near  the
;        edges of the images.
;    
;    verbose = yes
;        Print the centers, shifts, and trim section?
;
;   =irafdir
;        The IRAF home diretory.
;
;
; OUTPUTS:
;  Shifted images with the names in the "output" file.
;  =xoff  The offset in X between the original and shifted images
;         xorig = xshift + xoff
;  =yoff  The offset in Y between the original ahd shifted images
;         yorig = yshift + yoff
;  =trans The X/Y shifts in the format [nfiles,X/Y].
;  =trimsection  The trim section in the format [Xstart,Xend,Ystart,Yend]
;
; USAGE:
;  IDL>iraf_imalign,'@inlist','ref.fits','coordfile','@outlist'
;
;
;
;                                    I R A F  
;                     Image Reduction and Analysis Facility
; PACKAGE = immatch
;    TASK = imalign
;
; input   =         @fits_8.list  Input images
; referenc=        obj034_8.fits  Reference image
; coords  =          chip8.coord  Reference coordinates file
; output  =      @outfits_8.list  Output images
; (shifts =          chip8.shift) Initial shifts file
; (boxsize=                    7) Size of the small centering box
; (bigbox =                   11) Size of the big centering box
; (negativ=                   no) Are the features negative ?
; (backgro=                INDEF) Reference background level
; (lower  =                INDEF) Lower threshold for data
; (upper  =                INDEF) Upper threshold for data
; (niterat=                    5) Maximum number of iterations
; (toleran=                    0) Tolerance for convergence
; (maxshif=                INDEF) Maximum acceptable pixel shift
; (shiftim=                  yes) Shift the images ?
; (interp_=               linear) Interpolant
; (boundar=              nearest) Boundary type
; (constan=                   0.) Constant for constant boundary extension
; (trimima=                  yes) Trim the shifted images ?
; (verbose=                  yes) Print the centers, shifts, and trim section ?
; (list   =                     )
; (mode   =                   ql)
;
;
; By D. Nidever    February 2008 
;-

pro iraf_imalign,input,reference,coords,output,shifts=shifts,boxsize=boxsize,$
                bigbox=bigbox,negative=negative,background=background,$
                lower=lower,upper=upper,niterate=niterate,tolerance=tolerance,$
                maxshift=maxshift,shiftimages=shiftimages,interp_type=interp_type,$
                boundary_type=boundary_type,constant=constant,trimimages=trimimages,$
                verbose=verbose,xoff=xoff,yoff=yoff,stp=stp,error=error,trans=trans,$
                trimsection=trimsection,irafdir=irafdir

undefine,error

ninput = n_elements(input)
nreference = n_elements(reference)
ncoords = n_elements(coords)
noutput = n_elements(output)

; Not enough inputs
if ninput eq 0 or nreference eq 0 or ncoords eq 0 or noutput eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax -  iraf_imalign,input,reference,coords,output,shifts=shifts,boxsize=boxsize,'
  print,'                bigbox=bigbox,negative=negative,background=background,'
  print,'                lower=lower,upper=upper,niterate=niterate,tolerance=tolerance,'
  print,'                maxshift=maxshift,shiftimages=shiftimages,interp_type=interp_type,'
  print,'                boundary_type=boundary_type,constant=constant,trimimages=trimimages,'
  print,'                verbose=verbose'
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
if n_elements(shifts) eq 0 then shifts=''
if n_elements(boxsize) eq 0 then boxsize=7
if n_elements(bigbox) eq 0 then bigbox=11
if n_elements(negative) eq 0 then negative='no'
if n_elements(background) eq 0 then background='INDEF'
if n_elements(lower) eq 0 then lower='INDEF'
if n_elements(upper) eq 0 then upper='INDEF'
if n_elements(niterate) eq 0 then niterate=3
if n_elements(tolerance) eq 0 then tolerance=0
if n_elements(maxshift) eq 0 then maxshift='INDEF'
if n_elements(shiftimages) eq 0 then shiftimages='yes'
if n_elements(interp_type) eq 0 then interp_type='spline3'
if n_elements(boundary_type) eq 0 then boundary_type='constant'
if n_elements(constant) eq 0 then constant=0.
if n_elements(trimimages) eq 0 then trimimages='yes'
if n_elements(verbose) eq 0 then verbose='yes'




; Input strings
sinput = strtrim(input,2)
sreference = strtrim(reference,2)
scoords = strtrim(coords,2)
soutput = strtrim(output,2)
sshifts = strtrim(shifts,2)
sboxsize = strtrim(boxsize,2)
sbigbox = strtrim(bigbox,2)
snegative = strtrim(negative,2)
sbackground = strtrim(background,2)
slower = strtrim(lower,2)
supper = strtrim(upper,2)
sniterate = strtrim(niterate,2)
stolerance = strtrim(tolerance,2)
smaxshift = strtrim(maxshift,2)
sshiftimages = strtrim(shiftimages,2)
sinterp_type = strtrim(interp_type,2)
sboundary_type = strtrim(boundary_type,2)
sconstant = strtrim(constant,2)
strimimages = strtrim(trimimages,2)
sverbose = strtrim(verbose,2)

; Write IRAF script
push,cmd,'print("")'   ; first line will be ignored
push,cmd,'cd '+curdir
push,cmd,'immatch'
push,cmd,'imalign.boxsize = '+sboxsize
push,cmd,'imalign.bigbox = '+sbigbox
push,cmd,'imalign.negative = '+snegative
push,cmd,'imalign.background = '+sbackground
push,cmd,'imalign.lower = '+slower
push,cmd,'imalign.upper = '+supper
push,cmd,'imalign.niterate = '+sniterate
push,cmd,'imalign.tolerance = '+stolerance
push,cmd,'imalign.maxshift = '+smaxshift
push,cmd,'imalign.shiftimages = '+sshiftimages
push,cmd,'imalign.interp_type = "'+sinterp_type+'"'
push,cmd,'imalign.boundary_type = "'+sboundary_type+'"'
push,cmd,'imalign.constant = '+sconstant
push,cmd,'imalign.trimimages = '+strimimages
push,cmd,'imalign.verbose = '+sverbose
push,cmd,'imalign("'+sinput+'","'+sreference+'","'+scoords+'","'+soutput+'",shifts="'+sshifts+'")'
;push,cmd,'imalign("'+sinput+'","'+sreference+'","'+scoords+'","'+soutput+'",shifts="'+sshifts+$
;          '",boxsize='+sboxsize+',bigbox='+sbigbox+',negative='+snegative+',background='+$
;          sbackground+',lower='+slower+',upper='+supper+',niterate='+sniterate+$
;          ',tolerance='+stolerance+',maxshift='+smaxshift+',shiftimages='+sshiftimages+$
;          ',interp_type="'+sinterp_type+'",boundary_type="'+sboundary_type+'",constant='+$
;          sconstant+',trimimages='+strimimages+',verbose='+sverbose+')'
push,cmd,'logout'
;cmdfile = MAKETEMP('temp','.cl')
cmdfile = MKTEMP('temp')        ; absolute filename
WRITELINE,cmdfile,cmd


; Goto the IRAF directory
CD,current=curdir
CD,irafdir

; Running IRAF
undefine,out
;SPAWN,'cl < '+curdir+'/'+cmdfile,out,errout
SPAWN,'cl < '+cmdfile,out,errout

; IMALIGN can fail if it doesn't find any matching sources in
; the images.

; The output
lo = first_el(where(stregex(out,'Shifts',/boolean) eq 1))  ; where to start output
lo = lo > 0
printline,out[lo:*]

; Return to original directory
CD,curdir


; Erasing the temporary files
FILE_DELETE,cmdfile,/allow,/quiet


; What were the shifts
; this is not very robust yet.
LOADINPUT,sinput,infiles,count=ninfiles
trans = fltarr(ninfiles,2)
out_trans = out[lo+1:lo+1+ninfiles-1]
transarr = strsplitter(out_trans,' ',/extract)
names = reform(transarr[0,*])
xshift = float(reform(transarr[1,*]))
yshift = float(reform(transarr[3,*]))
trans[*,0] = xshift
trans[*,1] = yshift


; Getting trim section
; [Xstart,Xend,Ystart,Yend]
gdtrim = where(stregex(out,'#Trim_Section',/fold_case,/boolean) eq 1,ngtrim)
tarr = strsplitter(out[gdtrim[0]],' ',/extract)
tsection = reform(tarr[2])  ; i.e. [4:2046,1:2041]
trimsection = [4]
len = strlen(tsection)
tsection = strmid(tsection,0,len-1)
tsection = strmid(tsection,1,len-2)
tsecarr = strsplit(tsection,',',/extract)
trimsection = strsplitter(tsecarr,':',/extract)
trimsection = float(trimsection)
trimsection = (trimsection)(*)


; Because of the trimming this can cause shifts between the
; original and shifted/combined images
head1 = headfits(sreference)
crpix1a = sxpar(head1,'CRPIX1',/silent)
crpix2a = sxpar(head1,'CRPIX2',/silent)
; Image has WCS
if strtrim(crpix1a,2) ne '0' then begin
  base = file_basename(sreference,'.fits')
  head2 = headfits(base+'.shft.fits')
  crpix1b = sxpar(head2,'CRPIX1',/silent)
  crpix2b = sxpar(head2,'CRPIX2',/silent)
  xoff = crpix1a-crpix1b
  yoff = crpix2a-crpix2b
; No WCS, Use LTV1/LTV2
endif else begin
  base = file_basename(sreference,'.fits')
  head = headfits(base+'.shft.fits')
  ; Take negative
  ltv1 = sxpar(head,'LTV1',/silent)
  ltv2 = sxpar(head,'LTV2',/silent)
  ; Take negative
  ; xorig = xshift + xoff
  ; yorig = yshift + yoff
  xoff = float(-ltv1)
  yoff = float(-ltv2)

endelse

print,'Xoff=',strtrim(xoff,2)
print,'Yoff=',strtrim(yoff,2)

if keyword_set(stp) then stop

end
