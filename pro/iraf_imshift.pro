;+
;
; IRAF_IMSHIFT
; 
; PURPOSE:
;  This runs IRAF's IMSHIFT that shifts images
;
; INPUTS:
;    
;    input
;        List of images to be transformed.
;    
;    output
;        List of output images.
;    
;    xshift, yshift
;        Fractional  pixel shift in x and y such that xout = xin + xshift
;        and yout = yin + yshift.
;    
;    shifts_file = ""
;        The name of the text file containing the shifts for  each  input
;        image.  If  no file name is supplied each input image is shifted
;        by xshift and yshift. Shifts are listed in the text file, 1  set
;        of  shifts  per image, with the x and y shift in columns 1 and 2
;        respectively. The number of shifts in the file  must  equal  the
;        number of input images.
;    
;    interp_type = "linear"
;        The  interpolant  type use to computed the output shifted image.
;        The choices are the following:
;        
;        nearest
;            nearest neighbor.
;        
;        linear
;            bilinear interpolation in x and y.
;        
;        poly3
;            third order interior polynomial in x and y.
;        
;        poly5
;            fifth order interior polynomial in x and y.
;        
;        spline3
;            bicubic spline.
;        
;        sinc
;            2D  sinc  interpolation.  Users   can   specify   the   sinc 
;            interpolant   width  by  appending  a  width  value  to  the 
;            interpolant string, e.g. sinc51 specifies a 51 by  51  pixel
;            wide  sinc  interpolant.  The  sinc  width input by the user
;            will be rounded up to the nearest odd  number.  The  default
;            sinc width is 31 by 31.
;        
;        drizzle
;            2D  drizzle  resampling. Users can specify the drizzle pixel
;            fractions in x and y by appending  values  between  0.0  and
;            1.0  in  square  brackets  to  the  interpolant string, e.g.
;            drizzle[0.5]. The default value is 1.0.  The  value  0.0  is
;            increased to 0.001. Drizzle resampling with a pixel fraction
;            of 1.0 in x and y is identical to bilinear interpolation.
;    
;    boundary_type = "nearest"
;        The choices are:
;        
;        nearest
;            Use the value of the nearest boundary pixel.
;        
;        constant
;            Use a constant value.
;        
;        reflect
;            Generate value by reflecting about the boundary.
;        
;        wrap
;            Generate a value by wrapping around to the opposite side  of
;            the image.
;
;   =constant    The constant value to use for boundary_type="constant".
;   =irafdir
;        The IRAF home diretory.
;
;
; OUTPUTS:
;  Shifted images with the names in the "output" file.
;
; USAGE:
;  IDL>iraf_imshift,'@inlist','@outlist',shifts_file='file.shift'
;
;
;                                    I R A F  
;                     Image Reduction and Analysis Facility
; PACKAGE = immatch
;    TASK = imshift
; 
; input   =                       Input images to be fit
; output  =                       Output images
; xshift  =                       Fractional pixel shift in x
; yshift  =                       Fractional pixel shift in y
; (shifts_=                     ) Text file containing shifts for each image
; (interp_=               linear) Interpolant (nearest,linear,poly3,poly5,spline3,sinc,drizzle)
; (boundar=              nearest) Boundary (constant,nearest,reflect,wrap)
; (constan=                   0.) Constant for boundary extension
; (mode   =                   ql)
;
; By D. Nidever    February 2013
;-

pro iraf_imshift,input,output,xshift,yshift,shifts_file=shifts_file,interp_type=interp_type,$
                boundary_type=boundary_type,constant=constant,verbose=verbose,stp=stp,error=error,$
                irafdir=irafdir

undefine,error

ninput = n_elements(input)
noutput = n_elements(output)
nxshift = n_elements(xshift)
nyshift = n_elements(yshift)
nshifts_file = n_elements(shifts_file)

; Not enough inputs
if ninput eq 0 or noutput eq 0 or (nshifts_file eq 0 and (nxshift eq 0 or nyshift eq 0)) then begin
  error = 'Not enough inputs'
  print,'Syntax -  iraf_shift,input,output,xshift,yshift,shifts_file=shifts_file,interp_type=interp_type,'
  print,'                boundary_type=boundary_type,constant=constant,verbose=verbose,stp=stp,error=error,'
  print,'                irafdir=irafdir'
  return
endif

; xshift/yshift must be scalar
if nxshift gt 1 or nyshift gt 1 then begin
  error = 'XSHIFT/YSHIFT must be scalar'
  print,error
  return
endif

; Important directories
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
;  These are the defaults from IMALIGN
if n_elements(interp_type) eq 0 then interp_type='spline3'
if n_elements(boundary_type) eq 0 then boundary_type='constant'
if n_elements(constant) eq 0 then constant=0.

; Multiple files input
if n_elements(shifts_file) gt 0 then begin
  if strmid(input,0,1) eq '@' then readcol,strmid(input,1),inpnames,format='A',/silent
  if strmid(output,0,1) eq '@' then readcol,strmid(output,1),outnames,format='A',/silent
  readcol,shifts_file,xshift,yshift,format='F,F',/silent

; One file
endif else begin
  inpnames = input
  outnames = output
endelse
nfiles = n_elements(inpnames)

; Input strings
sinterp_type = strtrim(interp_type,2)
sboundary_type = strtrim(boundary_type,2)
sconstant = strtrim(constant,2)

; Write IRAF script
push,cmd,'print("")'   ; first line will be ignored
push,cmd,'cd '+curdir
push,cmd,'immatch'
push,cmd,'imshift.input=""'
push,cmd,'imshift.output=""'
push,cmd,'imshift.shifts_file = ""'  ; don't use shifts_file use x/yshifts
push,cmd,'imshift.interp_type = "'+sinterp_type+'"'
push,cmd,'imshift.boundary_type = "'+sboundary_type+'"'
push,cmd,'imshift.constant = '+sconstant

; Loop through files
for i=0,nfiles-1 do begin
  sinpname1 = strtrim(inpnames[i],2)
  soutname1 = strtrim(outnames[i],2)
  ; if shift is very small round to zero, otherwise can cause problems
  if abs(xshift[i]) lt 1e-4 then xshift1=0.0 else xshift1=xshift[i]
  if abs(yshift[i]) lt 1e-4 then yshift1=0.0 else yshift1=yshift[i]
  sxshift1 = strtrim(xshift1,2)
  syshift1 = strtrim(yshift1,2)
  push,cmd,'imshift("'+sinpname1+'","'+soutname1+'",xshift='+sxshift1+',yshift='+syshift1+')'
  push,cmd,'flprcache'
endfor

push,cmd,'logout'
cmdfile = MKTEMP('temp')        ; absolute filename
WRITELINE,cmdfile,cmd


; Goto the IRAF directory
CD,current=curdir
CD,irafdir

; Running IRAF
undefine,out
SPAWN,'cl < '+cmdfile,out,errout

if keyword_set(verbose) then print,out,errout

; Return to original directory
CD,curdir

; Erasing the temporary files
FILE_DELETE,cmdfile,/allow,/quiet

if keyword_set(stp) then stop

end
