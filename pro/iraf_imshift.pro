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

; Input strings
sinput = strtrim(input,2)
soutput = strtrim(output,2)
if n_elements(shifts_file) gt 0 then sshifts_file = strtrim(shifts_file,2) else sshifts_file=''
if n_elements(xshift) gt 0 then sxshift=strtrim(xshift,2)
if n_elements(yshift) gt 0 then syshift=strtrim(yshift,2)
sinterp_type = strtrim(interp_type,2)
sboundary_type = strtrim(boundary_type,2)
sconstant = strtrim(constant,2)

; Write IRAF script
push,cmd,'cd '+curdir
push,cmd,'immatch'
push,cmd,'imshift.shifts_file = "'+sshifts_file+'"'
push,cmd,'imshift.interp_type = "'+sinterp_type+'"'
push,cmd,'imshift.boundary_type = "'+sboundary_type+'"'
push,cmd,'imshift.constant = "'+sconstant+'"'
if n_elements(xshift) gt 0 and n_elements(yshift) gt 0 then begin
  push,cmd,'imshift("'+sinput+'","'+soutput+'",xshift='+sxshift+',yshift='+syshift+')'
endif else begin
  push,cmd,'imshift("'+sinput+'","'+soutput+'")'
endelse
push,cmd,'logout'
cmdfile = MKTEMP('temp')        ; absolute filename
WRITELINE,cmdfile,cmd


; Goto the IRAF directory
CD,current=curdir
CD,irafdir

; Running IRAF
undefine,out
SPAWN,'cl < '+cmdfile,out,errout

; Return to original directory
CD,curdir

; Erasing the temporary files
FILE_DELETE,cmdfile,/allow,/quiet

if keyword_set(stp) then stop

end
