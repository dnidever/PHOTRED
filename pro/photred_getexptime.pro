;+
;
; PHOTRED_GETEXPTIME
;
; This gets exposure time information for a FITS file.
; 
; INPUTS:
;  file      FITS filename
;  =head     Use this header array instead of reading FITS file.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The exposure time is output in seconds.  If there
;  is an error then -1.0 is output.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>exptime = photred_getexptime(file)
;
; By D.Nidever  February 2008
;-

function photred_getexptime,file,head=head,error=error,stp=stp

COMMON photred,setup

undefine,error

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - exptime = photred_getexptime(file,head=head,stp=stp)'
  return,-1.0
endif

;; Can't use input HEAD if multiple fits files input
if nfile gt 1 then undefine,head

; More than one filter name input
if nfile gt 1 then begin
  exptime = fltarr(nfile)
  for i=0,nfile-1 do exptime[i] = photred_getexptime(file[i])
  return,exptime
endif

;; No header input, read from fits file
if n_elements(head) eq 0 then begin
  ;; Check that the file exists
  test = file_test(file)
  if test eq 0 then begin
    error = file+' NOT FOUND'
    print,error
    return,-1.0
  endif
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then head=PHOTRED_READFILE(file,exten=1,/header) else $
    head = PHOTRED_READFILE(file,/header)
endif

exptime = SXPAR(head,'EXPTIME',/silent,count=nexptime)

; No EXPTIME
if nexptime eq 0 then exptime=-1.0
if exptime lt 0.0 then begin
  error = 'NO EXPOSURE TIME'
  print,error
  return,-1.0
endif

if keyword_set(stp) then stop

return,exptime

end
