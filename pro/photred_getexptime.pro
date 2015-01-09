function photred_getexptime,file,stp=stp


;+
;
; PHOTRED_GETEXPTIME
;
; This gets exposure time information for a FITS file.
; 
; INPUTS:
;  file      FITS filename
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The exposure time is output in seconds.  If there
;  is an error then -1.0 is output.
;
; USAGE:
;  IDL>exptime = photred_getexptime(file)
;
; By D.Nidever  February 2008
;-

COMMON photred,setup

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  print,'Syntax - exptime = photred_getexptime(file,stp=stp)'
  return,-1.0
endif

; More than one filter name input
if nfile gt 1 then begin
  exptime = fltarr(nfile)
  for i=0,nfile-1 do exptime[i] = photred_getexptime(file[i])
  return,exptime
endif

test = file_test(file)
if test eq 0 then begin
  print,file,' NOT FOUND'
  return,-1.0
endif

head = HEADFITS(file)
exptime = SXPAR(head,'EXPTIME')

; No EXPTIME
if strtrim(exptime,2) eq '0' then exptime=-1.0
if exptime lt 0.0 then begin
  print,'NO EXPOSURE TIME'
  return,-1.0
endif

if keyword_set(stp) then stop

return,exptime

end
