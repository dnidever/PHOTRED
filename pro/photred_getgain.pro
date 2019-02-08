;+
;
; PHOTRED_GETGAIN
;
; This gets the GAIN information from a FITS files
; 
; INPUTS:
;  file      FITS filename
;  =head     Use this header array instead of reading FITS file.
;  =keyword  The actual header keyword used.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The gain in e/ADU
;  If the gain is not found in the header then -1.0
;    is returned.
;  =error    The error message if one occurred.  
;
; USAGE:
;  IDL>gain = photred_getgain(file,keyword=keyword)
;
; By D.Nidever  May 2008
;-

function photred_getgain,file,head=head,keyword=keyword,error=error,stp=stp

COMMON photred,setup

undefine,gain,keyword,error

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - gain = photred_getgain(file,keyword=keywordstp=stp)'
  return,-1.0
endif

;; Can't use input HEAD if multiple fits files input
if nfile gt 1 then undefine,head

; More than one name input
if nfile gt 1 then begin
  gain = fltarr(nfile)
  keyword = strarr(nfile)
  for i=0,nfile-1 do begin
    gain[i] = photred_getgain(file[i],keyword=ikeyword)
    keyword[i] = ikeyword
  endfor
  return,gain
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

gain = SXPAR(head,'GAIN',count=ngain,/silent)
egain = SXPAR(head,'EGAIN',count=negain,/silent)          ; imacs

; Use GAIN
if ngain gt 0 then begin
  keyword = 'GAIN'
endif
; Use EGAIN
if ngain eq 0 and negain gt 0 then begin
  gain = egain
  keyword = 'EGAIN'
endif


; No GAIN
if (ngain eq 0 and negain eq 0) then begin
  error = 'NO GAIN'
  print,error
  gain = -1.0
endif

if keyword_set(stp) then stop

return,gain

end
