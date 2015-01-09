function photred_getgain,file,keyword=keyword,stp=stp


;+
;
; PHOTRED_GETGAIN
;
; This gets the GAIN information from a FITS files
; 
; INPUTS:
;  file      FITS filename
;  =keyword  The actual header keyword used.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The gain in e/ADU
;  If the gain is not found in the header then -1.0
;    is returned.
;
; USAGE:
;  IDL>gain = photred_getgain(file,keyword=keyword)
;
; By D.Nidever  May 2008
;-

COMMON photred,setup

undefine,gain,keyword

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  print,'Syntax - gain = photred_getgain(file,keyword=keywordstp=stp)'
  return,-1.0
endif

; More than one name input
if nfile gt 1 then begin
  gain = fltarr(nfile)
  keyword = strarr(nfile)
  for i=0,nfile-1 do begin
    gain[i] = photred_getgain(file[i],keyword=ikeyword)
    keyword[i] = ikeyword
  end
  return,gain
endif

test = file_test(file)
if test eq 0 then begin
  print,file,' NOT FOUND'
  return,-1.0
endif

head = HEADFITS(file)
gain = SXPAR(head,'GAIN',count=ngain)
egain = SXPAR(head,'EGAIN',count=negain)          ; imacs

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
  print,'NO GAIN'
  gain = -1.0
endif

if keyword_set(stp) then stop

return,gain

end
