function photred_getrdnoise,file,keyword=keyword,stp=stp


;+
;
; PHOTRED_GETRDNOISE
;
; This gets the READNOISE information from a FITS files
; 
; INPUTS:
;  file      FITS filename
;  =keyword  The actual header keyword used.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The READNOISE in electrons/read
;  If the readnoise is not found in the header then -1.0
;    is returned.
;
; USAGE:
;  IDL>rdnoise = photred_getrdnoise(file,keyword=keyword)
;
; By D.Nidever  May 2008
;-

COMMON photred,setup

undefine,rdnoise,keyword

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  print,'Syntax - rdnoise = photred_getrdnoise(file,keyword=keyword,stp=stp)'
  return,-1.0
endif

; More than one name input
if nfile gt 1 then begin
  rdnoise = fltarr(nfile)
  keyword = strarr(keyword)
  for i=0,nfile-1 do begin
    rdnoise[i] = photred_getrdnoise(file[i],keyword=ikeyword)
    keyword[i] = ikeyword
  end
  return,rdnoise
endif

test = file_test(file)
if test eq 0 then begin
  print,file,' NOT FOUND'
  return,-1.0
endif

head = HEADFITS(file)
rdnoise = SXPAR(head,'RDNOISE',count=nrdnoise)
readnois = SXPAR(head,'READNOIS',count=nreadnois)  ; Swope
enoise = SXPAR(head,'ENOISE',count=nenoise)        ; IMACS

; Use RDNOISE
if nrdnoise gt 0 then begin
  readnoise = rdnoise
  keyword = 'RDNOISE'
endif
; Use READNOIS
if nrdnoise eq 0 and nreadnois gt 0 then begin
  readnoise = readnois
  keyword = 'READNOIS'
endif
; Use ENOISE
if nrdnoise eq 0 and nreadnois eq 0 and nenoise gt 0 then begin
  readnoise = enoise
  keyword = 'ENOISE'
endif

; No Read Noise
if (nrdnoise eq 0 and nreadnois eq 0 and nenoise eq 0) then begin
  print,'NO READNOISE'
  readnoise = -1.0
endif

if keyword_set(stp) then stop

return,readnoise

end
