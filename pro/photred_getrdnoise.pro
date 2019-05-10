;+
;
; PHOTRED_GETRDNOISE
;
; This gets the READNOISE information from a FITS files
; 
; INPUTS:
;  file      FITS filename
;  =head     Use this header array instead of reading FITS file. 
;  =keyword  The actual header keyword used.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The READNOISE in electrons/read
;  If the readnoise is not found in the header then -1.0
;    is returned.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>rdnoise = photred_getrdnoise(file,keyword=keyword)
;
; By D.Nidever  May 2008
;-

function photred_getrdnoise,file,head=head,keyword=keyword,error=error,stp=stp

COMMON photred,setup

undefine,rdnoise,keyword,error

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - rdnoise = photred_getrdnoise(file,keyword=keyword,stp=stp)'
  return,-1.0
endif

;; Can't use input HEAD if multiple fits files input
if nfile gt 1 then undefine,head

; More than one name input
if nfile gt 1 then begin
  rdnoise = fltarr(nfile)
  keyword = strarr(keyword)
  for i=0,nfile-1 do begin
    rdnoise[i] = photred_getrdnoise(file[i],keyword=ikeyword)
    keyword[i] = ikeyword
  endfor
  return,rdnoise
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

rdnoise = SXPAR(head,'RDNOISE',count=nrdnoise,/silent)
readnois = SXPAR(head,'READNOIS',count=nreadnois,/silent)  ; Swope
enoise = SXPAR(head,'ENOISE',count=nenoise,/silent)        ; IMACS

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
  error = 'NO READNOISE'
  print,error
  readnoise = -1.0
endif

if keyword_set(stp) then stop

return,readnoise

end
