OB;+
;
; PHOTRED_GETCHIPNUM
;
; This gets the chip (amplifier) number for a given file.
; 
; INPUTS:
;  file      The FITS filename
;  imager    The imager structure with the chip number separation character.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The chip (amplifier) number for the file given the imager naming convention.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>chipnum = photred_getchipnum(file,imager)
;
; By D.Nidever  January 2016
;-

function photred_getchipnum,file,imager,stp=stp,error=error
  
COMMON photred,setup

undefine,error

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 or n_elements(imager) eq 0 then begin
  print,'Syntax - chipnum = photred_getchipnum(file,imager,stp=stp,error=error)'
  error = 'Not enough inputs' 
  return,-1
endif

tilesep = '+'

; Check the imager structure
if size(imager,/type) ne 8 then begin
  error = 'IMAGER must be a structure'
  print,error
  return,-1
endif
if tag_exist(imager,'separator') eq 0 or tag_exist(imager,'namps') eq 0 then begin
  error = 'IMAGER must have SEPARATOR and NAMPS tags.'
  print,error
  return,-1
endif

; Parse the filename
if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then base=file_basename(file,'.fits.fz') else $
  base = file_basename(file,'.fits')   
shfield = first_el(strsplit(base,'-',/extract))

; --- Multi AMP imager ---
if imager.namps gt 1 then begin
  if strpos(base,imager.separator) eq -1 then return,-1  ;; no amp/chip portion
  tmp = first_el(strsplit(base,'-',/extract),/last)
  expnum = first_el(strsplit(tmp,imager.separator,/extract))
  strchip = first_el(strsplit(tmp,imager.separator,/extract),/last)
  ;; Deal with TILE suffix
  if stregex(base,'\'+tilesep+'T',/boolean) eq 1 then strchip=(strsplit(strchip,tilesep+'T',/extract))[0]
  ; Check that it's an actual number/integer
  isnum = valid_num(strchip,chipnum,/integer)
  if isnum eq 0 then begin
    error = strchip+' IS NOT an INTEGER'
    print,error
    return,-1
  endif

; --- Single AMP imager ---
endif else begin
  chip = 1L
  expnum = first_el(strsplit(base,'-',/extract),/last)
  ;; Deal with TILE suffix
  if stregex(base,'\'+tilesep+'T',/boolean) eq 1 then expnum=(strsplit(expnum,tilesep+'T',/extract))[0]
endelse

if keyword_set(stp) then stop

return,chipnum

end
