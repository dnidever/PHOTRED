pro getpixscale,file,scale,stp=stp

;+
;
; GETPIXSCALE
;
; Get the pixel scale for an image
;
; INPUTS:
;  file   FITS filename
;  /stp   Stop at the end of the program.
;
; OUTPUTS:
;  scale  The pixel scale of the image in arcsec/pix.
;
; USAGE:
;  IDL>getpixscale,'ccd1001.fits',scale
;
; BY D. Nidever   February 2008
;-

scale = -1         ; bad until proven good

; Not enough inputs
nfile = n_elements(file)
if nfile eq 0 then begin
  print,'Syntax - getpixscale,file,scale'
  return
endif

test = file_test(file)
if test eq 0 then begin
  print,file,' NOT FOUND'
  return
endif

; Read the header
head = headfits(file)

; Does the image have a SCALE parameter
hscale = sxpar(head,'SCALE',/silent)
if strtrim(hscale,2) ne '0' then scale=hscale
; Does the image have a PIXSCALE parameter
if scale eq -1 then begin
  pixscale = sxpar(head,'PIXSCALE',/silent)
  if strtrim(pixscale,2) ne '0' then scale=pixscale
endif
; Does the image have a PIXSCALE1 parameter
if scale eq -1 then begin
  pixscale1 = sxpar(head,'PIXSCALE1',/silent)
  if strtrim(pixscale1,2) ne '0' then scale=pixscale1
endif

; Try the WCS
if scale eq -1 then begin
  EXTAST,head,astr
  nastr = n_elements(astr)

  ; The image has a WCS
  if nastr gt 0 then begin

    ; Get the coordinates for two positions
    ; separated by 1 pixel
    ;xy2ad,0.0,0.0,astr,ra1,dec1
    ;xy2ad,1.0,0.0,astr,ra2,dec2
    head_xyad,head,0.0,0.0,ra1,dec1,/degree
    head_xyad,head,1.0,0.0,ra2,dec2,/degree
    dist = sphdist(ra1,dec1,ra2,dec2,/deg)*3600.
    scale = dist

    if scale eq 0.0 then scale=99.9

    ;; astr.CD are how xi/eta vary with x/y
    ;scale = max(abs(astr.cd))*3600.0

  endif

endif

; Couldn't determine the pixel scale
if scale eq -1 then begin
  print,'WARNING! COULD NOT DETERMINE THE PIXEL SCALE'
endif

;; Are there any other "scale"-like values in the header
;gd = where(stregex(head,'scale',/boolean,/fold_case) eq 1,ngd)

if keyword_set(stp) then stop

end
