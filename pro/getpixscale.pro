;+
;
; GETPIXSCALE
;
; Get the pixel scale for an image.
;
; INPUTS:
;  file    FITS filename
;  =head   The image header for which to determine the pixel scale.
;  /stp    Stop at the end of the program.
;
; OUTPUTS:
;  scale   The pixel scale of the image in arcsec/pix.
;  =error  The error if one occurred.
;
; USAGE:
;  IDL>getpixscale,'ccd1001.fits',scale
;
; BY D. Nidever   February 2008
;-

pro getpixscale,file,scale,head=head,error=error,stp=stp

undefine,error
scale = -1         ; bad until proven good

; Not enough inputs
nfile = n_elements(file)
if nfile eq 0 and n_elements(head) eq 0 then begin
  error = 'Not enougn inputs'
  print,'Syntax - getpixscale,file,scale,head=head'
  return
endif

;; No header input, read from fits file
fpack = 0
if n_elements(head) eq 0 then begin
  ;; Check that the file exists
  test = file_test(file)
  if test eq 0 and n_elements(head) eq 0 then begin
    error = file+' NOT FOUND'
    print,error
    return
  endif

  ; Fpack or regular fits
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
    fpack = 1
    exten = 1
  endif else begin
    fpack = 0
    exten = 0
  endelse

  ; Read the header
  if n_elements(head) eq 0 then head = PHOTRED_READFILE(file,exten=exten,/header)

  ; Fix NAXIS1/2 in header
  if fpack eq 1 then begin
    sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
    sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
  endif
endif

; Does the image have a SCALE parameter
hscale = sxpar(head,'SCALE',count=nhscale,/silent)
if nhscale ne 0 then scale=hscale
; Does the image have a PIXSCALE parameter
if scale eq -1 then begin
  pixscale = sxpar(head,'PIXSCALE',count=npixscale,/silent)
  if npixscale ne 0 then scale=pixscale
endif
; Does the image have a PIXSCALE1 parameter
if scale eq -1 then begin
  pixscale1 = sxpar(head,'PIXSCALE1',count=npixscale1,/silent)
  if npixscale1 ne 0 then scale=pixscale1
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
  error = 'WARNING! COULD NOT DETERMINE THE PIXEL SCALE'
  print,error
endif

;; Are there any other "scale"-like values in the header
;gd = where(stregex(head,'scale',/boolean,/fold_case) eq 1,ngd)

;stop

if keyword_set(stp) then stop

end
