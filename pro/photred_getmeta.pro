;+
;
; PHOTRED_GETMETA
;
; Retrive meta-data from the header of a file
;
; INPUTS:
;  file    FITS filename for which to get the meta-data.
;
; OUTPUTS:
;  meta    Structure with meta-data.
;  =error  The error message if one occurred.
;
; USAGE:
;  IDL>meta = photred_getmeta('image.fits')
;
; By D. Nidever  Feb 2019
;-

function photred_getmeta,file,error=error

undefine,meta,error

;; Not enough inputs
if n_elements(file) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - meta = photred_getmeta(file,error=error)'
  return,-1
endif

;; Get the header
head = PHOTRED_READFILE(file,/header,error=error)
if n_elements(error) ne 0 then return,-1

;; Size of the image
naxis1 = sxpar(head,'NAXIS1')
naxis2 = sxpar(head,'NAXIS2')

;; Get chipnum
chipnum = PHOTRED_GETCHIPNUM(file,error=error)
if n_elements(error) ne 0 then return,-1

;; Get the AST structure
EXAST,head,ast,nparams

;; Get RA and DEC
sra = sxpar(head,'RA',count=nra,/silent)
if nra gt 0 then begin
  if strpos(sra,':') eq -1 then ra=double(sra) else ra=double(sexig2ten(sra))
endif else begin
  if nparams lt 1 then begin  ; no WCS information
    error = 'No coordinate/WCS information in header'
    return,-1
  endif
  HEAD_XYAD,head,naxis1/2,naxis2/2,ra,dec,/degree,error=error
  if n_elements(error) gt 0 then return,-1
endelse
sdec = sxpar(head,'DEC',count=ndec,/silent)
if ndec gt 0 then begin
  if strpos(sdec,':') eq -1 then dec=double(sdec) else dec=double(sexig2ten(sdec))
endif else begin
  if nparams lt 1 then begin  ; no WCS information
    error = 'No coordinate/WCS information in header'
    return,-1
  endif
  HEAD_XYAD,head,naxis1/2,naxis2/2,ra,dec,/degree,error=error
  if n_elements(error) gt 0 then return,-1
endelse

;; Get the filter
filter = PHOTRED_GETFILTER(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1

;; Get the exposure time
exptime = PHOTRED_GETEXPTIME(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1

;; Get data, uttime, date-obs
date = PHOTRED_GETDATE(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1
uttime = PHOTRED_UTTIME(file,head=head,,error=error)
if n_elements(error) ne 0 then return,-1
dateobs = date+'T'+uttime

;; Get MJD
mjd = date2jd(dateobs,/mjd)

;; Get the gain
gain = PHOTRED_GETGAIN(file,head=head,,error=error)
if n_elements(error) ne 0 then return,-1

;; Get the rdnoise
rdnoise = PHOTRED_GETRDNOISE(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1

;; Get the airmass
airmass = PHOTRED_GETAIRMASS(file,head=head,obs=observatory)
if n_elements(error) ne 0 then return,-1

;; Make the structure

stop

return,meta

end
