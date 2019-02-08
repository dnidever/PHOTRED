;+
;
; PHOTRED_GETMETA
;
; Retrive meta-data from the header of a file.
; PHOTRED_GATHERFILEINFO.PRO does a lot of this as well.
;
; INPUTS:
;  file      FITS filename for which to get the meta-data.
;  imager    The imager structure with the chip number separation character.
;  obs       The observatory name.
;
; OUTPUTS:
;  meta    Structure with meta-data.
;  =error  The error message if one occurred.
;
; USAGE:
;  IDL>meta = photred_getmeta('image.fits',imager,obs)
;
; By D. Nidever  Feb 2019
;-

function photred_getmeta,file,imager,obs,error=error

undefine,meta,error

;; Not enough inputs
if n_elements(file) eq 0 or n_elements(imager) eq 0 or n_elements(obs) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - meta = photred_getmeta(file,imager,obs,error=error)'
  return,-1
endif

;; Load the header
head = PHOTRED_READFILE(file,/header,error=error)
if n_elements(error) ne 0 then return,-1

;; Size of the image
naxis1 = sxpar(head,'NAXIS1')
naxis2 = sxpar(head,'NAXIS2')

;; Get chipnum
chipnum = PHOTRED_GETCHIPNUM(file,imager,error=error)
if n_elements(error) ne 0 then return,-1

;; Get the AST structure
EXTAST,head,ast,nparams
;; Get vertices
if nparams ge 1 then HEAD_XYAD,head,[0,naxis1-1,naxis1-1,0],[0,0,naxis2-1,naxis2-1],vra,vdec,/degree

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

;; Get the pixel scale
GETPIXSCALE,file,pixscale,head=head,error=error
if n_elements(error) ne 0 then return,-1

;; Get the exposure time
exptime = PHOTRED_GETEXPTIME(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1

;; Get data, uttime, date-obs
date = PHOTRED_GETDATE(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1
uttime = PHOTRED_GETUTTIME(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1
dateobs = date+'T'+uttime

;; Get MJD
mjd = date2jd(dateobs,/mjd)

;; Get the gain
gain = PHOTRED_GETGAIN(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1

;; Get the rdnoise
rdnoise = PHOTRED_GETRDNOISE(file,head=head,error=error)
if n_elements(error) ne 0 then return,-1

;; Get the airmass
airmass = PHOTRED_GETAIRMASS(file,head=head,obs=obs)
if n_elements(error) ne 0 then return,-1

;; Make the structure
meta = {file:strtrim(file,2),head:head,naxis1:long(naxis1),naxis2:long(naxis2)}
if nparams ge 0 then meta=CREATE_STRUCT(meta,'ast',ast,'vra',double(vra),'vdec',double(vdec))
meta = CREATE_STRUCT(meta,'ra',double(ra),'dec',double(dec),'filter',strtrim(filter,2),'pixscale',float(pixscale),$
                     'exptime',float(exptime),'utdate',strtrim(date,2),'uttime',strtrim(uttime,2),$
                     'dateobs',strtrim(dateobs,2),'mjd',double(mjd),'gain',float(gain),$
                     'rdnoise',float(rdnoise),'airmass',float(airmass))

return,meta

end
