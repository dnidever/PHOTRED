;+
;
; PHOTRED_MKSEXCONFIG
;
; This program creates a SExtractor configuration file
; for an inputs FITS file.
;
; INPUTS:
;  file        Filename of the FITS image
;  configfile  Filename for the SExtractor configuration file.
;                This is optional.  If not input then the
;                name will be FITSFILEBASE.sex
;  catfile     The output catalog name.  If not input then
;                FITSFILEBASE.cat will be used
;  /stp        Stop at the end of the program
;
; OUTPUTS:
;  The SExtractor config file will be written to CONFIGFILE.
;
; USAGE:
;  IDL>photred_mksexconfig,'F1-89230023_05.fits','F1-89230023_05.sex'
;
; By D. Nidever    January 2019
;-

pro photred_mksexconfig,file,configfile,catfile,stp=stp

undefine,error

nfile = n_elements(file)
if nfile eq 0 then begin
  print,'Syntax - photred_mksexconfig,fitsfile,configfile,catfile,stp=stp'
  return
endif

base = FILE_BASENAME(file,'.fits')

; Make sure file exists
if file_test(file) eq 0 then begin
  error = file+' NOT FOUND'
  print,error
  return
endif
if file_test(base+'.opt') eq 0 then begin
  error = bsae+'.opt NOT FOUND'
  print,error
  return
endif

; Read the .opt file
READLINE,base+'.opt',optlines
optarr = strsplitter(optlines,' ',/extract)
g = where(stregex(optlines,'HI =',/boolean) eq 1,ng)
satlevel = optarr[2,g[0]]
g = where(stregex(optlines,'GA =',/boolean) eq 1,ng)
gain = optarr[2,g[0]]
g = where(stregex(optlines,'FW =',/boolean) eq 1,ng)
fwhm = optarr[2,g[0]]

; Get pixel scale
GETPIXSCALE,file,scale
if scale eq 99.9 then scale=0.5   ; default

; Make customized SEXTRACTOR file
;--------------------------------
if n_elements(configfile) eq 0 then configfile=base+'.sex'
if n_elements(catfile) eq 0 then catfile=base+'.cat'
FILE_COPY,'default.sex',configfile,/overwrite,/allow
READLINE,configfile,sexlines
sexlines2 = sexlines
; CATALOG_NAME
g = where(stregex(sexlines2,'^CATALOG_NAME',/boolean) eq 1,ng)
;catfile = base+'.cat'
sexlines2[g[0]] = 'CATALOG_NAME    '+catfile+' # name of the output catalog'
; SATUR_LEVEL
g = where(stregex(sexlines2,'^SATUR_LEVEL',/boolean) eq 1,ng)
sexlines2[g[0]] = 'SATUR_LEVEL     '+satlevel+'         # level (in ADUs) at which arises saturation'
; GAIN
g = where(stregex(sexlines2,'^GAIN',/boolean) eq 1,ng)
sexlines2[g[0]] = 'GAIN            '+gain+'             # detector gain in e-/ADU.'
; PIXEL_SCALE
g = where(stregex(sexlines2,'^PIXEL_SCALE',/boolean) eq 1,ng)
sexlines2[g[0]] = 'PIXEL_SCALE     '+strtrim(scale,2)+'             # size of pixel in arcsec (0=use FITS WCS info).'
; SEEING_FWHM
g = where(stregex(sexlines2,'^SEEING_FWHM',/boolean) eq 1,ng)
fwhmas = float(fwhm)*float(scale)
sexlines2[g[0]] = 'SEEING_FWHM     '+strtrim(fwhmas,2)+'            # stellar FWHM in arcsec'
; MASK/WEIGHT file
if n_elements(maskfile) ne 0 then begin
  PUSH,sexlines2,'WEIGHT_IMAGE  '+maskfile
  PUSH,sexlines2,'WEIGHT_TYPE   MAP_WEIGHT'
endif
; Write the file
WRITELINE,configfile,sexlines2

if keyword_set(stp) then stop

end
