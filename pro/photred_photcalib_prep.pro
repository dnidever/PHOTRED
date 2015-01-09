pro photred_photcalib_prep,mchfile,apcor,outfile,silent=silent,$
                           error=error,observatory=observatory,$
                           photfile=photfile

;+
;
; PHOTRED_photcalib_PREP
;
; This makes the input file for PHOTCALIB.PRO
; This must be run from the directory that has
; the FITS and MCH files.
; 
; The input list needs to be of MCH files.
; The apcfile is the aperture correction file
;
; INPUTS:
;  mchfile   The MCH file
;  apcor     The aperture correction structure
;  outfile   The output filename
;  =observatory  The observatory name.  Might be needed for
;                  calculating airmass.
;  =photfile The name of the photometry file to use.  Normally
;              this is mchfile with a .ast ending.
; 
; OUTPUTS
;  The input file for PHOTCALIB.PRO
;  =errpr    The error message if there was an error.
;
; USAGE:
;  IDL>photred_photcalib_prep,'ccd1001.mch',apcorstr,'ccd1001.input',error=error
;
; By D.Nidever  March 2008
;-

undefine,error

; Not enough inputs
nmchfile = n_elements(mchfile)
napcor = n_elements(apcor)
noutfile = n_elements(outfile)
if nmchfile eq 0 or napcor eq 0 or noutfile eq 0 then begin
  print,'Syntax - photred_photcalib_prep,mchfiles,apcfile,outfile,error=error'
  error = 'Not enough inputs'
  return
endif

; Aperture arrays
apcnames = apcor.name
apcvalue = apcor.value
apcnames2 = repstr(apcnames,'a.del')    ; base names


; Load the MCH file
;READCOL,list[i],files,format='a',/silent
;READCOL,mchfile,files,format='a',/silent
LOADMCH,mchfile,files,trans

files = strtrim(files,2)
;files = repstr(files,"'")
nfiles = n_elements(files)

; Initializing arrays
filtarr = strarr(nfiles)
amarr = dblarr(nfiles)
exparr = dblarr(nfiles)
apcarr = dblarr(nfiles)

; Loop through the photometry files
for j=0,nfiles-1 do begin

  file = files[j]
  arr = strsplit(file,'.',/extract)
  narr = n_elements(arr)
  if narr ge 2 then filebase = strjoin(arr[0,narr-2]) else filebase = file

  fitsfile = filebase+'.fits'

  ; Reading header of FITS file
  undefine,errmsg
  head = HEADFITS(fitsfile,errmsg=errmsg)
  if (errmsg ne '') then begin
    print,'ERROR in getting header'
    error = errmsg
    goto,BOMB
  endif

  ; Getting the filter name
  filt = PHOTRED_GETFILTER(fitsfile)

  ; Getting the airmass
  am = PHOTRED_GETAIRMASS(fitsfile,obs=observatory,/update,/recalculate)
  ;am = sxpar(head,'AIRMASS')
  if (am lt 0.9) then begin
    print,'ERROR NO AIRMASS'
    error = 'ERROR NO AIRMASS'
    goto,BOMB
  end

  ; Getting the exposure time
  exp = PHOTRED_GETEXPTIME(fitsfile)
  if strtrim(exp,2) eq '-1' then begin
    print,'ERROR NO EXPTIME'
    error = 'ERROR NO EXPTIME'
    goto,BOMB
  endif

  ; Getting the aperture correction
  apcorr = 0.0
  gd = where(apcnames eq filebase,ngd)
  if (ngd gt 0) then apcorr = apcvalue[gd[0]]
  if apcorr lt 0.0 then apcorr=0.0        ; don't want negative correction

  ; Plugging into the arrays
  filtarr[j] = filt
  amarr[j] = am
  exparr[j] = exp
  apcarr[j] = apcorr

end ; looping through phot files

arr = strsplit(mchfile,'.',/extract)
narr = n_elements(arr)
if narr ge 2 then filebase = strjoin(arr[0,narr-2]) else filebase = mchfile


; OUTPUTTING
; RAW filename, Band1 name, Band1 airmass, Band1 exptime, Band1 aperture correction, Band2 ...
tab = '   '
;; DAOPHOT, use RAW file
;if (not keyword_set(allframe)) then begin
;  out = string(filebase,format='(A15)')+'.raw' + tab
;; ALLFRAME, use MAG file
;endif else begin
;  out = string(filebase,format='(A15)')+'.mag' + tab
;endelse

; Photfile input
if n_elements(photfile) gt 0 then begin
  out = string(photfile,format='(A19)') + tab
endif else begin
  out = string(filebase,format='(A15)')+'.ast' + tab
endelse

for j=0,nfiles-1 do begin
  out = out + string(filtarr[j],format='(A6)') + tab
  out = out + string(amarr[j],format='(F7.4)') + tab
  out = out + string(exparr[j],format='(F7.1)') + tab
  out = out + string(apcarr[j],format='(F7.4)') + tab
end


; Opening the output file
OPENW,unit,/get_lun,outfile

; Writing to the output file
PRINTF,unit,out

; Closing the output file
CLOSE,unit
FREE_LUN,unit


BOMB:

if keyword_set(stp) then stop

end
