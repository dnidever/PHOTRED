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
;  =imager   The imaging structure.  Needed to find the chip number.
;  =observatory  The observatory name.  Might be needed for
;                  calculating airmass.
;  =photfile The name of the photometry file to use.  Normally
;              this is mchfile with a .ast ending.
; 
; OUTPUTS
;  The input file for PHOTCALIB.PRO
;  =error    The error message if there was an error.
;
; USAGE:
;  IDL>photred_photcalib_prep,'ccd1001.mch',apcorstr,'ccd1001.input',error=error
;
; By D.Nidever  March 2008
;-

pro photred_photcalib_prep,mchfile,apcor,outfile,silent=silent,$
                           error=error,imager=imager,observatory=observatory,$
                           photfile=photfile

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

tilesep = '+'

; Aperture arrays
apcnames = apcor.name
apcvalue = apcor.value
apcnames2 = repstr(apcnames,'a.del','')    ; base names


; Load the MCH file
;READCOL,list[i],files,format='a',/silent
;READCOL,mchfile,files,format='a',/silent
LOADMCH,mchfile,files,trans

files = strtrim(files,2)
;files = repstr(files,"'",'')
nfiles = n_elements(files)

; Initializing structure
info = replicate({file:'',mjd:0L,chipnum:0L,filter:'',am:0.0d0,exptime:0.0d0,apcorr:0.0d0},nfiles)

; Loop through the photometry files
for j=0,nfiles-1 do begin

  file = files[j]
  arr = strsplit(file,'.',/extract)
  narr = n_elements(arr)
  if narr ge 2 then filebase = strjoin(arr[0,narr-2]) else filebase = file

  fitsfile = filebase+'.fits'
  if file_test(fitsfile) eq 0 then fitsfile=filebase+'.fits.fz'

  ; Get the chip/amp number
  chipnum = PHOTRED_GETCHIPNUM(fitsfile,imager,error=errorchipnum)
  if chipnum le 0 or n_elements(errorchipnum) gt 0 then begin
    error = 'ERROR NO CHIPNUM'
    print,error
    goto,BOMB
  endif
  
  ; Reading header of FITS file
  undefine,errmsg
  if strmid(fitsfile,6,7,/reverse_offset) eq 'fits.fz' then $
    head = PHOTRED_READFILE(fitsfile,exten=1,error=error,/header) else $
    head = PHOTRED_READFILE(fitsfile,error=error,/header)
  if n_elements(error) gt 0 then begin
    print,'ERROR in getting header'
    error = errmsg
    goto,BOMB
  endif

  ; Getting the filter name
  filt = PHOTRED_GETFILTER(fitsfile)

  ; Getting MJD
  mjd = PHOTRED_GETMJD(fitsfile,observatory,error=errormjd)
  if mjd lt 0 or n_elements(errormjd) gt 0 then begin
    error = 'ERROR NO MJD'
    print,error
    goto,BOMB
  endif
  
  ; Getting the airmass
  ;am = PHOTRED_GETAIRMASS(fitsfile,obs=observatory,/update,/recalculate)
  am = PHOTRED_GETAIRMASS(fitsfile,obs=observatory,/recalculate)
  ;am = sxpar(head,'AIRMASS')
  if (am lt 0.9) then begin
    error = 'ERROR NO AIRMASS'
    print,error
    goto,BOMB
  endif

  ; Getting the exposure time
  exp = PHOTRED_GETEXPTIME(fitsfile)
  if strtrim(exp,2) eq '-1' then begin
    error = 'ERROR NO EXPTIME'
    print,error
    goto,BOMB
  endif

  ; Getting the aperture correction
  apcorr = 0.0
  gd = where(apcnames eq filebase,ngd)
  if ngd eq 0 and stregex(filebase,'\'+tilesep+'T',/boolean) eq 1 then begin
    filebase1 = (strsplit(filebase,tilesep+'T',/extract))[0]
    gd = where(apcnames eq filebase1,ngd)
  endif
  if (ngd gt 0) then apcorr = apcvalue[gd[0]]
  if apcorr lt 0.0 then apcorr=0.0        ; don't want negative correction

  ; Plugging into the structure
  info[j].file = fitsfile
  info[j].mjd = mjd
  info[j].chipnum = chipnum
  info[j].filter = filt
  info[j].am = am
  info[j].exptime = exp
  info[j].apcorr = apcorr

endfor ; looping through phot files

arr = strsplit(mchfile,'.',/extract)
narr = n_elements(arr)
if narr ge 2 then filebase = strjoin(arr[0,narr-2]) else filebase = mchfile


; OUTPUTTING
;------------
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
  len = strlen(photfile)
  out = string(photfile,format='(A'+strtrim(len,2)+')') + tab
  ;out = string(photfile,format='(A19)') + tab
endif else begin
  len = strlen(filebase)
  out = string(filebase,format='(A'+strtrim(len,2)+')')+'.ast' + tab
  ;out = string(filebase,format='(A15)')+'.ast' + tab
endelse

; Keep all exposures/bands on ONE line for now
for j=0,nfiles-1 do begin
  out = out + string(info[j].mjd,format='(I8)') + tab
  out = out + string(info[j].chipnum,format='(I5)') + tab
  out = out + string(info[j].filter,format='(A6)') + tab
  out = out + string(info[j].am,format='(F7.4)') + tab
  out = out + string(info[j].exptime,format='(F7.1)') + tab
  out = out + string(info[j].apcorr,format='(F7.4)') + tab
endfor


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
