;+
;
; PHOTRED_GETAIRMASS
;
; This gets the AIRMASS for a FITS file.
; 
; INPUTS:
;  file      FITS filename
;  =head     Use this header array instead of reading FITS file. 
;  =obs      Observatory name.  This is needed if the airmas
;              needs to be recalculated.
;  /update   Update the FITS header with the airmass.
;  /recalculate  Recalculate the AIRMASS for the center of the FITS file.
;  /silent   Don't print anything to the screen.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The AIRMASS is output.  If there is an error
;  then -1.0 is output.
;  =error    The error message if an error occurred. 
;
; USAGE:
;  IDL>am = photred_getairmass(file,obs=obs,update=update,silent=silent,stp=stp)
;
; By D.Nidever  February 2008
;-

function photred_getairmass,file,head=head,obs=obs,update=update,silent=silent,$
         recalculate=recalculate,error=error,stp=stp

COMMON photred,setup

undefine,error
am = -1.0

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - am = photred_getairmass(file,head=head,obs=obs,update=update,recalculate=recalculate,error=error,silent=silent,stp=stp)'
  return,-1
endif

;; Can't use input HEAD if multiple fits files or filter names input
if nfile gt 1 then undefine,head

; More than one filter name input
if nfile gt 1 then begin
  am = fltarr(nfile)
  for i=0,nfile-1 do am[i] = photred_getairmass(file[i],obs=obs,update=update,recalculate=recalculate,silent=silent)
  return,am
endif

;; No header input, read from fits file
if n_elements(head) eq 0 then begin
  ;; Check that the file exists  
  test = file_test(file)
  if test eq 0 then begin
    error = file+' NOT FOUND'
    if not keyword_set(silent) then print,error
    return,-1
  endif

  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
    head = PHOTRED_READFILE(file,exten=1,/header)
    ; Fix the NAXIS1/NAXIS2 in the header
    orig_head = head
    sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
    sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
  endif else begin
    head = PHOTRED_READFILE(file,/header)
  endelse
endif

if n_elements(obs) eq 0 then obs=''

; Get AIRMASS from header
am = SXPAR(head,'AIRMASS',count=nam,/silent)


; Compute airmass with RA/DEC, DATE and observatory
if nam eq 0 or float(am) lt 0.9 or keyword_set(recalculate) then begin

  ; Try WCS
  EXTAST,head,astr
  if n_elements(astr) gt 0 then begin
    ;ra = astr.crval[0]
    ;dec = astr.crval[1]
    naxis1 = sxpar(head,'NAXIS1',/silent)
    naxis2 = sxpar(head,'NAXIS2',/silent)

    ; use center of image
    head_xyad,head,naxis1/2,naxis2/2,ra,dec,/degree

  ; NO WCS, try RA/DEC
  endif else begin

    ; Get RA/DEC
    ra = SXPAR(head,'RA',count=nra,/silent)
    if (nra gt 0 and (strpos(ra,':'))[0] ne -1) then ra=sexig2ten(ra)*15.0
    dec = SXPAR(head,'DEC',count=ndec,/silent)
    if (ndec gt 0 and (strpos(dec,':'))[0] ne -1) then dec=sexig2ten(dec)

    if (nra eq 0) then begin
      error = 'RA not in header'
      if not keyword_set(silent) then print,error
      return,-1.0
    endif
    if (ndec eq 0) then begin
      error = 'DEC not in header'
      if not keyword_set(silent) then print,error
      return,-1.0
    endif

    ; Need float
    ra = double(ra)
    dec = double(dec)

  endelse


  ; Get DATE, YYYY-MM-DD
  date = PHOTRED_GETDATE(file,head=head)
  if (date eq '') then begin
    error = 'DATE ERROR'
    if not keyword_set(silent) then print,error
    return,-1.0
  endif
  datearr = strsplit(date,'-',/extract)
  year = long(datearr[0])
  month = long(datearr[1])
  day = long(datearr[2])

  ; Get UT-Time, HH:MM:SS.SSS
  uttime = PHOTRED_GETUTTIME(file,head=head)
  if (uttime eq '') then begin
    error = 'UT-TIME ERROR'
    if not keyword_set(silent) then print,error
    return,-1.0
  endif
  timearr = strsplit(uttime,':',/extract)
  hour = long(timearr[0])
  min = long(timearr[1])
  sec = float(timearr[2])

  ; Get observatory
  ;if obs eq '' then read,'What observatory: ',obs
  if obs eq '' then begin
    error = 'Error OBS(ERVATORY) not input'
    if not keyword_set(silent) then print,error
    return,-1.0
  endif
  OBSERVATORY,obs,obs_struct
  if obs_struct.observatory eq '' then begin
    error = 'Observatory "'+obs+'" NOT FOUND'
    if not keyword_set(silent) then print,error
    return,-1.0
  endif
  lat = obs_struct.latitude
  lon = obs_struct.longitude

  ; Calculate JD
  jd = JULDAY(month,day,year,hour,min,sec)

  ; Calculate airmass
  am = AIRMASS(jd,ra,dec,lat,lon)
  am = am[0]

  ; Update the HEADER
  if (keyword_set(update) and am gt 0.9) then begin
    SXADDPAR,head,'AIRMASS',am,' Added on '+systime()
    if not keyword_set(silent) then print,'AIRMASS=',strtrim(am,2),' added to '+file
    undefine,errmsg

    ;; Resource file exists
    dir = file_dirname(file)
    rfile = dir+'/.'+file
    if file_test(rfile) eq 1 then begin
      FITS_WRITE_RESOURCE,file,0,head
    ;; No resource file
    endif else begin
      ; Fpack compressed FITS file
      if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
        ; Put the original NAXIS1/2 values back
        sxaddpar,head,'NAXIS1',sxpar(orig_head,'NAXIS1')
        sxaddpar,head,'NAXIS2',sxpar(orig_head,'NAXIS2')
        ; Create temporary symbolic link to make modfits.pro
        ; think this is an ordinary FITS file
        tempfile = MAKETEMP('temp')
        FILE_LINK,file,tempfile+'.fits'
        MODFITS,tempfile+'.fits',0,head,exten_no=1,errmsg=errmsg
        FILE_DELETE,[tempfile,tempfile+'.fits'],/allow  ; delete temporary files  
      ; Normal FITS
      endif else begin
        MODFITS,file,0,head,errmsg=errmsg
      endelse
    endelse  ; no resource file
    if n_elements(errmsg) gt 0 then print,errmsg
  endif

endif

; Bad AIRMASS
if am lt 0.9 then begin
  error = 'Error calculating AIRMASS for '+file
  if not keyword_set(silent) then print,error
  return,-1.0
endif

if keyword_set(stp) then stop

return,am

end
