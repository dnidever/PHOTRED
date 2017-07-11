function photred_getairmass,file,obs=obs,update=update,silent=silent,$
         recalculate=recalculate,stp=stp


;+
;
; PHOTRED_GETAIRMASS
;
; This gets the AIRMASS for a FITS file.
; 
; INPUTS:
;  file      FITS filename
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
;
; USAGE:
;  IDL>am = photred_getairmass(file,obs=obs,update=update,silent=silent,stp=stp)
;
; By D.Nidever  February 2008
;-

COMMON photred,setup

am = -1.0

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  print,'Syntax - am = photred_getairmass(file,obs=obs,update=update,recalculate=recalculate,silent=silent,stp=stp)'
  return,-1
endif

; More than one filter name input
if nfile gt 1 then begin
  am = fltarr(nfile)
  for i=0,nfile-1 do am[i] = photred_getairmass(file[i],obs=obs,update=update,recalculate=recalculate,silent=silent)
  return,am
endif

test = file_test(file)
if test eq 0 then begin
  if not keyword_set(silent) then $
    print,file,' NOT FOUND'
  return,-1
endif

if n_elements(obs) eq 0 then obs=''

head = HEADFITS(file)

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
      if not keyword_set(silent) then $
        print,'RA not in header'
      return,-1.0
    endif
    if (ndec eq 0) then begin
      if not keyword_set(silent) then $
        print,'DEC not in header'
      return,-1.0
    endif

    ; Need float
    ra = double(ra)
    dec = double(dec)

  endelse


  ; Get DATE, YYYY-MM-DD
  date = PHOTRED_GETDATE(file)
  if (date eq '') then begin
    if not keyword_set(silent) then $
      print,'DATE ERROR'
    return,-1.0
  endif
  datearr = strsplit(date,'-',/extract)
  year = long(datearr[0])
  month = long(datearr[1])
  day = long(datearr[2])

  ; Get UT-Time, HH:MM:SS.SSS
  uttime = PHOTRED_GETUTTIME(file)
  if (uttime eq '') then begin
    if not keyword_set(silent) then $
      print,'UT-TIME ERROR'
    return,-1.0
  endif
  timearr = strsplit(uttime,':',/extract)
  hour = long(timearr[0])
  min = long(timearr[1])
  sec = float(timearr[2])

  ; Get observatory
  ;if obs eq '' then read,'What observatory: ',obs
  if obs eq '' then begin
    if not keyword_set(silent) then $
      print,'Error OBS(ERVATORY) not input'
    return,-1.0
  endif
  OBSERVATORY,obs,obs_struct
  if obs_struct.observatory eq '' then begin
    if not keyword_set(silent) then $
      print,'Observatory "',obs,'" NOT FOUND'
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
    MODFITS,file,0,head,errmsg=errmsg
    if n_elements(errmsg) gt 0 then print,errmsg
  endif

endif

; Bad AIRMASS
if am lt 0.9 then begin
  if not keyword_set(silent) then $
    print,'Error calculating AIRMASS for ',file
  return,-1.0
endif

if keyword_set(stp) then stop

return,am

end
