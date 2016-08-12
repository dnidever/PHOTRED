;+
;
; PHOTRED_GETMJD
;
; This gets the MJD night number for a given file.  Also needs to know
; the observatory.
; 
; INPUTS:
;  file      FITS filename
;  obs       The observatory name.
;  
;  =dateobs  The DATE-OBS string.  This is used to determine the MJD.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The Modified Julian Date properly modified for the observatory
;  local time.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>mjd = photred_getmjd(file,'ctio')
;
; By D.Nidever  January 2016
;  much of it copied from APOGEE getmjd5.pro
;-
function photred_getmjd,file,obs,dateobs=dateobs,stp=stp,error=error

undefine,error

nfile = n_elements(file)
nobs = n_elements(obs)
ndateobs = n_elements(dateobs)
; Not enough inputs
if (nfile eq 0 and ndateobs eq 0) or nobs eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - mjd = photred_getmjd(file,obs,dateobs=dateobs,stp=stp,error=error)'
  return,-1
endif

; Get observatory
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

; More than one filename name input
if nfile gt 1 then begin
  mjd = lonarr(nfile)
  for i=0,nfile-1 do mjd[i] = photred_getmjd(file[i],obs)
  return,mjd
endif

test = file_test(file)
if test eq 0 and ndateobs eq 0 then begin
  error = file+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return,-1
endif

; Load file
if ndateobs eq 0 then head = HEADFITS(file)

; --- Get the DATE ---
if ndateobs eq 0 then begin
  date = photred_getdate(file)
endif else begin  ; use DATE-OBS
  dum = strsplit(dateobs,'T',/extract)
  date = dum[0]
endelse

; parse the date information
datearr = strsplit(date,'-',/extract)

if n_elements(datearr) lt 3 then begin
  error = 'DATE = '+date+' NOT in the correct format. YYYY-MM-DD'
  if not keyword_set(silent) then print,error
  return,-1
endif

; Check each number
year = datearr[0]
if VALID_NUM(year,/integer) eq 0 then begin
  error = 'YEAR = '+year+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
year_num = long(year)

month = datearr[1]
if VALID_NUM(month,/integer) eq 0 then begin
  error = 'MONTH = '+month+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
month_num = long(month)
if month_num lt 1 or month_num gt 12 then begin
  error = 'MONTH = '+month+' MUST be 1-12'
  if not keyword_set(silent) then print,error
  return,-1
endif

day = datearr[2]
if VALID_NUM(day,/integer) eq 0 then begin
  error = 'DAY = '+day+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
day_num = long(day)
if day_num lt 1 or day_num gt 31 then begin
  error = 'DAY = '+day+' Must be 1-31'
  if not keyword_set(silent) then print,error
  return,-1
endif

; --- Get the TIME ---
if ndateobs eq 0 then begin
  time = photred_getuttime(file)
endif else begin  ; use DATE-OBS
  dum = strsplit(dateobs,'T',/extract)
  time = dum[1]
endelse

; parse the time information
timearr = strsplit(time,':',/extract)

if n_elements(timearr) lt 3 then begin
  error = 'TIME = '+date+' NOT in the correct format. HH:MM:SS.S'
  if not keyword_set(silent) then print,error
  return,-1
endif

hour = timearr[0]
if VALID_NUM(hour,/integer) eq 0 then begin
  error = 'hour = '+hour+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
hour_num = long(hour)
if hour_num lt 0 or hour_num gt 24 then begin
  error = 'HOUR = '+hour+' Must be 1-24'
  if not keyword_set(silent) then print,error
  return,-1
endif

minute = timearr[1]
if VALID_NUM(minute,/integer) eq 0 then begin
  error = 'MIN = '+minute+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
minute_num = long(minute)
if minute_num lt 0 or minute_num gt 60 then begin
  error = 'MIN = '+minute+' Must be 1-60'
  if not keyword_set(silent) then print,error
  return,-1
endif

sec = timearr[2]
if VALID_NUM(sec) eq 0 then begin
  error = 'SEC = '+sec+' is NOT a valid number'
  if not keyword_set(silent) then print,error
  return,-1
endif
sec_num = float(sec)
if sec_num lt 0 or sec_num gt 60 then begin
  error = 'SEC = '+sec+' Must be 1-60'
  if not keyword_set(silent) then print,error
  return,-1
endif


; Getting Julian Date
;---------------------
jd = JULDAY(month_num, day_num, year_num, hour_num, minute_num, sec_num)

; Modified Julian Date
;----------------------
; The modified Julian date is MJD = JD - 2400000.5
; The Julian day starts at NOON, while MJD starts at midnight
MJD = JD - 2400000.5

; Correct for the local observatory time
;---------------------------------------
; correct to noon local time, MJD starts at midnight
MJD += 0.5    ; correct to noon
MJD -= obs_struct.tz/24.   ; time zone, number of hours *west* of Greenwich
;  correct to local time

; Convert to long, trim the decimal
MJD = long(MJD)  ; clip the decimals

if keyword_set(stp) then stop

return,mjd

end
