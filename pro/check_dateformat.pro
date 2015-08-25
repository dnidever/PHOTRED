function check_dateformat,dateobs

; This checks that the DATEOBS FITS header format is
; correct.  YYYY-MM-DDTHH:MM:SS.S


; Check "T" format, i.e. 2007-10-08T10:09:48.511
Tind = strpos(dateobs,'T')
if Tind[0] ne -1 then begin
  time = first_el(strsplit(dateobs,'T',/extract),/last)
  date = first_el(strsplit(dateobs,'T',/extract))
; No "T" found
endif else begin
  return,0
endelse

; Numbers, letters, dash, and slash
bindate = byte(date)  ; binary date
nchar = n_elements(bindate)
date_chars = strarr(nchar)
for i=0,nchar-1 do date_chars[i] = string(bindate[i])
isnum = intarr(nchar)
for i=0,nchar-1 do isnum[i] = VALID_NUM(date_chars[i])
; Letters are A-Z (65-90) and a-z (97-122)
isletter = where((bindate ge 65B and bindate le 90B) or $
                 (bindate ge 97B and bindate le 122B),nisletter)
dash = where(bindate eq 45B,ndash)    ; '-' = 45B
slash = where(bindate eq 47B,nslash)  ; '/' = 47B

date_out = ''

; No dash "-"
if ndash eq 0 then return,0

; Check date format
date = strtrim(date,2)
datearr = strsplit(date,'-',/extract)
if n_elements(datearr) ne 3 then return,0  ; need 3 elements
syear = datearr[0]
smonth = datearr[1]
sday = datearr[2]

year_isnum = VALID_NUM(syear,year)
month_isnum = VALID_NUM(smonth,month)
day_isnum = VALID_NUM(sday,day)

date_okay = 0   ; bad until proven good

; All valid numbers, day<=31,  month <=12
if (year_isnum eq 1 and month_isnum eq 1 and day_isnum eq 1) and $
   (day ge 1 and day le 31 and month ge 1 and month le 12) then begin

   date_okay = 1
endif



; No colon ":"
colon_ind = strpos(time,':')
if colon_ind[0] eq -1 then return,0

; Check time format
timarr = strsplit(time,':',/extract)
if n_elements(timarr) ne 3 then return,0  ; need 3 elements
shour = timarr[0]
sminute = timarr[1]
ssec = timarr[2]

hours=-1 & minute=-1 & sec=-1
hour_isnum = VALID_NUM(shour,hour)
minute_isnum = VALID_NUM(sminute,minute)
sec_isnum = VALID_NUM(ssec,sec)

time_okay = 0   ; bad until proven good

; All valid numbers, hour <=24,  min<=59,  sec < 60.0
if (hour_isnum eq 1 and minute_isnum eq 1 and sec_isnum eq 1) and $
   (hour ge 0 and hour lt 24 and minute ge 0 and minute le 59 and sec ge 0.0 and sec lt 60.0) then begin

   time_okay = 1
endif

okay = (date_okay AND time_okay)  ; just both be oky

return,okay
end
