function date2jd,input,mjd=mjd,sdss=sdss

;+
;
; DATE2JD
;
; Convert a date string to Julian date
;
; INPUTS:
;  input   The date string in YYYY-MM-DDTHH:MM:SS.S format.
;  /mjd    Output Modified Julian date
;  /sdss   use SDSS version of MJD
;
; OUTPUTS:
;  jd      Julian date for the input date string.
;
; USAGE:
;  IDL>jd = date2jd('2004-08-08T06:55:20')
;
; By D.Nidever  July 2010
;-

;date must be in this format
;  2004-08-08T06:55:20

; Not enough inputs
if n_elements(input) eq 0 then begin
  print,'Syntax - jd = date2jd(input,mjd=mjd,sdss=sdss)'
  return,-1
endif

in = input

if check_dateformat(input) eq 0 then return,-1

len = strlen(in)
icut = strpos(in,'T')
dat = strmid(in,0,icut)
tim = strmid(in,icut+1,len-icut-1)

; Getting date info
datarr = strsplit(dat,'-',/extract)
year = long(datarr[0])
month = long(datarr[1])
day = long(datarr[2])

; Getting time info
timarr = strsplit(tim,':',/extract)
hour = long(timarr[0])
min = long(timarr[1])
sec = double(timarr[2])

; Calculate JD
jd = JULDAY(month,day,year,hour,min,sec)

; MJD
if keyword_set(mjd) then jd-=2400000.5

; SDSS MJD
if keyword_set(mjd) and keyword_set(sdss) then jd+=0.3    ; use SDSS version

return,jd
end
