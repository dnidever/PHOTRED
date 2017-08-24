;+
;
; PHOTRED_GETDATE
;
; This gets the UT date for a FITS file.
; 
; INPUTS:
;  file      FITS filename
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The UT DATE is output in the format YYYY-MM-DD,
;  i.e. 2007-12-29.  If there is an error then
;  and empty string '' is output.
;
; USAGE:
;  IDL>ut = photred_getdate(file)
;
; By D.Nidever  May 2008
;-

function photred_getdate,file,head=head,stp=stp

COMMON photred,setup

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  print,'Syntax - uttime = photred_getdate(file,stp=stp)'
  return,''
endif

; More than one name input
if nfile gt 1 then begin
  date = strarr(nfile)
  for i=0,nfile-1 do date[i] = photred_getdate(file[i])
  return,date
endif

test = file_test(file)
if test eq 0 then begin
  print,file,' NOT FOUND'
  return,''
endif

if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then head=HEADFITS(file,exten=1) else $
  head = HEADFITS(file)

; Getting UT DATE
;----------------

; Different types of formats
;
; Swope
;DATE-OBS= '24Dec97'                                / UT DATE AT START OF FRAME  
;
;DATE-OBS= '10/02/00'                               / UT DATE AT START OF FRAME
;DATE-LCO= '10Feb00 '
;DATENEW = '2000-02-10'
;
;
; IMACS
;DATE-OBS= '2007-12-29'                 / UT date (start)
;UT-DATE = '2007-12-29'                 / UT date (start)
;
; LBC
;DATE_OBS= '2007-10-08T10:09:48.511' / Starting date of the observation
;
; MOSAIC
;DATE-OBS= '2007-08-17T08:03:58.2' / Date of observation start (UTC approximate)


; Try DATE-OBS
date = sxpar(head,'DATE-OBS',count=ndate,/silent)
if ndate eq 0 then $
; Try DATENEW
if ndate eq 0 then $
date = sxpar(head,'DATENEW',count=ndate,/silent)
; Try UT-DATE
if ndate eq 0 then $
date = sxpar(head,'UT-DATE',count=ndate,/silent)
; Try DATE_OBS
if ndate eq 0 then $
date = sxpar(head,'DATE_OBS',count=ndate,/silent)


; NO date found
if ndate eq 0 then begin
  print,'NO DATE IN HEADER'
  return,''
endif


; Check format
date_orig = date

; Check "T" format, i.e. 2007-10-08T10:09:48.511
Tind = strpos(date,'T')
if Tind[0] ne -1 then begin
  time = first_el(strsplit(date,'T',/extract),/last)
  date = first_el(strsplit(date,'T',/extract))
endif

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

; Check "-" format, i.e. 2007-12-29
;----------------------------------
; This is also the "output" format
if (ndash gt 0) then begin
  date = strtrim(date,2)
  datearr = strsplit(date,'-',/extract)
  syear = datearr[0]
  smonth = datearr[1]
  sday = datearr[2]

  year_isnum = VALID_NUM(syear,year)
  month_isnum = VALID_NUM(smonth,month)
  day_isnum = VALID_NUM(sday,day)

  ; All valid numbers, day<=31,  month <=12
  if (year_isnum eq 1 and month_isnum eq 1 and day_isnum eq 1) and $
     (day ge 1 and day le 31 and month ge 1 and month le 12) then begin

    date_out = string(year,format='(I4)')+'-'+string(month,format='(I02)')+$
               '-'+string(day,format='(I02)')
  endif

endif


; Check "/" format, i.e. 10/02/00, DD/MM/YY
;------------------------------------------
if (nslash gt 0) and (date_out eq '') then begin
  date = strtrim(date,2)
  datearr = strsplit(date,'/',/extract)

  sday   = datearr[0]
  smonth = datearr[1]
  syear = datearr[2]

  year_isnum = VALID_NUM(syear,year)
  month_isnum = VALID_NUM(smonth,month)
  day_isnum = VALID_NUM(sday,day)

  ; Year in two digit format
  if strlen(syear) eq 2 then begin

    ; Get current year
    CALDAT,systime(/julian),cmonth,cday,cyear
    cyear2 = cyear-2000L
    ; if "year" is great than current 2-digit year then
    ;   add 1900 else add 2000
    if year gt cyear2 then year=year+1900 else year=year+2000L
  endif

  ; All valid numbers, day<=31,  month <=12
  if (year_isnum eq 1 and month_isnum eq 1 and day_isnum eq 1) and $
     (day ge 1 and day le 31 and month ge 1 and month le 12) then begin

    date_out = string(year,format='(I4)')+'-'+string(month,format='(I02)')+$
               '-'+string(day,format='(I02)')
  endif
endif


; Check "DDFebYY" format, i.e. 10Feb00
;-------------------------------------
if (nisletter gt 0) and (date_out eq '') then begin
  date = strtrim(date,2)
  len = strlen(date)

  ; "DDMMMYY" format
  if (VALID_NUM(strmid(date,1,1)) eq 1) then begin
    sday = strmid(date,0,2)
    smonth_letter = strmid(date,2,3)
    syear = strtrim(strmid(date,5,4),2)

  ; "DMMMYY" format, one day digit
  endif else begin
    sday = strmid(date,0,1)
    smonth_letter = strmid(date,1,3)
    syear = strtrim(strmid(date,4,4),2)
  endelse

  ; Getting month
  months = ['jan','feb','mar','apr','may','jun','jul','aug','sep',$
            'oct','nov','dec']
  monind = where(strlowcase(smonth_letter) eq months,nmonind)
  month = -1
  if nmonind gt 0 then month=monind[0]+1

  year_isnum = VALID_NUM(syear,year)
  day_isnum = VALID_NUM(sday,day)

  ; Year in two digit format
  if (strlen(syear) eq 2 and year_isnum eq 1) then begin

    ; Get current year
    CALDAT,systime(/julian),cmonth,cday,cyear
    cyear2 = cyear-2000L
    ; if "year" is great than current 2-digit year then
    ;   add 1900 else add 2000
    if year gt cyear2 then year=year+1900 else year=year+2000L
  endif

  ; All valid numbers, day<=31,  month <=12
  if (year_isnum eq 1 and day_isnum eq 1) and $
     (day ge 1 and day le 31 and month ge 1 and month le 12) then begin

    date_out = string(year,format='(I4)')+'-'+string(month,format='(I02)')+$
               '-'+string(day,format='(I02)')
  endif
endif


; NO date found
if date_out eq '' then begin
  print,'NO DATE IN HEADER'
endif


if keyword_set(stp) then stop

return,date_out

end
