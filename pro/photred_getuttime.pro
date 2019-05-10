;+
;
; PHOTRED_GETUTTIME
;
; This gets the UT time for a FITS file.
; 
; INPUTS:
;  file      FITS filename
;  =head     Use this header array instead of reading FITS file.   
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The UT time is output in the format HH:MM:SS.SSS.
;  If there is an error then and empty string '' is output.
;  =error    The error message if one occurred. 
;
; USAGE:
;  IDL>ut = photred_getuttime(file)
;
; By D.Nidever  February 2008
;-

function photred_getuttime,file,head=head,error=error,stp=stp

COMMON photred,setup

undefine

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - uttime = photred_getuttime(file,stp=stp)'
  return,''
endif

;; Can't use input HEAD if multiple fits files input
if nfile gt 1 then undefine,head

; More than one filter name input
if nfile gt 1 then begin
  ut = strarr(nfile)
  for i=0,nfile-1 do ut[i] = photred_getuttime(file[i])
  return,ut
endif

;; Header not input, read from FITS file
if n_elements(head) eq 0 then begin
  ;; Check that the file exists
  test = file_test(file)
  if test eq 0 then begin
    error = file+' NOT FOUND'
    print,error
    return,''
  endif

  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then head=PHOTRED_READFILE(file,exten=1,/header) else $
    head = PHOTRED_READFILE(file,/header)
endif

; Getting UT TIME
;----------------
; Try TIME-OBS
ut = sxpar(head,'TIME-OBS',/silent)
ind = strpos(ut,':')
if ind[0] eq -1 then ut='0'

; Try UT
if strtrim(ut,2) eq '0' then begin
  ut = sxpar(head,'UT',/silent)
  ind = strpos(ut,':')
  if ind[0] eq -1 then ut='0'
endif

; Try DATE-OBS
if strtrim(ut,2) eq '0' then begin
  ut = sxpar(head,'DATE-OBS',/silent)
  ind = strpos(ut,':')
  if ind[0] eq -1 then ut='0'
  ; Is this in the format YYYY-MM-DDTHH:MM:SS.SSS ?
  if strpos(ut,'T') ne -1 then ut = first_el(strsplit(ut,'T',/extract),/last)
endif

; Try DATE
if strtrim(ut,2) eq '0' then begin
  ut = sxpar(head,'DATE',/silent)
  ind = strpos(ut,':')
  if ind[0] eq -1 then ut='0'
  ; Is this in the format YYYY-MM-DDTHH:MM:SS.SSS ?
  if strpos(ut,'T') ne -1 then ut = first_el(strsplit(ut,'T',/extract),/last)
endif

; Try DATE_OBS (LBC data)
if strtrim(ut,2) eq '0' then begin
  ut = sxpar(head,'DATE_OBS',/silent)
  ind = strpos(ut,':')
  if ind[0] eq -1 then ut='0'
  ; Is this in the format YYYY-MM-DDTHH:MM:SS.SSS ?
  if strpos(ut,'T') ne -1 then ut = first_el(strsplit(ut,'T',/extract),/last)
endif

; Try UT-TIME
if (strtrim(ut,2) eq '0') then begin
  ut = SXPAR(head,'UT-TIME',/silent)
  ut = strtrim(ut,2)
  ind = strpos(ut,':')
  if ind[0] eq -1 then ut='0'
endif

; Try UTSTART (old Swope data)
if strtrim(ut,2) eq '0' then begin
  utstart = sxpar(head,'UTSTART',/silent)
  utstart = strtrim(utstart,2)
  ind = strpos(ut,':')
  if ind[0] eq -1 then ut = REPSTR(utstart,' ',':')
endif

; NO time found
if strtrim(ut,2) eq '0' then begin
  print,'NO UT TIME IN HEADER'
  return,''
endif

time = ut

if keyword_set(stp) then stop

return,time

end
