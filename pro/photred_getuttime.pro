function photred_getuttime,file,head=head,stp=stp


;+
;
; PHOTRED_GETUTTIME
;
; This gets the UT time for a FITS file.
; 
; INPUTS:
;  file      FITS filename
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The UT time is output in the format HH:MM:SS.SSS.
;  If there is an error then and empty string '' is output.
;
; USAGE:
;  IDL>ut = photred_getuttime(file)
;
; By D.Nidever  February 2008
;-

COMMON photred,setup

nfile = n_elements(file)
; Not enough inputs
if nfile eq 0 then begin
  print,'Syntax - uttime = photred_getuttime(file,stp=stp)'
  return,''
endif

; More than one filter name input
if nfile gt 1 then begin
  ut = strarr(nfile)
  for i=0,nfile-1 do ut[i] = photred_getuttime(file[i])
  return,ut
endif

test = file_test(file)
if test eq 0 then begin
  print,file,' NOT FOUND'
  return,''
endif

head = HEADFITS(file)

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
