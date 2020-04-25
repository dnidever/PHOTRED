;+
;
; FILE_WAIT
;
; Wait until a file exists.
;
; INPUTS:
;  file      Filename to check.
;  =wait     Wait time between checks.  Default is 5 sec.
;  =timeout  Stop trying and throw an error after this time.  Default
;               is 600 sec.
;  /silent   Don't print anything to the screen
;
; OUTPUTS:
;  None
;
; USAGE:
;  IDL>file_wait,'file.txt'
;
; By D. Nidever  April 2020
;-

pro file_wait,file,wait=wait,timeout=timeout,silent=silent

;; Not enough inputs
if n_elements(file) eq 0 then begin
  print,'Syntax - file_wait,file,wait=wait,timeout=timeout,silent=silent'
  return
endif

;; Defaults
if n_elements(wait) eq 0 then wait=5
if n_elements(timeout) eq 0 then timeout=600

;; While loop
t0 = systime(1)
while (file_test(file) eq 0) do begin
  if not keyword_set(silent) then message,file+' NOT FOUND.  Waiting '+strtrim(wait,2)+' sec.',/informational
  wait,wait
  if systime(1)-t0 gt timeout then message,'Timeout ('+strtrim(timeout,2)+' sec) reached on '+file
endwhile

end
