pro touchzero,files,stp=stp

;+
;
; TOUCHZERO
;
; This creates an empty file
;
; INPUTS:
;  files   Files to create
;  /stp    Stop at the end of the program.
;
; OUTPUTS:
;  The emtpy files
;
; USAGE:
;  IDL>touchzero,'file.txt'
;
; By D.Nidever  March 2008
;-

; No inputs
nfiles = n_elements(files)
if (nfiles eq 0) then begin
  print,'Syntax - touchzero,files,stp=stp'
  return
endif

; Loop through the files
FOR i=0,nfiles-1 do begin

  OPENW,unit,/get_lun,files[i]
  CLOSE,unit
  FREE_LUN,unit

END


;stop

if keyword_set(stp) then stop

end
