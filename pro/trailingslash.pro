;+
;
; TRAILINGSLASH
;
; Add a trailing slash to a directory string, if needed.
;
; INPUTS:
;  dir     An input directory string.
;
; OUTPUTS:
;  dir     The directory string with a trailing slash added, if needed.
;
; USAGE:
;  IDL>dir = trailingslash(dir)
;
; By D. Nidever  Feb 2019
;-

function trailingslash,dir

;; Not enough inputs
if n_elements(dir) eq 0 then begin
  print,'Syntax - dir = trailingslash(dir)'
  return,''
endif
;; Not a string
if size(dir,/type) ne 7 then begin
  print,'DIR must be a STRING'
  return,''
endif

;; Get the last character
last = strmid(dir,0,1,/reverse_offset)
bd = where(last ne '/',nbd)
if nbd gt 0 then dir[bd] += '/'

return,dir

end
