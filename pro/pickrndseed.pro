;+
;
; PIKCRNDSEED
;
; This program picks a random number generator
; seed in a reproducible way given an input
; number (i.e. based on same data).
;
; INPUTS:
;  data   Some data that is used to the create
;          the seed.  Only the first element is
;          used.
;
; OUTPUTS:
;  seed   The seed to use for a random number
;           generator (e.g. randomn or randomu).
;           Currently this will be a number
;           between 1 and 1000.
;
; USAGE:
;  IDL>seed = pickrndseed(data)
;
; By D. Nidever  Sep 2017
;-

function pickrndseed,data

; No value input, use seed=1
if n_elements(data) eq 0 then return,1L

; Initialize seed
;  Use the first value of DATA to set seed
;  so it's always the same and reproducible for
;  the same data set.
val = data[0]
type = size(val,/type)
; Check if we have a number that we can use
dum = where(type eq [1,2,3,4,5,6,9,12,13,14,15],isnum)
if type eq 7 then isnum=valid_num(val,val)
; Number
if isnum eq 1 then begin
  ; Convert all types to float
  val = fix(val,type=4)
  ; It's a finite number
  if finite(val) eq 1 then begin
     frac = abs(alog10(abs(val)))-floor(abs(alog10(abs(val))))
     seed = round(frac*1000) > 1  ; number between 1 and 1000
  endif else seed=1
; Non-number
endif else begin
  seed = 1
endelse

return,seed

end
