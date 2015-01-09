function scale,arr,oldrange,newrange

;+
; This function maps an array or image onto a new
; scale given two points on the old scale and
; the corresponding points on the new scale.
; The array is converted to double type.
; It's similar to BYTSCL.PRO except that you
; can set the bottom value as well.
; The ranges can be increasing or decreasing.
;
; INPUTS:
;  arr      The array of values to be scaled
;  oldrange Two-element array specifiying The original range which
;           will be scaled to newrange.
;  newrange Two-element array specifiying The new range which
;           the oldrange will be scaled to.
;
; OUTPUTS:
;  narr     The new scaled array
;
; USAGE:
;  IDL>arr2 = scale(arr,[0,1],[150,2000])
;
; By D.Nidever   March 2007
;-

; Not enough inputs
if n_elements(arr) eq 0 or n_elements(oldrange) lt 2 or n_elements(newrange) lt 2 then begin
  print,'Syntax - narr = scale(arr,oldrange,newrange)'
  return,-1
endif

; Does it flip around
signchange=1.0
if signs(oldrange[1]-oldrange[0]) ne signs(newrange[1]-newrange[0]) then signchange=-1.0 

; Scale
narr = range(newrange) * signchange*(double(arr)-oldrange[0])/range(oldrange) + newrange[0]
;narr = abs(newrange[1]-newrange[0]) * (double(arr)-oldrange[0])/abs(oldrange[1]-oldrange[0]) + newrange[0]

return,narr

end
