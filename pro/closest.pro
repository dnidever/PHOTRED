function closest,value,arr,ind=ind,upper=upper,lower=lower,stp=stp

;+
;
; CLOSEST
;
; This function finds the array element that is closest to value.
;
; INPUTS:
;  value   The value that is desired
;  arr     The array to check
;  /upper  Want the element closest to value, but less than it
;  /lower  Want the element closest to value, but less than it
;  /stp    Stop at the end of the program
;
; OUTPUTS:
;  The value of the array element closest to the value input
;  =ind  The index of the array element closest to the value input
;
; USAGE:
;  IDL>val = closest(1.0,array,ind=ind)
;
; By D.Nidever   2006
;-

; Not enough inputs
if n_params() lt 2 then begin
  print,'Syntax - out = closest(value,arr,ind=ind,upper=upper,lower=lower,stp=stp)'
  return,-1
end

; Checking the parameters
;zparcheck,'CLOSEST',value,1,[2,3,4,5],0,'VALUE not correct'
;zparcheck,'CLOSEST',arr,2,[2,3,4,5],[0,1],'ARR not correct'

; Normal
if not keyword_set(upper) and not keyword_set(lower) then begin
  ind = first_el(minloc(abs(value-arr)))

; /upper or /lower set
endif else begin

  ; /upper
  if keyword_set(upper) then begin
    gd = where(arr ge value,ngd)
    if ngd gt 0 then ind = gd(first_el(minloc(abs(value-arr(gd)))))
    if ngd eq 0 then ind = first_el(minloc(abs(value-arr)))
  endif

  ; /lower
  if keyword_set(lower) then begin
    gd = where(arr le value,ngd)
    if ngd gt 0 then ind = gd(first_el(minloc(abs(value-arr(gd)))))
    if ngd eq 0 then ind = first_el(minloc(abs(value-arr)))
  endif
endelse

if keyword_set(stp) then stop

return,arr(ind)

end
