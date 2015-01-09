function odd,num

;+
;
; ODD
;
; Returns 1 if number is odd and 0 if number is even or zero
;
; INPUTS:
;  num     The number(s) to test for odd-ness.
;
; OUTPUTS:
;   1       Num is odd
;   0       Num is even
;  -1       A problem
;
; USAGE:
;  odd = odd(5)
;
; By D.Nidever  Feb 2008
;-

nnum = n_elements(num)
if nnum eq 0 then begin
  print,'Syntax - odd = odd(number)'
  return,-1
end
oddeven = fix(long64(num)/2LL ne float(num)*0.5)
nodd = total(oddeven)

return,oddeven
end
