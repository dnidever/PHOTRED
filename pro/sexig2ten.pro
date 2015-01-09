; It appears that this program is totally
; unnecessary since IRAF will automatically
; convert a sexigesimal string to a base ten
; digit number if you use real().
;
; cl> print(real(-15:44:09.24))
; -15.7359
;

; The IDL program ten.pro would also do this
; pretty easily:  print,ten(-15,44,09.24)
;                 -15.7359
; even with sexigesimal string:  print,ten("-0:23:34")
;                                -0.39277778
; tenv.pro can be used with a vector of sexigesimal values


function sexig2ten,ra,stp=stp

nra = n_elements(ra)

; array input
if nra gt 1 then begin

  tenarr = dblarr(nra)
  for i=0.,nra-1 do tenarr(i) = sexig2ten(ra(i))

endif else  begin

  rgt = strtrim(ra,2)
  arr = strsplit(rgt,':',/extract)

  ; Still only one element.  Maybe try spaces
  if n_elements(arr) eq 1 then begin
    arr = strsplit(rgt,' ',/extract)
  endif

  ; Only one element, assume it's a number
  if n_elements(arr) eq 1 then begin
    tenarr = double(rgt)
    return,tenarr
  endif

  tenarr = double(ten(double(arr)))

  if strmid(rgt,0,1) eq '-' then tenarr=-abs(tenarr)

endelse ; not an array input

if keyword_set(stp) then stop

return,tenarr

end
