pro setdisp,colorin,silent=silent

; Set the display properly

if n_elements(colorin) gt 0 then clr=colorin[0] else clr=39
if clr lt 0 or clr gt 40 then begin
  print,'Color table must be from 0 to 40 - Using 39'
  clr = 39
endif

if not keyword_set(silent) then $
  print,'Loading color table '+strtrim(long(clr),2)
loadct,clr,/silent
if not keyword_set(silent) then $
  print,'Setting DECOMPOSED=0'
device,decomp=0
psym8

end
