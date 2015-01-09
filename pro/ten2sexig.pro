function ten2sexig,ra,stp=stp

nra = n_elements(ra)

; array input
if (nra gt 1) then begin
  sexig = strarr(nra)
  for i=0.,nra-1 do sexig(i) = ten2sexig(ra(i))
endif else begin

rgt = double(ra)
;rgt = float(ra)
sign = signs(rgt)
;sign = rgt/abs(rgt)
rgt = abs(rgt)

; getting hour
ihour = long(rgt)
;ihour = ihour * sign

; getting minute
rgt = rgt-long(rgt)
imin = long(rgt*60.0d0)

; getting second
rgt = rgt*60.0d0-long(rgt*60.0d0)
rsec = rgt*60.0d0

; putting it together
ssign = ''
if sign eq -1 then ssign = '-'
;ssign = strmid(strtrim(sign,2),0,1)
shour = strtrim(ihour,2)
if (abs(ihour) lt 10) then shour = '0' + shour
smin = strtrim(imin,2)
if (imin lt 10) then smin = '0' + smin
ssec = strtrim(string(rsec,format='(F30.20)'),2)
;ssec = strtrim(rsec,2)
if (rsec lt 10.) then ssec = '0' + ssec
ssec = strmid(ssec,0,6)

sexig = ssign + shour + ':' + smin + ':' + ssec

endelse ; not an array input

if keyword_set(stp) then stop

return,sexig

end
