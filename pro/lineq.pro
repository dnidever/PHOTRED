pro lineq,x1,y1,x2,y2,m,b,noprint=noprint,curp=curp

if n_params() eq 0 and n_elements(curp) eq 0 then begin
  print,'Syntax - lineq,x1,y1,x2,y2,m,b,noprint=noprint,curp=curp'
  return
endif

; y = mx + b

if keyword_set(curp) then begin
  print,'Click on screen'
  cursor,x1,y1
  wait,0.1
  print,'Click on screen again'
  cursor,x2,y2
endif

m = (y2-y1)/(x2-x1)

b = y1 - m*x1

if not keyword_set(noprint) then print,'y = ',strtrim(m,2),'*x + ',strtrim(b,2)

end
