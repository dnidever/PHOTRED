function slope, y, x, acc=acc, error=error

;+
;
; SLOPE
;
; This program calculate the slope of an array
;
; INPUTS:
;  y       The array of y values
;  x       The array of x values
;  /acc    Calculate it more accurately
;
; OUTPUTS:
;  result  The slope of the array.
;  =error  The error message if one occured.
;
; USAGE:
;  IDL> result=slope(Yarray [,Xarray])
;
; By D. Nidever   2005?
;-

ny = n_elements(y)

if N_params() LT 1 or ny eq 0 then begin   
  print, 'Syntax - result=slope(Yarray [,Xarray])'
  print, 'dY(i) = Y(i+1)-Y(i)'
  error = 'Not enough inputs'
  return,-1
endif

n=n_elements(y)

if n eq 1 then begin
  error = 'Need more than 1 point'
  return,-1
endif

if not keyword_set(acc) then begin
  dely=y[1:n-1]-y[0:n-2]

  if keyword_set(x) then dely=dely/(x[1:n-1]-x[0:n-2])
endif else begin

  xx = dindgen(n)
  yy = spline(xx,y,xx+0.01)
  dely = (yy-y)/0.01

  if keyword_set(x) then begin
    x2 = spline(xx,x,xx+0.01)
    delx = (x2-x)/0.01
    dely = dely/delx
  endif
endelse

return,dely

end
