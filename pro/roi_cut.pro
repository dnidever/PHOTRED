;+
;
; ROI_CUT
;
; Use cuts in a 2D plane to select points from arrays.
; You can use CLICKER.PRO to obtain the xcut and ycut
; arrays by clicking on the screen.
;
; Uses polyfillv.pro.  Must be on a regular 2D grid.
; So I'll have to round the veloities, but no big deal
; because the resolution is ~1 km/s anyway.
;
; You can use clicker.pro to get the xcut/ycut arrays
;
; INPUTS:
;  xcut  Array of x-values for the cut
;  ycut  Array of y-values for the cut
;  x     Array of x-values that should be cut
;  y     Array of y-values that should be cut
;  =fac  A factor to multiply the values by to improve
;         the cut.
;  /silent  Don't print anything to the screen.
;
; OUTPUTS:
;  ind     The indices of values OUTSIDE the cut
;  cutind  The indices of values INSIDE the cut
;
; USAGE:
;  roi_cut,xcut,ycut,x,y,ind,cutind,fac=fac
;
; By D.Nidever 2006
;-

function polyfillv_exists

Catch, theError
IF theError NE 0 THEN BEGIN
  Catch, /CANCEL
  RETURN,0
ENDIF

xc = [1.0,4.0,4.0,1.0]
yc = [1.0,1.0,4.0,4.0]
out = polyfillv(xc,yc,10,10)

return,1

end

;--------

pro roi_cut,xcut,ycut,x,y,ind,cutind,fac=fac,silent=silent

; Not enough inputs
if n_params() lt 4 then begin
  print,'Syntax - roi_cut,xcut,ycut,x,y,ind,cutind,fac=fac,silent=silent'
  return
endif

xmin = min(xcut)
xmax = max(xcut)
ymin = min(ycut)
ymax = max(ycut)

nx = n_elements(x)

; Make a first cut with the X/Y minima/maxima
gd = where(x ge xmin and x le xmax and y ge ymin and y le ymax,ngd)
bd = where(x lt xmin or x gt xmax or y lt ymin or y gt ymax,nbd)
if ngd eq 0 then begin
  ind = findgen(nx)
  undefine,cutind
  return
end
x2 = x[gd]
y2 = y[gd]

; Making temporary copies
xcut2 = xcut
ycut2 = ycut

; Multiplying all numbers by a large factor
if (xmax-xmin) lt 100. or (ymax-ymin) lt 100. or n_elements(fac) gt 0 then begin
  if n_elements(fac) eq 0 then fac = round(100./(xmax-xmin))
  x2 = x2*fac
  y2 = y2*fac
  xcut2 = xcut2*fac
  ycut2 = ycut2*fac
  if not keyword_set(silent) then $
    print,'Multiplying numbers by ',strtrim(fac,2)
endif


; Now we need to round them to the nearest integer
xcut2 = round(xcut2)
ycut2 = round(ycut2)
x3 = round(x2)
y3 = round(y2)

xmin2 = min(xcut2)
xmax2 = max(xcut2)
ymin2 = min(ycut2)
ymax2 = max(ycut2)

nx = xmax2-xmin2+1
ny = ymax2-ymin2+1

; Getting the indices of the points on the 2D grid.
xind = x3-xmin2
yind = y3-ymin2

; Using polyfillv to get the good indices of the cut
; using positions relative to (xmin2,ymin2)
if polyfillv_exists() eq 1 then begin
  ind1d = polyfillv(xcut2-xmin2,ycut2-ymin2,nx,ny)
endif else begin
  ; Make the arrays
  xarr = (findgen(nx)#replicate(1,ny))(*)
  yarr = (replicate(1,nx)#findgen(ny))(*)
  ind1d = inside(xarr,yarr,xcut2-xmin2,ycut2-ymin2)
endelse
ind2d = array_indices(fltarr(nx,ny),ind1d)
indx = reform(ind2d(0,*)) 
indy = reform(ind2d(1,*))

; Looping through the 
mask = fltarr(nx,ny)
(mask)(ind1d) = 1

; Getting the 1D indices for the input points
pind2d = fltarr(2,ngd)
pind2d(0,*) = xind
pind2d(1,*) = yind
pind1d = array_indices2(fltarr(nx,ny),pind2d)

; Getting the points within the cut region
out2 = where((mask)(pind1d) eq 0,nout2)
in2 = where((mask)(pind1d) eq 1,nin2)

; The indices of stars inside the region
; ind - outside the cut
; cutind - inside the cut
undefine,ind,cutind
if nbd gt 0 then push,ind,bd
if nout2 gt 0 then push,ind,gd[out2]
if nin2 gt 0 then push,cutind,gd[in2]
;ind = [bd,gd[out2]]   ; outside 
;cutind = gd[in2]      ; inside

if n_elements(ind) eq 0 then ind=-1
if n_elements(cutind) eq 0 then cutind=-1

; Sorting the indices
if n_elements(ind) gt 0 then begin
  si = sort(ind)
  ind = ind[si]
endif
if n_elements(cutind) gt 0 then begin
  si = sort(cutind)
  cutind = cutind[si]
endif

;stop

end
