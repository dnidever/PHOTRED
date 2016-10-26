pro curp,x,y,xr=xr,yr=yr,position=position

;; Make the plot to create the coordinate system
;if keyword_set(xr) and keyword_set(yr) then begin
;  xmult=1
;  ymult=1
;  if keyword_set(xflip) then xmult=-1
;  if keyword_set(yflip) then ymult=-1
;
;  ; Doing the plot
;  pxr = xr
;  pyr = yr
;  if keyword_set(xflip) then pxr=-reverse(xr)
;  if keyword_set(yflip) then pyr=-reverse(yr)
;  plot,dist(10),/nodata,/noerase,xr=pxr,yr=pyr,xs=1,ys=1,pos=position
;endif

; cursor and print

cursor,x,y

;if keyword_set(xr) and keyword_set(yr) then begin
;  plot,[x],[y],/noerase,xr=pxr,yr=pyr,xs=1,ys=1,pos=position,ps=1
;  if keyword_set(xflip) then x=-x
;  if keyword_set(yflip) then y=-y
;endif
print,x,y



end
