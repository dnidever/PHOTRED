pro psym8,size,fill=fill,thick=thick,square=square,diamond=diamond,$
          vline=vline,hline=hline,npts=npts
;+
;
; psym8,size,fill=fill,thick=thick,square=square,diamond=diamond,
;          vline=vline,hline=hline,npts=npts
;
;-

if n_elements(fill) eq 0 then fill=1
if not keyword_set(thick) then thick=1
if not keyword_set(size) then size=1

if keyword_set(diamond) then begin
  x=size*[-1.,0.,1.,0.,-1]
  y=size*[0.,1.,0.,-1.,0.]
  usersym,x,y,fill=fill,thick=thick
endif else if keyword_set(square) then begin
  x=size*[-1.,1.,1.,-1.,-1]
  y=size*[1.,1.,-1.,-1.,1.]
  usersym,x,y,fill=fill,thick=thick
endif else if keyword_set(vline) then begin
  x=size*[0.,0.]
  y=size*[-1.,1.]
  usersym,x,y,fill=fill,thick=thick
endif else if keyword_set(hline) then begin
  x=size*[-1.,1.]
  y=size*[0.,0.]
  usersym,x,y,fill=fill,thick=thick
endif else begin
  if keyword_set(npts) then n=npts else n=32
  if n gt 49 then n=49
  a=dindgen(n)*(!dpi*2./(n-1.))	;defining psym 8
  usersym,size*cos(a),size*sin(a),fill=fill,thick=thick
endelse


end
