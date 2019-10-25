;+
;
; ALLFRAME_ADXYINTERP
;
; Instead of transforming the entire large 2D RA/DEC
; arrays to X/Y do a sparse grid and perform linear
; interpolation to the full size.                                
;
; INPUTS:
;  head    The FITS header with the WCS.
;  rr      The 2D RA array.
;  dd      The 2D DEC array.
;  =nstep  The undersampling to use, default nstep=10.
;  /xyad   Perform X/Y->RA/DEC conversion instead.
;            In this case the meaning of the coordinate
;            arrays are flipped, i.e. rr/dd<->xx/yy
;
; OUTPUTS:
;  xx      The 2D X array.
;  yy      The 2D Y array.
;
; USAGE:
;  IDL>allframe_adxyinterp,head,ra,dec,xx,yy,nstep=10
;
; By D. Nidever  Oct. 2016
;-

pro allframe_adxyinterp,head,rr,dd,xx,yy,nstep=nstep,xyad=xyad

undefine,xx
undefine,yy

; Not enough inputs
if n_elements(head) eq 0 or n_elements(rr) eq 0 or n_elements(dd) eq 0 then begin
  print,'Syntax - allframe_adxyinterp,head,rr,dd,xx,yy,nstep=nstep,xyad=xyad'
  return
endif

if n_elements(nstep) eq 0 then nstep = 10
sz = size(rr)
nx = sz[1]
ny = sz[2]
nxs = (nx-1)/nstep + 1
nys = (ny-1)/nstep + 1

;; Small dimensions, just transform the full grid
if nx le nstep or ny le nstep then begin
  if not keyword_set(xyad) then begin
    HEAD_ADXY,head,rr,dd,xx,yy,/deg
  endif else begin
    HEAD_XYAD,head,rr,dd,xx,yy,/deg
  endelse
  return
endif

; Subsample the RA/DEC arrays
rrs = rr[0:nx-1:nstep,0:ny-1:nstep]
dds = dd[0:nx-1:nstep,0:ny-1:nstep]
if not keyword_set(xyad) then begin
  HEAD_ADXY,head,rrs,dds,xxs,yys,/deg
endif else begin
  HEAD_XYAD,head,rrs,dds,xxs,yys,/deg
endelse

; Start final arrays
xx = dblarr(nx,ny)
yy = dblarr(nx,ny)

; Use CONGRID to perform the linear interpolation
;   congrid normally does nx_out/nx_in scaling
;   if /minus_one set it does (nx_out-1)/(nx_in-1)
ixx = CONGRID(xxs,(nxs-1)*nstep+1,(nys-1)*nstep+1,/interp,/minus_one)
iyy = CONGRID(yys,(nxs-1)*nstep+1,(nys-1)*nstep+1,/interp,/minus_one)
xx[0:(nxs-1)*nstep,0:(nys-1)*nstep] = ixx
yy[0:(nxs-1)*nstep,0:(nys-1)*nstep] = iyy

; Deal with right edge
if (nxs-1)*nstep+1 lt nx then begin
  ; Make a short grid in X at the right edge
  rrs_rt = rr[nx-nstep-1:nx-1:nstep,0:ny-1:nstep]
  dds_rt = dd[nx-nstep-1:nx-1:nstep,0:ny-1:nstep]
  if not keyword_set(xyad) then begin
    HEAD_ADXY,head,rrs_rt,dds_rt,xxs_rt,yys_rt,/deg
  endif else begin
    HEAD_XYAD,head,rrs_rt,dds_rt,xxs_rt,yys_rt,/deg
  endelse
  ixx_rt = CONGRID(xxs_rt,nstep+1,(nys-1)*nstep+1,/interp,/minus_one)
  iyy_rt = CONGRID(yys_rt,nstep+1,(nys-1)*nstep+1,/interp,/minus_one)
  xx[nx-nstep-1:nx-1,0:(nys-1)*nstep] = ixx_rt
  yy[nx-nstep-1:nx-1,0:(nys-1)*nstep] = iyy_rt
endif
; Deal with top edge
if (nys-1)*nstep+1 lt ny then begin
  ; Make a short grid in Y at the top edge
  rrs_tp = rr[0:nx-1:nstep,ny-nstep-1:ny-1:nstep]
  dds_tp = dd[0:nx-1:nstep,ny-nstep-1:ny-1:nstep]
  if not keyword_set(xyad) then begin
    HEAD_ADXY,head,rrs_tp,dds_tp,xxs_tp,yys_tp,/deg
  endif else begin
    HEAD_XYAD,head,rrs_tp,dds_tp,xxs_tp,yys_tp,/deg
  endelse
  ixx_tp = CONGRID(xxs_tp,(nxs-1)*nstep+1,nstep+1,/interp,/minus_one)
  iyy_tp = CONGRID(yys_tp,(nxs-1)*nstep+1,nstep+1,/interp,/minus_one)
  xx[0:(nxs-1)*nstep,ny-nstep-1:ny-1] = ixx_tp
  yy[0:(nxs-1)*nstep,ny-nstep-1:ny-1] = iyy_tp
endif
; Deal with top/right corner
if (nxs-1)*nstep+1 lt nx and (nys-1)*nstep+1 lt ny then begin
  ; Make a short grid in X and Y at the top-right corner
  rrs_tr = rr[nx-nstep-1:nx-1:nstep,ny-nstep-1:ny-1:nstep]
  dds_tr = dd[nx-nstep-1:nx-1:nstep,ny-nstep-1:ny-1:nstep]
  if not keyword_set(xyad) then begin
    HEAD_ADXY,head,rrs_tr,dds_tr,xxs_tr,yys_tr,/deg
  endif else begin
    HEAD_XYAD,head,rrs_tr,dds_tr,xxs_tr,yys_tr,/deg
  endelse
  ixx_tr = CONGRID(xxs_tr,nstep+1,nstep+1,/interp,/minus_one)
  iyy_tr = CONGRID(yys_tr,nstep+1,nstep+1,/interp,/minus_one)
  xx[nx-nstep-1:nx-1,ny-nstep-1:ny-1] = ixx_tr
  yy[nx-nstep-1:nx-1,ny-nstep-1:ny-1] = iyy_tr
endif

end
