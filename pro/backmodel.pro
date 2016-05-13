;+
;
; BACKMODEL
;
; Determine 2D backgroun model for an image.
;
; INPUTS:
;  im       The 2D image to determine the background for.
;  =nbins   Number of spatial bins to use in each dimension.
;             The default is 12.
;  =ngrow   The number of pixels to grow the bad pixel mask.
;             The default is 2.
;  =maxiter The maximum number of outer rejection iterations
;             to use.  The default is 2.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  model    The 2D background model
;  =medstr  The structure with the median values in the spatial bins.
;
; USAGE:
;  IDL>backmodel,im,model
;
; By D.Nidever  May 2016
;-

pro backmodel_medim,im,medstr,nbins=nbins

if n_elements(nbins) eq 0 then nbins=12

sz = size(im)
nx = sz[1]
ny = sz[2]
xstep = round(nx/nbins)
ystep = round(ny/nbins)
medim = fltarr(nbins,nbins)
medstr = replicate({x:0.0,y:0.0,value:0.0,sig:0.0,xlo:0L,xhi:0L,ylo:0L,yhi:0L},nbins,nbins)
;medstr = replicate({x:0.0,y:0.0,value:0.0,xlo:0L,xhi:0L,ylo:0L,yhi:0L},nbins,nbins)

for i=0,nbins-1 do begin
  for j=0,nbins-1 do begin
    xlo = i*xstep
    xhi = (xlo+xstep-1) < (nx-1)
    ylo = j*ystep
    yhi = (ylo+ystep-1) < (ny-1)
    medstr[i,j].xlo = xlo
    medstr[i,j].xhi = xhi
    medstr[i,j].ylo = ylo
    medstr[i,j].yhi = yhi
    medstr[i,j].x = mean([xlo,xhi])
    medstr[i,j].y = mean([ylo,yhi])
    subim = im[xlo:xhi,ylo:yhi]
    medstr[i,j].value = median([subim],/even)
    medstr[i,j].sig = mad([subim])
  endfor
endfor

end

;----------

pro backmodel_fit2d,medstr,pars,npars=npars

if n_elements(npars) eq 0 then npars=8

; 2D polynomial fit
func = 'func_poly2d'
; 6:  a = p[0] + p[1]*x + p[2]*x^2 + p[3]*x*y + p[4]*y + p[5]*y^2
initpars = fltarr(npars)  ; 6, 11
initpars[0] = median(medstr.value)
;coef1 = robust_poly_fitq(medstr.x,medstr.value,1)
;initpars[1] = coef1[1]
;coef2 = robust_poly_fitq(medstr.y,medstr.value,1)
;initpars[4] = coef2[1]
pars = MPFIT2DFUN(func,medstr.x,medstr.y,medstr.value,medstr.sig,initpars,yfit=yfit,status=status,/quiet)

end

;----------


pro backmodel,im,model,medstr=medstr,nbins=nbins,ngrow=ngrow,maxiter=maxiter,stp=stp

undefine,model,medstr

; Not enough inputs
if n_elements(im) eq 0 then begin
  print,'Syntax - backmodel,im,model,medstr=medstr,nbins=nbins,ngrow=ngrow,stp=stp'
  return
endif

sz = size(im)
if sz[0] ne 2 then begin
  print,'Error.  Image must have 2 dimensions'
  return
endif
nx = sz[1]
ny = sz[2]

; Defaults
if n_elements(nbins) eq 0 then nbins=12
if n_elements(ngrow) eq 0 then ngrow=2
if n_elements(maxiter) eq 0 then maxiter=2

; Outlier rejection loop
tempim = im
lastmodel = im*0
niter = 0
flag = 0
WHILE (flag eq 0) do begin

  ; Determine medians in large bins
  BACKMODEL_MEDIM,tempim,medstr1,nbins=nbins

  ; 2D linear interpolation
  model1 = CONGRID(medstr1.value,nx,ny,/interp,/center)

  ; Get the residuals and their scatter
  BACKMODEL_MEDIM,tempim-model1,medstr2,nbins=nbins
  ; use the sigma of the residuals
  sig = medstr1.sig < medstr2.sig  ; use the lower of the two values
  sig1 = CONGRID(sig,nx,ny,/interp,/center)

  ; Remove outliers and fit again
  ;  measure scatter with a sample of the pixels
  ;nsample = 1000 > round(0.05*nx*ny) < (nx*ny)
  ;rnd = sort(randomu(seed,nsample))
  ;sig = mad((im-model1)(rnd))
  tempim = im
  bdpix = where( abs(im-model1) gt 2.5*sig1,nbdpix)
  if nbdpix gt 0 then begin
    ; Grow the bad pixel regions
    if ngrow gt 0 then begin
      mask = bytarr(nx,ny)
      mask[bdpix] = 1
      kernel = lonarr(2*ngrow+1,2*ngrow+1)+1
      mask = convol(mask,kernel)
      mask = mask/(mask>1)
      bdpix = where(mask eq 1,nbdpix)
    endif
    tempim[bdpix] = !values.f_nan

    ; Determine medians in large bins again
    BACKMODEL_MEDIM,tempim,medstr,nbins=nbins

    ; Final 2D linear interpolation
    model = CONGRID(medstr.value,nx,ny,/interp,/center)

  ; Nothing to mask
  endif else begin
    model = model1
    medstr = medstr1
  endelse

  ; Ready to stop?
  diff = lastmodel-model
  threshold = 50
  if nbdpix eq 0 or max(abs(diff)) lt threshold or niter ge maxiter then flag=1

  lastmodel = model
  niter++
ENDWHILE

if keyword_set(stp) then stop

end
