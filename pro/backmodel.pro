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
;             to use.  The default is 1.
;  =nsigrej The sigma clipping value to use for outlier
;             rejection. The default is 2.5.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  model    The 2D background model
;  =medstr  The structure with the median values in the spatial bins.
;  =maskim  The masked image.  Masked pixels are set to NAN.
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
medstr = replicate({x:0.0,y:0.0,med:0.0,sig:0.0,mode:0.0,value:0.0,$
                    xlo:0L,xhi:0L,ylo:0L,yhi:0L},nbins,nbins)
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
    med = median([subim],/even)
    if finite(med) eq 1 then begin
      medstr[i,j].med = med
      sig = mad([subim-med],/zero) > 0.1
      medstr[i,j].sig = sig
      ; Use histogram to get Mode
      bin1 = sig/3.
      hist1 = histogram(subim,bin=bin1,min=med-2.5*sig,max=med+2.5*sig,locations=xhist1)
      xhist1 += 0.5*bin1
      maxind1 = first_el(maxloc(hist1))
      mode1 = xhist1[maxind1]
      ; A second time around the mode to refine the value
      bin2 = sig/6.
      hist = histogram(subim,bin=bin2,min=mode1-1.5*sig,max=med+1.5*sig,locations=xhist)
      xhist += 0.5*bin2
      maxind = first_el(maxloc(hist))
      mode = xhist[maxind]
    endif else mode=!values.f_nan
    medstr[i,j].mode = mode
    medstr[i,j].value = mode
    ;MMM,subim,skymod,sigma,skew
    ;medstr[i,j].mode = skymod
    ;medstr[i,j].sigma = sigma
    ;medstr[i,j].skew = skew
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


pro backmodel,im,model,medstr=medstr,nbins=nbins,ngrow=ngrow,maxiter=maxiter,$
              nsigrej=nsigrej,maskim=maskim,stp=stp

undefine,model,medstr

; Not enough inputs
if n_elements(im) eq 0 then begin
  print,'Syntax - backmodel,im,model,medstr=medstr,nbins=nbins,ngrow=ngrow,nsigrej=nsigrej,maskim=maskim,stp=stp'
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
if n_elements(maxiter) eq 0 then maxiter=1
if n_elements(nsigrej) eq 0 then nsigrej=2.5

; Initial rough masking
maskim = im
med = median(im)
sig = mad(im-med,/zero)
bdpix = where( abs(im-med) gt 3*sig,nbdpix)
if nbdpix gt 0 then maskim[bdpix] = !values.f_nan

; Determine medians in large bins
BACKMODEL_MEDIM,maskim,medstr,nbins=nbins

; 2D linear interpolation
model = CONGRID(medstr.value,nx,ny,/interp,/center)

; Outlier rejection loop
;maskim = im
lastmodel = model
threshold = 0.01*median(medstr.value)  ; 1% of the background
niter = 0
flag = 0
WHILE (flag eq 0) do begin

  ;; Determine medians in large bins
  ;BACKMODEL_MEDIM,maskim,medstr1,nbins=nbins
  ;
  ;; 2D linear interpolation
  ;model1 = CONGRID(medstr1.value,nx,ny,/interp,/center)

  ; Get the residuals and their scatter
  BACKMODEL_MEDIM,maskim-model,medstr2,nbins=nbins
  ; use the sigma of the residuals
  sig = medstr.sig < medstr2.sig  ; use the lower of the two values
  sig1 = CONGRID(sig,nx,ny,/interp,/center)

  ; Remove outliers and fit again
  ;  measure scatter with a sample of the pixels
  ;nsample = 1000 > round(0.05*nx*ny) < (nx*ny)
  ;rnd = sort(randomu(seed,nsample))
  ;sig = mad((im-model1)(rnd))
  maskim = im
  bdpix = where( abs(im-model) gt nsigrej*sig1,nbdpix)
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
    maskim[bdpix] = !values.f_nan

    ; Determine medians in large bins again
    BACKMODEL_MEDIM,maskim,medstr,nbins=nbins

    ; Final 2D linear interpolation
    model = CONGRID(medstr.value,nx,ny,/interp,/center)
  endif

  ; Ready to stop?
  diff = lastmodel-model
  if nbdpix eq 0 or max(abs(diff)) lt threshold or niter ge maxiter then flag=1

  lastmodel = model
  niter++
ENDWHILE

;stop

if keyword_set(stp) then stop

end
