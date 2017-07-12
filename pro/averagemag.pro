;+
;
; AVERAGEMAG
;
; This program averages photometry using a flux
; weighting method.
;
; INPUTS:
;  mag       2D Array of magnitudes, [Nstar,Nobs]
;  err       2D array of magnitude uncertainties
;              that goes with MAG.
;  /robust   Reject outliers.
; 
; OUTPUTS:
;  newmag    The flux-weighted magnitudes [Nstar].
;  newerr    The uncertainties of NEWMAG.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>averagemag,mag,err,newmag,newerr,/robust
;
; By D. Nidever    break out of photcalib.pro  July 2017
;-

PRO averagemag,mag,err,newmag,newerr,robust=robust,error=error

; Average magnitudes
; /robust  Outlier rejection

undefine,newmag
undefine,newerr
undefine,error

; Not enough inputs
if n_elements(mag) eq 0 or n_elements(err) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - averagemag,mag,err,newmag,newerr,robust=robust,error=error'
  return
endif
  
sz = size(mag)
szerr = size(err)

; Check the sizes
if sz[0] ne 2 or szerr[0] ne 2 then begin
  error = 'MAG and ERR must be 2D'
  print,error
  return
endif
if sz[1] ne szerr[1] or sz[2] ne szerr[2] then begin
  error = 'MAG and ERR must be the same size'
  print,error
  return
endif

numstar = sz[1]
numobs = sz[2]

; Starting arrays
finalflux = fltarr(numstar)+99.9999
newmag = fltarr(numstar)+99.9999
newerr = fltarr(numstar)+9.9999

tmag = mag
terr = err
tflux = 2.5118864d0^tmag
twt = 1.0/(terr^2.0)
bd = where(tmag gt 50,nbd)
if nbd gt 0 then begin
  tmag[bd] = !values.f_nan
  terr[bd] = !values.f_nan
  tflux[bd] = !values.f_nan
  twt[bd] = !values.f_nan
endif

; Outlier rejection
if numobs gt 2 and keyword_set(robust) then begin

  ; Number of detections per star
  ngood = total(finite(tmag),2)
  ngood2 = ngood#replicate(1,numobs)

  ; Get average sigma value and median difference
  medmag = median(tmag,dim=2)  ; no /even, want an actual value
  diffmag = tmag-medmag#replicate(1,numobs)
  sigmag = mad(diffmag,dim=2,/zero)
  gdsig = where(finite(sigmag) eq 1,ngdsig,comp=bdsig,ncomp=nbdsig)
  if ngdsig gt 1 then maxsig=max(sigmag[gdsig]) else maxsig=0.5
  if nbdsig gt 0 then sigmag[bdsig]=maxsig

  ; Get average error
  avgerr = total(terr,2,/nan)/(ngood>1)
  gderr = where(finite(avgerr) eq 1,ngderr,comp=bderr,ncomp=nbderr)
  if ngderr gt 1 then maxerr=max(avgerr[gderr]) else maxerr=0.5
  if nbderr gt 0 then avgerr[bderr]=maxerr

  ; Need to set threshold on a star-by-star basis since the scatter
  ; will vary with magnitude
  
  ; Old method
  ;gdsig = where(finite(sigmag) eq 1 and sigmag gt 0.0,ngdsig,comp=bdsig,ncomp=nbdsig)
  ;if ngdsig gt 0 then avgsig=median([sigmag[gdsig]]) else avgsig=0.2  ; average sigma value

  ; Rejected "bad" values, only for stars detection in 3 or more frames
  ;bdval = where(abs(diffmag) gt 3*(0.2>avgsig) or finite(tmag) eq 0 and ngood2 ge 3,nbdval)
  threshold = 5.0 * (sigmag > avgerr > 0.05)
  threshold2 = threshold#replicate(1,numobs)
  bdval = where(abs(diffmag) gt threshold2 or finite(tmag) eq 0 and ngood2 ge 3,nbdval)
  if nbdval gt 0 then begin
    tflux[bdval] = !values.f_nan
    twt[bdval] = !values.f_nan
  endif

endif ; outlier rejection

; Now do flux-weighted average for the leftover "good" values
totalwt = total(twt,2,/nan)
totalfluxwt = total(twt*tflux,2,/nan)
bdwt = where(finite(totalwt) eq 0,nbdwt)
if nbdwt gt 0 then totwt[bdwt]=0.0

; Calculate final mags and errors
fgd = where(totalwt gt 0.0,nfgd)  ; stars with at least one good mag
if (nfgd gt 0) then begin
  finalflux[fgd] = totalfluxwt[fgd]/totalwt[fgd]
  newmag[fgd] = 2.50*alog10(finalflux[fgd])
  newerr[fgd] = sqrt(1.0/totalwt[fgd])
endif

; OLD METHOD

; ; Starting arrays
; colmag = fltarr(numstar)+99.9999
; colerr = fltarr(numstar)+9.9999
; totalwt = fltarr(numstar)
; totalfluxwt = fltarr(numstar)
; finalflux = fltarr(numstar)
;
; ; The equations are:
; ; wt = 1.0/(err^2.0)
; ; totalwt = (wt1 + wt2 + ...)
; ; finalflux = (flux1*wt1 + flux2*wt2 + ...)/totalwt
; ; finalmag = 2.50*alog10(finalflux)
; ; finalerr = sqrt(1.0/totalwt)
; ;
; ; Add totalwt and (flux1*wt1 + ...) as you go only for good mags
; ; Then calculate finalmag and finalerr at the end
;
; ; Looping through the multiple exposures
; for k=0,nind-1 do begin
;
;   mag = tempmag[*,ind[k]]
;   ;err = temperr[*,ind[k]]
;   err = inerr[*,ind[k]]            ; use the instrumental errors
;   gd = where(mag lt 50.0,ngd)
;
;   flux = mag*0.0
;   wt = mag*0.0
;   if (ngd gt 0) then begin
;     flux[gd] = 2.511864^mag[gd]
;     wt[gd] = 1.0/(err[gd]^2.0)
;     totalwt[gd] = totalwt[gd] + wt[gd]
;     totalfluxwt[gd] = totalfluxwt[gd] + flux[gd]*wt[gd]
;   endif
; end ; multiple bands loop
;
; ; Calculate final mags and errors
; fgd = where(totalwt gt 0.0,nfgd)  ; stars with at least one good mag
; if (nfgd gt 0) then begin
;   finalflux[fgd] = totalfluxwt[fgd]/totalwt[fgd]
;   colmag[fgd] = 2.50*alog10(finalflux[fgd])
;   colerr[fgd] = sqrt(1.0/totalwt[fgd])
; endif

end
