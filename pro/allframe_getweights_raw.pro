;+
;
; ALLFRAME_GETWEIGHTS_RAW
;
; This calculates weights, scales and sky values from the arrays from
; the DAOPHOT .raw file.  Called by allframe_getweights.pro
;
; INPUTS:
;  str         Input structure (one element per file) giving
;               MAG[Nstars], ERR[Nstars], FWHM, RDNOISE, MEDSKY
;
; OUTPUTS:
;  outstr      Output structure similar to input structure but
;                with WEIGHTS and SCALES added.
;
; USAGE:
;  IDL>allframe_getweights_raw,str,outstr
;
; By D.Nidever   Jan 2019, broke out functionality from allframe_getweights.pro
;-

pro allframe_getweights_raw,str,outstr,stp=stp

undefine,outstr
nstr = n_elements(str)
if nstr eq 0 then begin
  print,'Syntax - allframe_getweights_raw,str,outstr,stp=stp'
  return
endif

mag = transpose(str.mag)
err = transpose(str.err)
fwhm = str.fwhm
rdnoise = str.rdnoise
medsky = str.medsky

if size(mag,/n_dim) eq 2 then nstars = n_elements(mag[0,*]) else nstars=n_elements(mag)
if size(mag,/n_dim) eq 2 then nfiles = n_elements(mag[*,0]) else nfiles=1

; Getting the reference sources
totstars = total(mag lt 50,1)
si = reverse(sort(totstars))   ; get the stars with the most detections
;gdrefstars = si[0:(99<(nstars-1))]
gdrefstars = si[0:(49<(nstars-1))]
nrefstars = n_elements(gdrefstars)
; Getting the "good" frames
totframe = total(mag[*,gdrefstars] lt 50,2)  ; # of ref stars detected per frame
gdframe = where(totframe eq nrefstars,ngdframe,comp=bdframe,ncomp=nbdframe)
; No good frames, lower number of reference stars
if ngdframe eq 0 then begin
  gdrefstars = si[0:(29<(nstars-1))]
  nrefstars = n_elements(gdrefstars)
  totframe = total(mag[*,gdrefstars] lt 50,2)
  gdframe = where(totframe eq nrefstars,ngdframe,comp=bdframe,ncomp=nbdframe)
endif
; No good frames again, say the frame with the most detections is "good"
;   get weights relative to that one for the others
if ngdframe eq 0 then begin
  ; say the frame with the most detections is "good" and
  ;  the rest are bad, get weights relative to this one frame
  totstarsframe = total(mag lt 50,2)
  gdframe = first_el(maxloc(totstarsframe))
  bdframe = lindgen(nfiles,1)
  remove,gdframe,bdframe
  nbdframe = n_elements(bdframe)
  ; Get stars that are good in this frames and in ALL of the others
  gdrefstars = where(reform(mag[gdframe,*]) lt 50 and totstars eq nfiles,nrefstars)
  ;  lower threshold, at least half
  if nrefstars eq 0 then gdrefstars = where(reform(mag[gdframe,*]) lt 50 and totstars gt 0.5*nfiles,nrefstars)
  ;  just the good ones
  if nrefstars eq 0 then begin
    gdrefstars = where(reform(mag[gdframe,*]) lt 50,nrefstars)
    si = reverse(sort(totstars[gdrefstars]))           ; order by how many other frames they are detected in
    gdrefstars = gdrefstars[si[0:(49<(nrefstars-1))]]  ; only want 50
    nrefstars = n_elements(gdrefstars)
  endif
endif

; Calculate the weights
weights = fltarr(nfiles)
scales = fltarr(nfiles)
mag2 = mag[gdframe,*] & mag2 = mag2[*,gdrefstars]
err2 = err[gdframe,*] & err2 = err2[*,gdrefstars]
ALLFRAME_CALCWEIGHTS,mag2,err2,fwhm[gdframe],rdnoise[gdframe],medsky[gdframe],$
         weights1,scales1
weights[gdframe] = weights1
scales[gdframe] = scales1

; If there are some "bad" frames calculate their weights
;  relative to some of the "good" ones
for i=0,nbdframe-1 do begin

  iframe = bdframe[i]

  ; First round of "good" stars
  totgdstars = total(mag[gdframe,*] lt 50,1)  ; stars in good frames
  igdrefstars1 = where(mag[iframe,*] lt 50 and totgdstars ge 5,nigdrefstars1)
  if nigdrefstars1 lt 3 then $
    igdrefstars1 = where(mag[iframe,*] lt 50 and totgdstars ge 1,nigdrefstars1)
  if nigdrefstars1 lt 2 then goto,BOMB

  totgdstars1 = total( (mag[gdframe,*])[*,igdrefstars1] lt 50,1)
  si1 = reverse(sort(totgdstars1))  ; get the best ones
  igdrefstars = igdrefstars1[si1[0:(49<(nigdrefstars1-1))]]
  nirefstars = n_elements(igdrefstars)

  itotframe = total( (mag[gdframe,*])[*,igdrefstars] lt 50,2)
  igdframe1 = where( itotframe eq nirefstars,nigdframe1)

  ; Whittle down to the best stars/frames
  if nigdframe1 eq 0 then begin

    igdframe = gdframe
    igdrefstars = igdrefstars

    ; whittle down to best stars and frames
    count = 0
    endflag = 0
    while endflag eq 0 do begin

      ; sum over frames
      itot1 = total( (mag[igdframe,*])[*,igdrefstars] lt 50,1)
      si1 = sort(itot1)
      ; remove worst 5% stars
      if n_elements(igdrefstars) gt 3 then begin
        bd1 = si1[0:round(n_elements(igdrefstars)*0.05)]
        remove,bd1,igdrefstars
      endif

      ; sum over stars
      itot2 = total( (mag[igdframe,*])[*,igdrefstars] lt 50,2)
      si2 = sort(itot2)
      ; remove worst 5% frames
      if n_elements(igdframe) gt 1 then begin
        bd2 = si2[0:round(n_elements(igdframe)*0.05)]
        remove,bd2,igdframe
      endif     

      ; Testing if we should end
      itotframe = total( (mag[igdframe,*])[*,igdrefstars] lt 50,2)
      igdframe1 = where( itotframe eq n_elements(igdrefstars),nigdframe1)
      if nigdframe1 gt 0 then endflag=1

      if endflag eq 1 then igdframe=igdframe[igdframe1]
      if count gt 50 then goto,BOMB  ; not converging, end it

      count++
    endwhile

  endif else igdframe=gdframe[igdframe1]

  ; Get weights relative to some "good" frames
  indframes = [igdframe,iframe]
  mag3 = (mag[indframes,*])[*,igdrefstars]
  err3 = (err[indframes,*])[*,igdrefstars]
  ALLFRAME_CALCWEIGHTS,mag3,err3,fwhm[indframes],rdnoise[indframes],medsky[indframes],$
                       weights3,scales3

  ; Scale these relative to the original ones
  weights3a = weights[igdframe]         ; original
  weights3b = weights3[0:nigdframe1-1]  ; new
  wtfrac = median(weights3a/weights3b)
  scales3a = scales[igdframe]               ; original
  scales3b = scales3[0:nigdframe1-1]        ; new
  sclfrac = median(scales3a/scales3b)
  new_weights = weights3[nigdframe1] * wtfrac
  new_scale = scales3[nigdframe1] * sclfrac
  weights[iframe] = new_weights
  scales[iframe] = new_scale

  ;print,iframe,new_weights,new_scale

  BOMB:

endfor

;; Fix images with bad weights likely due to no overlap
;;  use the FLUX values to get a weight
bdweights = where(weights le 0.0,nbdweights,comp=gdweights,ncomp=ngdweights)
if nbdweights gt 0 then begin
  ;; Use fluxrate10 to get weights and flux10 to get scales
  ;; relative to "good" frame
  if ngdweights gt 0 then begin
    weights[bdweights] = str[bdweights].fluxrate10 * median([weights[gdweights]/str[gdweights].fluxrate10])
    scales[bdweights] = str[bdweights].flux10 * median([scales[gdweights]/str[gdweights].flux10])
  ;; all bad
  endif else begin
    weights = str.fluxrate10
    scales = str.flux10
  endelse
endif

; Normalize the weights
weights /= total(weights)

; Rescale SCALES so the reference frames has SCALE=1.0
if scales[0] gt 0.0 then scales /= scales[0] else scales/=max(scales)

;; Create the output structure
schema = str[0]
struct_assign,{dum:''},schema
if tag_exist(schema,'weight') eq 0 then schema = create_struct(schema,'weight',0.0)
if tag_exist(schema,'scale') eq 0 then schema = create_struct(schema,'scale',0.0)
outstr = replicate(schema,nfiles)
struct_assign,str,outstr
outstr.weight = weights
outstr.scale = scales

if keyword_set(stp) then stop

end
