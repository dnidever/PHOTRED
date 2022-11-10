;+
;
; DAOPHOT_IMPREP
;
; Use the CP Instcal flux and mask files to prepare
; an image for DAOPHOT
;
; INPUTS:
;  fluxfile  The filename for the CP Instcal flux file.
;  maskfile  The filename for the CP Instcal mask file.
;  /header   Return the header only.
;
; OUTPUTS:
;  im      The DAOPHOT-prepared image array
;  meta    The header/metadata for IM.
;  =error  The error if one occurred.
;
; USAGE:
;  IDL>daophot_imprep,fluxfile,maskfile,im,meta
;
; By D. Nidever  Feb 2019
;-

pro daophot_imprep,fluxfile,maskfile,im,meta,header=header,error=error

undefine,meta
im = 0

; Not enough inputs
if n_elements(fluxfile) eq 0 or n_elements(maskfile) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - daophot_imprep,fluxfile,maskfile,im,meta,header=header,error=error'
  return
endif

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   error = !ERROR_STATE.MSG
   PHOTRED_ERRORMSG,logfile=logf
   CATCH, /CANCEL 
   return
endif

;; Create the DAOPHOT file
;;   taken from smashred_imprep_single.pro
if not keyword_set(header) then begin
  FITS_READ,fluxfile,fim,fhead,/no_abort,message=ferror
  if ferror ne '' then begin
    error = ferror
    print,error
    return
  endif
  FITS_READ,maskfile,mim,mhead,/no_abort,message=merror
  if merror ne '' then begin
    error = merror
    print,error
    return
  endif
endif else fhead=HEADFITS(fluxfile)
ccdnum = sxpar(fhead,'CCDNUM')

;; --- Prepare the header ---

;; add gain, rdnoise, saturation
meta = fhead
if strmid(meta[0],0,5) eq 'XTENS' then meta[0]='SIMPLE  =                    T / file does conform to FITS standard             '

;gain = (arr[ccd-1].gaina+arr[ccd-1].gainb)*0.5
;rdnoise = (arr[ccd-1].rdnoisea+arr[ccd-1].rdnoiseb)*0.5
;gain = sxpar(fhead,'ARAWGAIN')
gainA = sxpar(fhead,'GAINA')
gainB = sxpar(fhead,'GAINB')
gain = (gainA+gainB)*0.5
rdnoiseA = sxpar(fhead,'RDNOISEA')
rdnoiseB = sxpar(fhead,'RDNOISEB')
rdnoise = (rdnoiseA+rdnoiseB)*0.5
sxaddpar,meta,'GAIN',gain
sxaddpar,meta,'RDNOISE',rdnoise

; REMOVE DUPLICATE KEYWORDS!!  They cause lots of annoying errors
; EXTVER, CHECKSUM, DATASUM
bd = where(strmid(meta,0,6) eq 'EXTVER',nbd)
if nbd gt 1 then remove,bd[1:*],meta
bd = where(strmid(meta,0,8) eq 'CHECKSUM',nbd)
if nbd gt 1 then remove,bd[1:*],meta
bd = where(strmid(meta,0,7) eq 'DATASUM',nbd)
if nbd gt 1 then remove,bd[1:*],meta

;; Add "COMMENT " before "BEGIN EXTENSION HEADER ---", it causes problems in daophot
bd = where(strmid(meta,0,5) eq 'BEGIN',nbd)
if nbd gt 0 then meta[bd]='COMMENT '+meta[bd]

if keyword_set(header) then return

;; --- Prepare the image ---

med1 = median(fim,dim=1)
med1slp = slope(med1)

im = fim

;; DAOPHOT cannot handle DOUBLE arrays (BITPIX=-64)
if sxpar(meta,'bitpix') eq -64 then begin
  sxaddpar,meta,'bitpix',-32
  im = float(im)
endif

; Check for differences in amp background levels
med1 = median(im[800:1023,*])
med2 = median(im[1024:1200,*])
err1 = mad(im[800:1023,*])/sqrt(n_elements(im[800:1023,*]))
err2 = mad(im[1024:1200,*])/sqrt(n_elements(im[1024:1200,*]))
err = sqrt(err1^2 + err2^2)

;; Set bad pixels to saturation value
;; --DESDM bit masks (from Gruendl):
;; BADPIX_BPM 1          /* set in bpm (hot/dead pixel/column)        */
;; BADPIX_SATURATE 2     /* saturated pixel                           */
;; BADPIX_INTERP 4
;;     /* interpolated pixel                        */
;; BADPIX_LOW     8      /* too little signal- i.e. poor read         */
;; BADPIX_CRAY   16      /* cosmic ray pixel                          */
;; BADPIX_STAR   32      /* bright star pixel                         */
;; BADPIX_TRAIL  64      /* bleed trail pixel                         */
;; BADPIX_EDGEBLEED 128  /* edge bleed pixel                          */
;; BADPIX_SSXTALK 256    /* pixel potentially effected by xtalk from super-saturated source */
;; BADPIX_EDGE   512     /* pixel flagged to exclude CCD glowing edges */
;; BADPIX_STREAK 1024    /* pixel associated with satellite (airplane/meteor) streak     */
;; BADPIX_FIX    2048    /* a bad pixel that was fixed                */
;; --CP bit masks, Pre-V3.5.0 (PLVER)
;; Bit   DQ Type  PROCTYPE
;; 1  detector bad pixel          InstCal
;; 1  detector bad pixel/no data  Resampled
;; 1  No data                     Stacked
;; 2  saturated                   InstCal/Resampled
;; 4  interpolated                InstCal/Resampled
;; 16  single exposure cosmic ray InstCal/Resampled
;; 64  bleed trail                InstCal/Resampled
;; 128  multi-exposure transient  InstCal/Resampled
;; --CP bit masks, V3.5.0 on (after ~10/28/2014), integer masks
;;  1 = bad (in static bad pixel mask)
;;  2 = no value (for stacks)
;;  3 = saturated
;;  4 = bleed mask
;;  5 = cosmic ray
;;  6 = low weight
;;  7 = diff detect
;; You can't have combinations but the precedence as in the order
;; of the list (which is also the order in which the processing
;; discovers them).  So a pixel marked as "bad" (1) won't ever be
;; flagged as "diff detect" (7) later on in the processing.
;;
;; "Turn off" the "difference image masking", clear the 8th bit
;; 128 for Pre-V3.5.0 images and set 7 values to zero for V3.5.0 or later.
if n_elements(nodiffmaskflag) eq 0 then nodiffmaskflag = 1  ; set by default
if keyword_set(nodiffmaskflag) then begin
  ;print,'Turning off the CP difference image masking flags'
  plver = sxpar(fhead,'plver',count=nplver)  ; DESDM doesn't have this
  plver = strtrim(plver,2)
  if nplver gt 0 then begin  ; CP data
    ; DES, didn't have this flag, so skip
    if strmid(plver,0,3) eq 'DES' then goto,SKIP

    ; V3.5.0 and on, Integer masks
    versnum = long(strsplit(strmid(plver,1),'.',/extract))
    if versnum[0] gt 3 or (versnum[0] eq 3 and versnum[1] ge 5) then begin
      bdpix = where(mim eq 7,nbdpix)
      if nbdpix gt 0 then mim[bdpix]=0

    ; Pre-V3.5.0, Bitmasks
    endif else begin
      bdpix = where( (mim and 2^7) eq 2^7,nbdpix)
      if nbdpix gt 0 then mim[bdpix]-=128   ; clear 128
    endelse
    ;print,strtrim(nbdpix,2),' pixels cleared of difference image mask flag'
  endif
endif
SKIP:

;; Add background back in for DES SV data
skybrite = sxpar(meta,'skybrite',count=nskybrite)
skysigma = sxpar(meta,'skybrite',count=nskysigma)
bunit = strtrim(sxpar(meta,'bunit',count=nbunit),2)
if nbunit eq 0 then bunit='adu'
if bunit eq 'electrons' and nskybrite gt 0 then begin
  saturate = sxpar(meta,'saturate',count=nsaturate)
  if nsaturate eq 0 then saturate=65000.0
  gdpix = where(mim eq 0 and im lt saturate,ngdpix)
  medim = median(im[gdpix])
  if medim lt 100 then begin
    ;; Add sky background back in so DAOPHOT can create a proper
    ;; noise model
    ;; SKYBRITE is given in electrons as well
    if skybrite ne 0 then begin
      im[gdpix] += skybrite
    endif else begin
      ;; Sometimes SKYBRITE=0, in that case use SKYSIGMA^2
      ;;  measure sigma ourselves if SKYSIGMA is not given
      if nskysigma gt 0 and skysigma gt 0 then im[gdpix]+=skysigma^2 else im[gdpix]+=mad(im[gdpix])^2
    endelse
  endif
endif

;; Mask bad half of DECam chip 31
dateobs = sxpar(meta,'DATE-OBS',count=ndateobs)
if ndateobs gt 0 then mjd=date2jd(dateobs,/mjd) else mjd=58000.0d0  ;; assume bad if no date
if strtrim(sxpar(meta,'INSTRUME'),2) eq 'DECam' and sxpar(meta,'CCDNUM') eq 31 and mjd gt 56660 then begin
  print,'Masking bad half of DECam chip 31'
  ;; X: 1-1000 okay
  ;; X: 1000-2049 bad
  fim[1000:*,*] = 1e6
endif

;; Set saturated pixels to 65000.0
if bunit ne 'electrons' then begin
  saturate = sxpar(meta,'saturate',count=nsaturate)
  if nsaturate gt 0 then saturate<=64500.0 else saturate=64500.0  ; set it slightly lower than 65000 for DAOPHOT
  sxaddpar,meta,'saturate',saturate
  bdpix = where(mim gt 0.0 or fim gt 65000.0,nbdpix)
  if nbdpix gt 0 then im[bdpix]=65000.0
;; allow saturation value to be larger for DES SV data in electrons
endif else begin
  saturate = sxpar(meta,'saturate',count=nsaturate)
  if nsaturate eq 0 then begin
    saturate = 64500.0  ; set it slightly lower than 65000 for DAOPHOT
    bdpix = where(mim gt 0.0 or fim gt 65000.0,nbdpix)
    if nbdpix gt 0 then im[bdpix]=65000.0
    sxaddpar,meta,'saturate',saturate
  endif else begin
    bdpix = where(mim gt 0.0 or fim gt saturate,nbdpix)
    if nbdpix gt 0 then im[bdpix]=saturate*1.01  ; set slightly higher for DAOPHOT
    sxaddpar,meta,'saturate',saturate
  endelse
endelse

end
