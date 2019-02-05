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
;
; OUTPUTS:
;  im      The DAOPHOT-prepared image array
;  head    The header for IM.
;  =error  The error if one occurred.
;
; USAGE:
;  IDL>daophot_imprep,fluxfile,maskfile,im,head
;
; By D. Nidever  Feb 2019
;-

pro daophot_imprep,fluxfile,maskfile,im,head,error=error

; Not enough inputs
if n_elements(fluxfile) eq 0 or n_elements(maskfile) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - daophot_imprep,fluxfile,maskfile,im,head,error=error'
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
   return,-1
endif

;; Create the DAOPHOT file
;;   taken from smashred_imprep_single.pro
FITS_READ,fluxfile,fim,fhead
FITS_READ,maskfile,mim,mhead
ccdnum = sxpar(fhead,'CCDNUM')

; Need to add gain, rdnoise, saturation
med1 = median(fim,dim=1)
med1slp = slope(med1)

newim = fim

; Check for differences in amp background levels
med1 = median(newim[800:1023,*])
med2 = median(newim[1024:1200,*])
err1 = mad(newim[800:1023,*])/sqrt(n_elements(newim[800:1023,*]))
err2 = mad(newim[1024:1200,*])/sqrt(n_elements(newim[1024:1200,*]))
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

bdpix = where(mim gt 0.0,nbdpix)
if nbdpix gt 0 then newim[bdpix]=6e4

;; add gain, rdnoise, saturation
newhead = fhead
if strmid(newhead[0],0,5) eq 'XTENS' then newhead[0]='SIMPLE  =                    T / Fits standard'

;gain = (arr[ccd-1].gaina+arr[ccd-1].gainb)*0.5
;rdnoise = (arr[ccd-1].rdnoisea+arr[ccd-1].rdnoiseb)*0.5
;gain = sxpar(fhead,'ARAWGAIN')
gainA = sxpar(fhead,'GAINA')
gainB = sxpar(fhead,'GAINB')
gain = (gainA+gainB)*0.5
rdnoiseA = sxpar(fhead,'RDNOISEA')
rdnoiseB = sxpar(fhead,'RDNOISEB')
rdnoise = (rdnoiseA+rdnoiseB)*0.5
sxaddpar,newhead,'GAIN',gain
sxaddpar,newhead,'RDNOISE',rdnoise

; REMOVE DUPLICATE KEYWORDS!!  They cause lots of annoying errors
; EXTVER, CHECKSUM, DATASUM
bd = where(strmid(newhead,0,6) eq 'EXTVER',nbd)
if nbd gt 1 then remove,bd[1:*],newhead
bd = where(strmid(newhead,0,8) eq 'CHECKSUM',nbd)
if nbd gt 1 then remove,bd[1:*],newhead
bd = where(strmid(newhead,0,7) eq 'DATASUM',nbd)
if nbd gt 1 then remove,bd[1:*],newhead

;; Add "COMMENT " before "BEGIN EXTENSION HEADER ---", it causes problems in daophot
bd = where(strmid(newhead,0,5) eq 'BEGIN',nbd)
if nbd gt 0 then newhead[bd]='COMMENT '+newhead[bd]

;;; Put in FPACK parameters
;if keyword_set(fpack) then begin
;  ; Remove all past FZ parameters
;  bd = where(strmid(newhead,0,2) eq 'FZ',nbd)
;  if nbd gt 0 then REMOVE,bd,newhead
;  sxaddpar,newhead,'FZALGOR','RICE_1'
;  sxaddpar,newhead,'FZQMETHD','SUBTRACTIVE_DITHER_1'
;  sxaddpar,newhead,'FZQVALUE',8
;  sxaddpar,newhead,'FZDTHRSD','CHECKSUM'
;endif

end
