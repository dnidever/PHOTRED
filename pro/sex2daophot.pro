;+
;
; SEX2DAOPHOT
;
; This program converts a SExtractor output file
; to DAOPHOT format.  Currently only coo.
;
; INPUTS:
;  catfile   The SExtractor catalog filename.
;  fitsfile  The name of the associated FITS file.
;  daofile   The name of the output DAOPHOT file.
;
; OUTPUTS:
;  The catalog is written to DAOFILE.
;  =error    The error if one occurred.
;
; USAGE:
;  IDL>sex2daophot,'F1-12340056_01.cat','F1-12340056_01.fits','F1-12340056_01.coo'
;
; By D. Nidever  Jan 2019
;-

pro sex2daophot,catfile,fitsfile,daofile,error=errror

;; Not enough inputs
if n_elements(catfile) eq 0 or n_elements(fitsfile) eq 0 or n_elements(daofile) eq 0 then begin
  print,'Syntax - sex2daophot,catfile,fitsfile,daofile,error=error'
  error = 'Not enough inputs'
  return
endif

;; Check that the needed files exist
if file_test(catfile) eq 0 then begin
  error = catfile+' NOT FOUND'
  print,error
  return
endif
if file_test(fitsfile) eq 0 then begin
  error = fitsfile+' NOT FOUND'
  print,error
  return
endif

;-------------------------------------
; Load sextractor output file
; default.param specifies the output columns
if file_isfits(catfile) eq 0 then begin
  fields = ['ID','X','Y','MAG','ERR','FLAGS','STAR']
  sex = IMPORTASCII(catfile,fieldnames=fields,/noprint)
endif else sex=MRDFITS(catfile,1,/silent)
nsex = n_elements(sex)

;-------------------------------------
; Get meta-data from the FITS file
head = headfits(fitsfile)
naxis1 = sxpar(head,'NAXIS1')
naxis2 = sxpar(head,'NAXIS2')
saturate = sxpar(head,'SATURATE')
rdnoise = PHOTRED_GETRDNOISE(fitsfile)
gain = PHOTRED_GETGAIN(fitsfile)
lowbad = 1.0
thresh = 20.0

; Header values:  this information comes from daophot2.pdf pg.69
; NL: Originally meant "number of lines" but not anymore
; NX: size of X-dimension of image in pixels
; NY: size of Y-dimension of image in pixels
; LOWBAD: lower good data limit, calculated by FIND
; HIGHBAD: upper good data limit, specified in option file
; THRESH: threshold calculated by FIND
; AP1: radius (pixels) of the first aperture used by PHOTOMETRY
; PH/ADU: gain in photons/ADU used when running FIND
; RDNOISE: rdnoise (ADU) used when running FIND
; FRAD: value of fitting radius


;; Making ALS structure for new SEX sources
schema = {ID:0L,X:0.0,Y:0.0,MAG:0.0,ERR:0.0,SKY:0.0,ITER:0.0,CHI:0.0,SHARP:0.0}
dao = replicate(schema,nsex)
dao.id = sex.id
dao.x = sex.x
dao.y = sex.y
dao.mag = sex.mag
dao.err = sex.err
dao.sky = 0.0
dao.iter = 1
dao.chi = 1.0
dao.sharp = 0.0

;;NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
;  1  2046  4094  1472.8 38652.0   80.94    0.00    3.91    1.55    3.90
; 
;      1  1434.67    15.59   -0.045    0.313    0.873    1.218
;      2   233.85    18.42   -0.018    0.218   -0.781    1.433
;    ID      X         Y       MAG     SHARP    ROUND    ROUND2
openw,unit,daofile,/get_lun
printf,unit," NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD"
printf,unit,format='(I3,I6,I6,F8.1,F8.1,F8.2,F8.2,F8.2,F8.2,F8.2)',1,naxis1,naxis2,lowbad,$
       saturate,thresh,3.0,gain,rdnoise/gain,3.9
printf,unit,''
; Write the data
for i=0,nsex-1 do $
   printf,unit,format='(I7,2F9.2,4F9.3)',dao[i].id,dao[i].x,dao[i].y,dao[i].mag,0.6,0.0,0.0
close,unit
free_lun,unit

end
