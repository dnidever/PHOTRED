;+
;
; ALLFRAME_GETWEIGHTS
;
; This calculates weights, scales and sky values from DAOPHOT
; photometry files for images combination.  Used by ALLFRAME.PRO and
; other programs.
;
; You need to have run daophot, allstar, daomatch and daomaster
; already.  There need als, mch and raw files.
;
; INPUTS:
;  mchfile     The MCH filename
;  /stp        Stop at the end of the program
;
; OUTPUTS:
;  actweight   The array of weights.
;  scales      The array of scales.
;  medsky      The array of skys.
;  =raw2       The RAW photometry structure.
;
; USAGE:
;  IDL>allframe_getweights,'ccd1001.mch',weights,scales,sky,raw2=raw2
;
;
; By D.Nidever   February 2008, broken out into its own program 4/8/15
;-


pro allframe_getweights,mchfile,actweight,scales,medsky,raw=raw,stp=stp

; OUTPUTS:
;  actweight  The weight for each frame
;  scales     The scale for each frame
;  medsky     The sky value for each frame

nmch = n_elements(mchfile)
if nmch eq 0 then begin
  print,'Input MCH file'
  return
endif

; Load the MCH file
LOADMCH,mchfile,files,trans
nfiles = n_elements(files)

;-----------------------------------
; Computs Weights
; getweights.sh does this

;grep RE ${image1}_1.opt > weights_1.inp
;grep FW ${image1}_1.opt >> weights_1.inp
;grep Clipped ${image1}_1.log >> weights_1.inp

; Load the opt files
rdnoise = fltarr(nfiles)
fwhm = fltarr(nfiles)
mnsky = fltarr(nfiles)
medsky = fltarr(nfiles) 
for i=0,nfiles-1 do begin
  dir = file_dirname(mchfile)
  base = file_basename(files[i],'.als')
  optfile = dir+'/'+base+'.opt'
  logfile = dir+'/'+base+'.log'

  READCOL,optfile,name,dum,value,format='A,A,F',/silent
  name = strtrim(strupcase(name),2)


  ind_re = where(name eq 'RE',nind_re)
  rdnoise[i] = value[ind_re[0]]
  ind_fw = where(name eq 'FW',nind_fw)
  fwhm[i] = value[ind_fw[0]]

  SPAWN,'grep Clipped '+logfile,out,errout
  ;              Clipped mean and median =  187.442  187.215

  ; daophot.sh log files are clipped on Tortoise for some reason
  ;  Get mean/median sky level
  if n_elements(out) eq 1 and out[0] eq '' then begin
    print,'Getting mean/median sky levels for ',base
    if file_test('daophot.opt') eq 0 then file_copy,base+'.opt','daophot.opt'
    undefine,cmdlines
    PUSH,cmdlines,'#!/bin/sh'
    PUSH,cmdlines,'export image=${1}'
    PUSH,cmdlines,'daophot << END_DAOPHOT >> ${image}.find.log'
    PUSH,cmdlines,'OPTIONS'
    PUSH,cmdlines,'${image}.opt'
    PUSH,cmdlines,' '
    PUSH,cmdlines,'ATTACH ${image}.fits'
    PUSH,cmdlines,'FIND'
    PUSH,cmdlines,'1,1'
    PUSH,cmdlines,'${image}.find.temp'
    PUSH,cmdlines,'y'
    PUSH,cmdlines,'EXIT'
    PUSH,cmdlines,'END_DAOPHOT'
    tempscript = MKTEMP('dfind')   ; absolute filename
    WRITELINE,tempscript,cmdlines
    FILE_CHMOD,tempscript,'755'o
    FILE_DELETE,base+'.find.temp',/allow
    FILE_DELETE,base+'.find.log',/allow
    ; Run DAOPHOT/FIND
    SPAWN,tempscript+' '+base,out1,errout1

    logfile2 = base+'.find.log'
    SPAWN,'grep Clipped '+logfile2,out,errout
    ;              Clipped mean and median =  187.442  187.215

    ; Delete temporary files
    FILE_DELETE,base+'.find.temp',/allow
    FILE_DELETE,tempscript,/allow
  endif

  arr = strsplit(out[0],' ',/extract)
  mnsky[i] = float(arr[5])
  medsky[i] = float(arr[6])

  ;stop

endfor

    
;      program getweight
;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
;C     Program computes weights according to signal to noise
;C     Input in aperture photometry at FWHM, the readnoise,
;C     and sky, and the program computes the S/N
;C        S/N = I/(Pi*FWHM**2) / sqrt(RN**2 + sky + I/(Pi*FWHM**2) )
;C     This S/N is scaled for each set of stars, such that the maximum is 1
;C      then the "scaled" S/N are added up, to get a "total" S/N for the frame
;C      (dividing by n would give you the average, "scaled" S/N), lastly
;C      all of the "total" S/N for each frame are summed, and that sum is
;C      then divided into each "total" S/N to give a weight, good for use
;C      in imcombine.
;C      JCO -- 2000-ish
;C
;C      Use companion shell scripts to generate input files:
;C       getweights.sh
;C       rmweights.sh
;C      Altered fortran (I hope) to read inpfile in format created by
;C       these scripts.
;C      RLB -- 06/08/2007
;C
;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

; Load the RAW file
dir = file_dirname(mchfile)
base = file_basename(mchfile,'.mch')
rawfile = dir+'/'+base+'.raw'
LOADRAW,rawfile,raw,rawhead
  
; Making mag and err arrays
nstars = n_elements(raw)
tags = tag_names(raw)
magind = where(stregex(tags,'MAG',/boolean) eq 1,nmagind)
errind = where(stregex(tags,'ERR',/boolean) eq 1,nerrind)
mag = fltarr(nfiles,nstars)
err = fltarr(nfiles,nstars)
for i=0,nfiles-1 do mag[i,*] = raw.(magind[i])
for i=0,nfiles-1 do err[i,*] = raw.(errind[i])

; Getting the reference sources
totstars = total(mag lt 50,1)
si = reverse(sort(totstars))
;gdrefstars = si[0:(99<(nstars-1))]
gdrefstars = si[0:(49<(nstars-1))]
nrefstars = n_elements(gdrefstars)
; Getting the "good" frames
totframe = total(mag[*,gdrefstars] lt 50,2)
gdframe = where(totframe eq nrefstars,ngdframe,comp=bdframe,ncomp=nbdframe)
; no good frames, lower number of reference stars
if ngdframe eq 0 then begin
  gdrefstars = si[0:(29<(nstars-1))]
  nrefstars = n_elements(gdrefstars)
  totframe = total(mag[*,gdrefstars] lt 50,2)
  gdframe = where(totframe eq nrefstars,ngdframe,comp=bdframe,ncomp=nbdframe)
endif

; Calculate the weights
actweight = fltarr(nfiles)
scales = fltarr(nfiles)
mag2 = mag[gdframe,*] & mag2 = mag2[*,gdrefstars]
err2 = err[gdframe,*] & err2 = err2[*,gdrefstars]
allframe_calcweights,mag2,err2,fwhm[gdframe],rdnoise[gdframe],medsky[gdframe],$
         actweight1,scales1
actweight[gdframe] = actweight1
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
  allframe_calcweights,mag3,err3,fwhm[indframes],rdnoise[indframes],medsky[indframes],$
                       actweight3,scales3

  ; Scale these relative to the original ones
  actweight3a = actweight[igdframe]         ; original
  actweight3b = actweight3[0:nigdframe1-1]  ; new
  wtfrac = median(actweight3a/actweight3b)
  scales3a = scales[igdframe]               ; original
  scales3b = scales3[0:nigdframe1-1]        ; new
  sclfrac = median(scales3a/scales3b)
  new_actweight = actweight3[nigdframe1] * wtfrac
  new_scale = scales3[nigdframe1] * sclfrac
  actweight[iframe] = new_actweight
  scales[iframe] = new_scale

  ;print,iframe,new_actweight,new_scale

  BOMB:

endfor

; Normalize the weights
actweight /= total(actweight)

;; Getting only stars that are well-measured in all frames
;gd = where(max(mag,dim=1) lt 50. and max(err,dim=1) lt 0.1,ngd)
;if ngd lt 10 then $
;  gd = where(max(mag,dim=1) lt 50. and max(err,dim=1) lt 0.3,ngd)
;if ngd lt 10 then $
;  gd = where(max(mag,dim=1) lt 50. and max(err,dim=1) lt 0.5,ngd)
;if ngd lt 10 then $
;  gd = where(max(mag,dim=1) lt 50.,ngd)
;
;;if ngd eq 0 then stop,'no good stars'
;
;raw2 = raw[gd]
;mag2 = mag[*,gd]
;err2 = err[*,gd]

; Rescale SCALES so the reference frames has SCALE=1.0
scales /= scales[0]

print,'Files: ',files
print,'Weights: ',actweight
print,'Scales: ',scales
print,'Sky: ',medsky

if keyword_set(stp) then stop

end
