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
;  =imager     Imager structure with basic information
;  =logfile    A logfile to print to output to.
;  /silent     Don't print anything to the screen.
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
; By D.Nidever   February 2008, broken out into its own program 4/8/15
;-


pro allframe_getweights,mchfile,actweight,scales,medsky,imager=imager,logfile=logfile,raw=raw,silent=silent,stp=stp

COMMON photred,setupOA

; OUTPUTS:
;  actweight  The weight for each frame
;  scales     The scale for each frame
;  medsky     The sky value for each frame

tilesep = '+'
;tilesep = '.'
btilesep = long(byte(tilesep))

nmch = n_elements(mchfile)
if nmch eq 0 then begin
  print,'Syntax - allframe_getweights,mchfile,actweight,scales,medsky,imager=imager,logfile=logfile,raw=raw,silent=silent,stp=stp'
  return
endif

; MCH file not found
if file_test(mchfile) eq 0 then begin
  print,mchfile,' NOT FOUND'
  return
endif

if n_elements(logfile) eq 0 then logfile=-1

; Load the MCH file
LOADMCH,mchfile,files,trans
nfiles = n_elements(files)

;-----------------------------------
; Get the information that we need

; Load the opt files
info = replicate({name:'',exists:0,filter:'',exptime:0.0,fwhm:0.0,rdnoise:0.0,$
                  mnsky:0.0,medsky:0.0,mag10:99.99,flux10:0.0,fluxrate10:0.0,$
                  weight:0.0,scale:0.0},nfiles)
info.name = files
for i=0,nfiles-1 do begin
  dir = file_dirname(mchfile)
  base = file_basename(files[i],'.als')
  optfile = dir+'/'+base+'.opt'
  logfile1 = dir+'/'+base+'.log'
  
  ;; Check that the FITS file exists
  info[i].exists = 1
  fitsfile = base+'.fits'
  if file_test(fitsfile) eq 0 then fitsfile+='.fz'
  if file_test(fitsfile) eq 1 then begin
    info1 = file_info(fitsfile)
    if info1.size le 1 then info[i].exists = 0
  endif else info[i].exists = 0
  if info[i].exists eq 0 then print,fitsfile+' NOT FOUND'

  READCOL,optfile,name,dum,value,format='A,A,F',/silent
  name = strtrim(strupcase(name),2)

  ind_re = where(name eq 'RE',nind_re)
  info[i].rdnoise = value[ind_re[0]]
  ind_fw = where(name eq 'FW',nind_fw)
  info[i].fwhm = value[ind_fw[0]]

  SPAWN,['grep','Clipped',logfile1],out,errout,/noshell
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
    SPAWN,['grep','Clipped',logfile2],out,errout,/noshell
    ;              Clipped mean and median =  187.442  187.215

    ; Delete temporary files
    FILE_DELETE,base+'.find.temp',/allow
    FILE_DELETE,tempscript,/allow
  endif

  arr = strsplit(out[0],' ',/extract)
  if n_elements(arr) ge 7 then begin
    info[i].mnsky = float(arr[5])
    info[i].medsky = float(arr[6])
  endif
  ;; Get sky from ALS file
  if n_elements(arr) lt 7 or info[i].medsky le 0.0 then begin
    alsfile = base+'.als'
    LOADALS,alsfile,als
    info[i].mnsky = mean(als.sky)
    info[i].medsky = median(als.sky)
  endif

  ;; Get exptime and filter
  info[i].exptime = PHOTRED_GETEXPTIME(fitsfile)
  info[i].filter = PHOTRED_GETFILTER(fitsfile)
endfor  ;; file loop

;; Only ONE file, return 1s
if nfiles eq 1 then begin
  actweight = 1.0
  scales = 1.0
  medsky = info[0].medsky
  return
endif

    
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

;; Calculate the magnitude and flux at the 10sigma magnitude
for i=0,nfiles-1 do begin
  gd = where(mag[i,*] lt 50,ngd)
  if ngd gt 0 then begin
    mag1 = mag[i,gd]
    snr1 = 1.087/err[i,gd]
    gdsnr = where(abs(snr1-10.0) lt 1,ngdsnr)
    if ngdsnr lt 5 then gdsnr = where(abs(snr1-10.0) lt 2,ngdsnr)
    if ngdsnr lt 5 then begin
      si = sort(abs(snr1-10.0))
      gdsnr = si[0:99<(ngd-1)]
      ngdsnr = n_elements(gdsnr)
    endif
    mag10 = median([mag1[gdsnr]])
    info[i].mag10 = mag10
    info[i].flux10 = 10.0^( (mag10-25.0)/(-2.5) )          ; total counts
    info[i].fluxrate10 = info[i].flux10 / info[i].exptime  ; flux rate = counts / sec
  endif  
endfor


;; Using TILES
;;--------------
;; We are using TILES and have multiple chips/amps
;;   F1-00507800_39+T2.als, '+T' and two dots
if n_elements(imager) gt 0 then namps=imager.namps else namps=1
if total(stregex(files,'\'+tilesep+'T',/boolean)) eq nfiles and $
   total(long(byte(files[0])) eq btilesep) ge 2 and namps gt 1 then begin
  usetiles = 1

  ;; Number of unique exposures
  expname = strarr(nfiles)
  chip = strarr(nfiles)
  for i=0,nfiles-1 do begin
    base1 = file_basename(files[i],'.als')           ; F1-00507800_39+T2
    field1 = (strsplit(base1,'-',/extract))[0]       ; F1
    expchptile = (strsplit(base1,'-',/extract))[1]   ; 00507800_39+T2
    expchp = (strsplit(expchptile,tilesep,/extract))[0]  ; 00507800_39
    expname[i] = (strsplit(expchp,imager.separator,/extract))[0]
    chip[i] = (strsplit(expchp,imager.separator,/extract))[1]
  endfor
  ;; Unique exposures
  uiexp = uniq(expname,sort(expname))
  uexpname = expname[uiexp]
  nexp = n_elements(uexpname)

  ;; Combine all the catalogs for a given exposure
  ;; Calculate the weights
  expstr = replicate({mag:fltarr(nstars),err:fltarr(nstars),nfiles:0L,index:lonarr(nfiles),$
                      exptime:0.0,filter:'',fwhm:0.0,rdnoise:0.0,medsky:0.0,mag10:99.99,flux10:0.0,fluxrate10:0.0},nexp)
  expmag = fltarr(nexp,nstars)
  experr = fltarr(nexp,nstars)
  for i=0,nexp-1 do begin
    ind = where(expname eq uexpname[i],nind)
    expstr[i].nfiles = nind
    expstr[i].index[0:nind-1] = ind
    expstr[i].filter = info[ind[0]].filter
    expstr[i].exptime = info[ind[0]].exptime
    expstr[i].fwhm = median([info[ind].fwhm])
    expstr[i].rdnoise = median([info[ind].rdnoise])
    expstr[i].medsky = median([info[ind].medsky])
    expstr[i].mag10 = median([info[ind].mag10])
    expstr[i].flux10 = median([info[ind].flux10])
    expstr[i].fluxrate10 = median([info[ind].fluxrate10])
    ;; Combine the photometry
    mag1 = mag[ind,*]
    err1 = err[ind,*]
    ;; Multiple chips
    ;;   they shouldn't overlap, so just use the mean/median
    ;;   and that should pick up the detections
    if nind gt 1 then begin
      bd = where(mag1 gt 50,nbd)
      if nbd gt 0 then mag1[bd]=!values.f_nan
      if nbd gt 0 then err1[bd]=!values.f_nan
      expmag[i,*] = median(mag1,dim=1)
      experr[i,*] = median(err1,dim=1)
    endif else begin
      expmag[i,*] = mag1
      experr[i,*] = err1
    endelse
  endfor
  ;; Replace NANs with 99.9999
  bdmag = where(finite(expmag) eq 0,nbdmag)
  expmag[bdmag] = 99.99
  experr[bdmag] = 9.99
  expstr.mag = transpose(expmag)
  expstr.err = transpose(experr)
  ;; Perform the weights and scales calculations
  ALLFRAME_GETWEIGHTS_RAW,expstr,outexpstr
  ;; Give each chip the weight of its respective exposure
  for i=0,nexp-1 do begin
    info[expstr[i].index[0:expstr[i].nfiles-1]].weight = outexpstr[i].weight
    info[expstr[i].index[0:expstr[i].nfiles-1]].scale = outexpstr[i].scale
  endfor

;; REGULAR Method
;;---------------
Endif else begin
  str = replicate({mag:fltarr(nstars),err:fltarr(nstars),$
                   exptime:0.0,filter:'',fwhm:0.0,rdnoise:0.0,medsky:0.0,flux10:0.0,fluxrate10:0.0},nfiles)
  struct_assign,info,str
  str.mag = transpose(mag)
  str.err = transpose(err)
  ;; Perform the weights and scales calculations
  ALLFRAME_GETWEIGHTS_RAW,str,outstr
  info.weight = outstr.weight
  info.scale = outstr.scale
Endelse


;; Print out the information
if not keyword_set(silent) then begin
  printlog,logfile,'        FILE         FILTER EXPTIME FWHM RDNOISE MEDSKY WEIGHT SCALE'
  for i=0,nfiles-1 do begin
    printlog,logfile,info[i].name,info[i].filter,info[i].exptime,info[i].fwhm,info[i].rdnoise,$
             info[i].medsky,info[i].weight,info[i].scale,format='(A-23,A4,F6.1,F6.2,F6.2,F8.1,F7.3,F7.3)'
  endfor
endif

;; Final output
actweight = info.weight
scales = info.scale
medsky = info.medsky

if keyword_set(stp) then stop

end
