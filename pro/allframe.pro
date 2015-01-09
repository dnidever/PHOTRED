;+
;
; ALLFRAME
;
; This runs ALLFRAME on images
;
; You need to have run daophot, allstar, daomatch and daomaster
; already.  There needs to be a fits, opt, als.opt, ap and als
; file for each file in the MCH file.  There also needs to be an
; associated RAW file for the MCH file.
;
; INPUTS:
;  file           The MCH filename
;  =finditer      The maximum number of iterations to use when finding
;                   sources with SExtractor/ALLSTAR.  The default is 2,
;                   and the maximum allowed it 10.
;  =detectprog    The program to use to detect sources.  Either
;                   'sextractor' or 'daophot'.  'sextractor' is the
;                   default.  'daophot' is better for VERY crowded
;                   regions.
;  =nocmbimscale  Don't scale the images when combining them.  Not
;                   recommended, but the old way of doing it.  Bright
;                   stars can be missed this way.
;  /combtrim      Trim the combined images to the overlapping region.
;                   This used to be the default, but now the default
;                   is to keep the entire original region.
;  =scriptsdir    The directory that contains all of the necessary scripts.
;  =irafdir       The IRAF home directory.
;  =logfile       A logfile to print to output to.
;  /stp           Stop at the end of the program
;
; OUTPUTS:
;  The final allframe output file name is filebase+'.mag'
;  =error  The error message, if there was one, else undefined
;
; USAGE:
;  IDL>allframe,'ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/',finditer=2
;
;
; By D.Nidever   February 2008 
; Automation of steps and scripts by J.Ostheimer and Rachael Beaton
;-

pro allframe_getweights,mchfile,actweight,scales,medsky,raw2=raw2

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

end

    
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

print,'Files: ',files
print,'Weights: ',actweight
print,'Scales: ',scales
print,'Sky: ',medsky

if keyword_set(stp) then stop

end

;####################################################################################

pro allframe,file,stp=stp,scriptsdir=scriptsdir,detectprog=detectprog,$
             error=error,logfile=logfile,finditer=finditer0,$
             irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,trimcomb=trimcomb

COMMON photred,setup

undefine,error

; Not enough inputs
nfile = n_elements(file)
if (nfile eq 0) then begin
  print,'Syntax - allframe,file,stp=stp,scriptsdir=scriptsdir,finditer=finditer,satlevel=satlevel,'
  print,'                  detectprog=detectprog,nocmbimscale=nocmbimscale,error=error,logfile=logfile,'
  print,'                  irafdir=irafdir,trimcomb=trimcomb'
  return
endif

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'ALLFRAME ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return
endif

; How many FIND iterations
if n_elements(finditer0) eq 0 then finditer=2 else finditer=finditer0
finditer = finditer < 10  ; maximum 10.

; Saturation level
if n_elements(satlevel) eq 0 then satlevel=6e4

; Scaling of the images to be combined
if n_elements(nocmbimscale) eq 0 then nocmbimscale=0

; Getting scripts directory and iraf directory
nsetup = n_elements(setup)
if nsetup gt 0 then begin
  scriptsdir = READPAR(setup,'SCRIPTSDIR')
  irafdir = READPAR(setup,'IRAFDIR')
endif


; No irafdir
if n_elements(scriptsdir) eq 0 then begin
  printlog,logf,'IRAFDIR NOT INPUT'
  error = 'IRAFDIR NOT INPUT'
  return
endif

; Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL
;---------------------------------------------------------------
CHECK_IRAF,iraftest,irafdir=irafdir
if iraftest eq 0 then begin
  print,'IRAF TEST FAILED.  EXITING'
  return
endif


; No scriptsdir
if n_elements(scriptsdir) eq 0 then begin
  printlog,logf,'SCRIPTSDIR NOT INPUT'
  error = 'SCRIPTSDIR NOT INPUT'
  return
endif
; Check if the scripts exist in the current directory
scripts = ['getpsf.sh','photo.opt','apcor.opt','lstfilter','goodpsf.pro','allframe.opt',$
           'default.sex','default.param','default.nnw','default.conv','makemag']
nscripts = n_elements(scripts)
; Loop through the scripts
for i=0,nscripts-1 do begin
  info = FILE_INFO(scriptsdir+'/'+scripts[i])
  curinfo = FILE_INFO(scripts[i])

  ; No file
  if info.exists eq 0 or info.size eq 0 then begin
    printlog,logf,scriptsdir+'/'+scripts[i],' NOT FOUND or EMPTY'
    error = scriptsdir+'/'+scripts[i]+' NOT FOUND or EMPTY'
    return
  endif

  ; Check if the two files are the same size, if not copy it
  if info.size ne curinfo.size then begin
    FILE_COPY,info.name,curinfo.name,/overwrite
  endif
end ; scripts loop


; Check that the ALLFRAME program exists
SPAWN,'which allframe',out,errout
allframefile = FILE_SEARCH(out,count=nallframefile)
;allframefile = FILE_SEARCH('/net/halo/bin/allframe.2008',count=nallframefile)
;allframefile = FILE_SEARCH('/net/halo/bin/allframe.2004.fixed',count=nallframefile)
if (nallframefile eq 0) then begin
  ;printlog,logf,'/net/halo/bin/allframe.2004.fixed NOT FOUND'
  ;printlog,logf,'/net/halo/bin/allframe.2008 NOT FOUND'
  printlog,logf,allframefile+'NOT FOUND'
  return
endif


printlog,logf,''
printlog,logf,'====================================='
printlog,logf,'RUNNING ALLFRAME on ',file
printlog,logf,'====================================='
printlog,logf,''


; FILENAME
mchfile = file_basename(file)
mchdir = file_dirname(file)
mchbase = file_basename(file,'.mch')


; CD to the directory
cd,current=curdir
cd,mchdir


; Check that the mch, als, and opt files exist
mchtest = file_test(mchfile)
if mchtest eq 0 then begin
  printlog,logf,mchfile,' NOT FOUND'
  return
endif



;###########################################
; CHECK NECESSARY FILES

; Load the MCH file
LOADMCH,mchfile,files,trans

; Check that the fits, als, and opt files exist
nfiles = n_elements(files)
for i=0,nfiles-1 do begin
  ;dir = file_dirname(mchfile)
  base = file_basename(files[i],'.als')
  
  ; Checking FITS file
  fitstest = file_test(base+'.fits')
  if fitstest eq 0 then begin
    printlog,logf,base+'.fits NOT FOUND'
    return
  endif

  ; Checking OPT file
  opttest = file_test(base+'.opt')
  if opttest eq 0 then begin
    printlog,logf,base+'.opt NOT FOUND'
    return
  endif

  ; Checking ALS.OPT file
  alsopttest = file_test(base+'.als.opt')
  if alsopttest eq 0 then begin
    printlog,logf,base+'.als.opt NOT FOUND'
    return
  endif

  ; Checking AP file
  aptest = file_test(base+'.ap')
  if aptest eq 0 then begin
    printlog,logf,base+'.ap NOT FOUND'
    return
  endif

  ; Checking ALS file
  alstest = file_test(base+'.als')
  if alstest eq 0 then begin
    printlog,logf,base+'.als NOT FOUND'
    return
  endif

  ; Checking LOG file
  logtest = file_test(base+'.log')
  if logtest eq 0 then begin
    printlog,logf,base+'.log NOT FOUND'
    return
  endif

  ; Checking RAW file
  rawtest = file_test(base+'.raw')
  if logtest eq 0 then begin
    printlog,logf,base+'.raw NOT FOUND'
    return
  endif

  ; REMOVE ALF if it exists
  if FILE_TEST(base+'.alf') then FILE_DELETE,base+'.alf',/allow,/quiet

end

; REMOVE the .mag file if it exists
if FILE_TEST(mchbase+'.mag') then FILE_DELETE,mchbase+'.mag',/allow,/quiet




;###########################################
; STEP 1: IMALIGN PREP

printlog,logf,'------------------------'
printlog,logf,'STEP 1: GETTING WEIGHTS'
printlog,logf,'------------------------'

;-----------------------------------
; Computs Weights
ALLFRAME_GETWEIGHTS,mchfile,weights,scales,sky,raw2=raw2
invscales = 1.0/scales
bdscale = where(scales lt 1e-5 or invscales gt 900,nbdscale)
if nbdscale gt 0 then begin
  scales[bdscale] = 1.0
  invscales[bdscale] = 1.0
  weights[bdscale] = 0.0
endif
weightfile = mchbase+'.weights'
WRITECOL,weightfile,weights,fmt='(F10.6)'
scalefile = mchbase+'.scale'
WRITECOL,scalefile,invscales,fmt='(F10.5)'  ; want to scale it UP
zerofile = mchbase+'.zero'
WRITECOL,zerofile,-sky,fmt='(F10.2)'  ; want to remove the background, set to 1st frame


;---------------------------------------
; Get X/Y translations using DAOMASTER
;  NO ROTATION ONLY SHIFTS
;  Need these shifts for IMSHIFT
print,'Measuring X/Y shifts'
shiftmch = mchbase+'_shift'
FILE_COPY,mchbase+'.mch',shiftmch+'.mch',/overwrite,/allow
; Make the DAOMASTER script
undefine,cmdlines
PUSH,cmdlines,'#!/bin/csh'
PUSH,cmdlines,'set input=${1}'
PUSH,cmdlines,'daomaster << DONE'
PUSH,cmdlines,'${input}.mch'
PUSH,cmdlines,'1,1,1'
PUSH,cmdlines,'99.'
PUSH,cmdlines,'2'
PUSH,cmdlines,'10'
PUSH,cmdlines,'5'
PUSH,cmdlines,'4'
PUSH,cmdlines,'3'
PUSH,cmdlines,'2'
PUSH,cmdlines,'1'
PUSH,cmdlines,'0'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'y'
PUSH,cmdlines,''
PUSH,cmdlines,''
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'DONE'
tempscript = MKTEMP('daomaster')   ; absolute filename
WRITELINE,tempscript,cmdlines
FILE_CHMOD,tempscript,'755'o
; Run DAOMASTER                                                                                                                           
cmd2 = tempscript+' '+shiftmch
SPAWN,cmd2,out2,errout2
; Remove temporary DAOMASTER script
FILE_DELETE,tempscript,/allow_non
LOADMCH,shiftmch+'.mch',files2,trans2

xshift = reform(trans2[*,0])
yshift = reform(trans2[*,1])
xyshifts = [[xshift],[yshift]]
printlog,logf,'Image shifts'
for i=0,nfiles-1 do printlog,logf,files[i],xshift[i],yshift[i]


;-----------------------------------
; Create imalign prep files
; This is done by preimalign_k.sh
; Need an input list of fits files
; Need an output list of fits files
; Shift file
base = file_basename(files,'.als')
fitsfiles = base+'.fits'
outfiles = base+'.shft.fits'
infile = mchbase+'.inlist'
outfile = mchbase+'.outlist'
;WRITELINE,infile,fitsfiles   ; this is done below now with the temp files
WRITELINE,outfile,outfiles
; Remove outfiles
FILE_DELETE,outfiles,/allow
; shift list
shiftfile = mchbase+'.shift'
WRITECOL,shiftfile,xshift,yshift,fmt='(2F15.4)'


; Make temporary files for bad pixel fixing and combining
;  FIXIMAGES doesn't work properly on the shifted images
;  because the interpolation can bring the bad pixel values down
;  below the saturation threshold, and we don't want to touch
;  the original images.
tempfits = base+'.temp.fits'
FILE_COPY,fitsfiles,tempfits,/overwrite,/allow
WRITELINE,infile,tempfits



;###########################################
; STEP 2: FIX BAD PIXELS
printlog,logf,'---------------------------'
printlog,logf,'STEP 2: FIXING BAD PIXELS'
printlog,logf,'---------------------------'
FIXIMAGES,'@'+infile,satlevel=satlevel  ;6e4
; This also makes the FILE.mask.fits files for each image

; Find the maximum saturation level
satlevelarr = fltarr(nfiles)
for i=0,nfiles-1 do begin
  ;head = headfits(base[i]+'.fits')
  FITS_READ,base[i]+'.fits',im,head,/no_abort
  saturate = sxpar(head,'SATURATE',count=nsaturate)
  if nsaturate eq 0 then saturate=max(im)-1000.
  satlevelarr[i] = saturate
endfor
maxsatlevel = max(satlevelarr)


;###########################################
; STEP 3: IMALIGN
;  This figures out the X/Y-shifts between the images
;  and creates the shifted images (".shft.fits")
printlog,logf,'-----------------'
printlog,logf,'STEP 3: IMALIGN'
printlog,logf,'-----------------'
reffile = mchbase+'.fits'

; IMALIGN basically is a script that runs:
;  IMCENTROID - to compute the shifts
;  IMSHIFT - shifts the images
;  IMCOPY - trims the images


; First, shift the images
printlog,logf,'Shifting the images'
IRAF_IMSHIFT,'@'+infile,'@'+outfile,shifts_file=shiftfile,interp_type='linear',$
             boundary_type='constant',constant=0,irafdir=irafdir,error=imshifterror
if n_elements(imshifterror) ne 0 then begin
  printlog,logf,'ERROR in IRAF_IMSHIFT'
  printlog,logf,imshifterror
  error = imshifterror
  return
endif

; Trim the images
if keyword_set(trimcomb) then begin

  ; Calculate the trim section
  hd = headfits(reffile)
  xsize = lonarr(nfiles)+sxpar(hd,'NAXIS1')
  ysize = lonarr(nfiles)+sxpar(hd,'NAXIS2')
  IA_TRIM,xshift,yshift,xsize,ysize,trimsection
  xoff = trimsection[0]-1
  yoff = trimsection[2]-1

  ; Trim the shifted images
  printlog,logf,'Trimming the shifted images'
  xstart = trimsection[0]-1
  xstop = trimsection[1]-1
  ystart = trimsection[2]-1
  ystop = trimsection[3]-1

  for i=0,nfiles-1 do begin
    FITS_READ,outfiles[i],im,head
    newim = im[xstart:xstop,ystart:ystop]
    MWRFITS,newim,outfiles[i],head,/create,/silent
  end
  ; could also use IRAF_IMCOPY here instead


; Don't trim the images
endif else begin
  hd = headfits(reffile)
  xsize = sxpar(hd,'NAXIS1')
  ysize = sxpar(hd,'NAXIS2')
  trimsection = [1,xsize,1,ysize]
  xoff = 0
  yoff = 0
endelse

; Delete the temporary FITS files
FILE_DELETE,tempfits,/allow



;###########################################
; STEP 4: MAKE BAD PIXEL/WEIGHT MAP
; in how many images does the pixel need to be bad??
; 1. shift the masks (created by FIXIMAGES.PRO)
;     using the shifts from IRAF_IMALIGN
; 2. trim the masks
; 3. combine the masks
printlog,logf,'-------------------------------'
printlog,logf,'STEP 4: Making BAD PIXEL MASK'
printlog,logf,'-------------------------------'

; Make lists
maskfiles = FILE_DIRNAME(tempfits)+'/'+FILE_BASENAME(tempfits,'.fits')+'.mask.fits'
outmaskfiles = FILE_DIRNAME(tempfits)+'/'+FILE_BASENAME(tempfits,'.fits')+'.mask.shft.fits'
maskinfile = mchbase+'.maskinlist'
maskoutfile = mchbase+'.maskoutlist'
maskshiftsfile = mchbase+'.maskshifts'
WRITELINE,maskinfile,maskfiles
WRITELINE,maskoutfile,outmaskfiles
FILE_DELETE,outmaskfiles,/allow
strxyshifts = strarr(nfiles)
for i=0,nfiles-1 do strxyshifts[i] = strjoin(reform(xyshifts[i,*]),'  ')
WRITELINE,maskshiftsfile,strxyshifts

; Run IMSHIFT
;  set boundary to 0=BAD
printlog,logf,'Shifting masks'
undefine,iraflines
push,iraflines,'cd '+curdir
push,iraflines,'images'
push,iraflines,'imgeom'
push,iraflines,'imshift("@'+maskinfile+'","@'+maskoutfile+'",shifts_file="'+maskshiftsfile+'",'+$
               'interp_type="linear",boundary_typ="constant",constant=0)'
push,iraflines,'logout'
imshiftscript = curdir+'/'+mchbase+'.imshift'
WRITELINE,imshiftscript,iraflines
IRAF_RUN,imshiftscript,irafdir,out=out,/silent,error=iraferror
if n_elements(iraferror) ne 0 then begin
  printlog,logf,'ERROR in running IMSHIFT with IRAF_RUN'
  printlog,logf,iraferror
  error = iraferror
  return
endif

; Trim
if keyword_set(trimcomb) then begin
  printlog,logf,'Trimming masks'
  xstart = trimsection[0]-1     ; should be same as xoff
  xstop = trimsection[1]-1
  ystart = trimsection[2]-1     ; should be same as yoff
  ystop = trimsection[3]-1

  for i=0,nfiles-1 do begin
    FITS_READ,outmaskfiles[i],im,head
    sz = size(im)
    newim = im[xstart:xstop,ystart:ystop]
    ; Add LTV1/LTV2 to the header
    ;  these are IRAF keywords to convert from logical to physical coords
    ltv1 = sxpar(head,'LTV1')  ; 0 if not found
    ltv2 = sxpar(head,'LTV2')  ; 0 if not found
    sxaddpar,head,'LTV1',ltv1-xstart
    sxaddpar,head,'LTV2',ltv2-ystart
    MWRFITS,newim,outmaskfiles[i],head,/create,/silent
  endfor
endif

; Combining masks
printlog,logf,'Combining masks'
undefine,bpm
for i=0,nfiles-1 do begin
  FITS_READ,outmaskfiles[i],im,head
  if i eq 0 then begin
    bpm = im
    whead = head
  endif else begin
    bpm = bpm+im
  endelse
end
;bpm = bpm/float(nfiles)

;; masks have 0-bad, 1-good.
;; anything with less than 1.0 is considered bad
;; weight map, -1 is bad, +1 is good
;weightmap = -2.0*float(bpm lt nfiles) + 1.
;combweightfile = mchbase+'_comb.mask.fits'
;FITS_WRITE,combweightfile,weightmap,whead
;
; THIS IS NOW DONE BELOW AFTER THE IMAGE IS COMBINED
; DEPENDING ON IF THE IMAGES ARE SCALED OR NOT!!!


;###########################################
; STEP 5: COMBINE IMAGES
printlog,logf,'-------------------'
printlog,logf,'STEP 5: IMCOMBINE'
printlog,logf,'-------------------'


; SCALE the images for combining
;-------------------------------
if not keyword_set(nocmbimscale) then begin

  ; Put BPM mask names in file headers
  ;  these will be used by IMCOMBINE
  for i=0,nfiles-1 do begin
    head = headfits(outfiles[i])
    sxaddpar,head,'BPM',outmaskfiles[i]
    modfits,outfiles[i],0,head
  end

  ; Combine the frames WITH scaling/offset/masking, for the bright stars
  ;printlog,logf,'Creating SCALED image'
  combfile = mchbase+'_comb.fits'
  FILE_DELETE,combfile,/allow
  FILE_DELETE,mchbase+'_comb.bpm.pl',/allow
  IRAF_IMCOMBINE,'@'+outfile,combfile,combine='average',reject='avsigclip',$
                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',$
                 irafdir=irafdir,error=imcombineerror2,scale='@'+scalefile,zero='@'+zerofile,$
                 masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm'

  if n_elements(imcombineerror2) ne 0 then begin
    printlog,logf,'ERROR in IRAF_IMCOMBINE'
    printlog,logf,imcombineerror2
    error = imcombineerror2
    return
  endif

  ; Convert BPM mask from PL to FITS
  FILE_DELETE,mchbase+'_comb.bpm.fits',/allow
  undefine,lines
  cd,current=curdir
  push,lines,'cd '+curdir
  push,lines,'imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits'
  push,lines,'logout'
  tempfile = mktemp('tiraf')
  WRITELINE,tempfile,lines
  IRAF_RUN,tempfile,irafdir,silent=silent,out=out,error=error

  ; Delete temporary scripts and PL file
  FILE_DELETE,[tempfile,mchbase+'_comb.bpm.pl'],/allow


  ; Fix the rdnoise and background/sky level and saturate
  ;  the bad pixels for DAOPHOT
  ;------------------------------------------------------

  ; 10/02/12
  ; THE IMAGES ARE (1) ZERO-SUBTRACTED, (2) SCALED, AND (3) WEIGHT AVERAGED
  ; The algorithm is:
  ; 1.) add zero-level correction.  im = im+zero
  ; 2.) scale the images.  im = im*scale
  ; 3.) take weighted average.  combim=total(weight*im)
  ;      there is also clipping that takes place during the averaging
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2))
  ; A gain that changes from frame to frame could be problematic,
  ; but this shouldn't happen since it's the same chip from the same night.

  ; IMCOMBINE wants rdnoise in electrons and gain in electrons/DN.
  ; DAOPHOT expects rdnoise in DN.  That's why mkopt converts
  ;  it with the gain.  So we are fine.  The header should have
  ;  rdnoise in ELECTRONS.

  ; page 63-65 of daophot2.pdf shows how you need to modify rdnoise/gain
  ; when averaging/summing frames. in observing/mosaic/.

  ; Load the IMCOMBINE output combined file and BPM
  FITS_READ,combfile,combim,combhead
  FITS_READ,mchbase+'_comb.bpm.fits',badmask,maskhead  ; 0-good, 1-bad


  ; Fix the gain
  ; For N averaged frames gain(N)=N*gain(1)
  ; Leave the gain as is!  We are scaling everything to the reference
  ; and using its gain.  It's nearly impossible to figure out the real
  ; gain since we are scaling the images and then taking a weighted
  ; average with outlier rejection.  Find a gain that properly
  ; describes/follows Poisson errors for the final combined image is
  ; difficult/impossible.  But that's okay.  This is just for source
  ; detection and DAOPHOT FIND just cares about the noise in the
  ; background.  We just need to ensure that the sky and rdnoise
  ; are correct.

  ; Fix the rdnoise
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2))
  rdnoisearr = fltarr(nfiles)
  for i=0,nfiles-1 do rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits')
  ;  the "scales" array here is actually 1/scales used by IMCOMBINE.
  rdnoise = sqrt(total((weights*rdnoisearr/scales)^2))
  dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey) ; get keyword
  sxaddpar,combhead,rdnoisekey,rdnoise

  ; Fix the sky
  ; DAOPHOT FIND computes the random error per pixel in ADU as
  ; noise = sqrt( sky level/gain + rdnoise^2)
  ; So it assumes that the noise in the background is sqrt(sky/gain)
  ; in ADU.  We need to set the sky level so this is correct.
  ; The final noise should be 
  ; final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
  ; So the final sky level should be
  ; final sky = total((weights*scale*sqrt(sky/gain))^2)*gain
  gain = PHOTRED_GETGAIN(combfile,keyword=gainkey)
  comb_sky = total((weights*sqrt((sky>0)/gain)/scales)^2)*gain
  ; the "scales" array here is actually 1/scale
  combim += float(comb_sky)  ; keep it float


  ; set the maximum to a "reasonable" level
  ; Rescale the image and increase the gain
  if max(combim) gt 50000 then begin
    rescale = 50000./max(combim)
    combim = combim*rescale
    sxaddpar,combhead,gainkey,gain/rescale
    ; rdnoise does NOT get modified since it's in electrons
    ; we just need to modify the gain which takes you from ADUs to electrons
  endif

  maskdatalevel = max(combim) + 10000       ; set "bad" data level above the highest "good" value
  combim2 = combim*(1-badmask) + maskdatalevel*badmask    ; set bad pixels to maskdatalevel
  MWRFITS,combim2,combfile,combhead,/create  ; fits_write can create an empty PDU

  ; Create the weight map for Sextractor using the BPM output by IMCOMBINE
  ;  bad only if bad in ALL images
  weightmap = -2.0*float(badmask eq 1) + 1.0
  combweightfile = mchbase+'_comb.mask.fits'
  MWRFITS,weightmap,combweightfile,whead,/create

; NO SCALING of the images for combining
;---------------------------------------
Endif else begin

  combfile = mchbase+'_comb.fits'
  FILE_DELETE,combfile,/allow
  IRAF_IMCOMBINE,'@'+outfile,combfile,combine='average',reject='avsigclip',$
                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',$
                 irafdir=irafdir,error=imcombineerror2

  if n_elements(imcombineerror2) ne 0 then begin
    printlog,logf,'ERROR in IRAF_IMCOMBINE'
    printlog,logf,imcombineerror2
    error = imcombineerror2
    return
  endif

  ; Fix the rdnoise and background/sky level and saturate
  ;  the bad pixels for DAOPHOT
  ;------------------------------------------------------

  ; See the explanations for all these steps above!!

  ; Load the IMCOMBINE output combined file and BPM
  FITS_READ,combfile,combim,combhead
  FITS_READ,mchbase+'_comb.bpm.fits',badmask,maskhead  ; 0-good, 1-bad

  ; Fix the rdnoise
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2))
  rdnoisearr = fltarr(nfiles)
  for i=0,nfiles-1 do rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits')
  rdnoise = sqrt(total((weights*rdnoisearr)^2))
  dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey) ; get keyword
  sxaddpar,combhead,rdnoisekey,rdnoise

  ; Fix the sky
  ; So the final sky level should be
  ; final sky = total((weights*scale*sqrt(sky/gain))^2)*gain
  gain = PHOTRED_GETGAIN(combfile,keyword=gainkey)
  comb_sky = total((weights*sqrt(sky/gain))^2)*gain
  ; the "scales" array here is actually 1/scale
  combim += comb_sky


  ; set the maximum to a "reasonable" level
  ; Rescale the image and increase the gain
  if max(combim) gt 50000 then begin
    rescale = 50000./max(combim)
    combim = combim*rescale
    sxaddpar,combhead,gainkey,gain/rescale
    ; rdnoise does NOT get modified since it's in electrons
    ; we just need to modify the gain which takes you from ADUs to electrons
  endif


  ; Making Sextractor "weight" map file
  ;------------------------------------
  ; masks have 0-bad, 1-good.
  ; anything with less than 1.0 is considered bad
  ; weight map, -1 is bad, +1 is good
  ; "bpm" is the SUM of the bad pixel masks
  ; consider a pixel bad that is bad in ANY image
  weightmap = -2.0*float(bpm lt nfiles) + 1.
  combweightfile = mchbase+'_comb.mask.fits'
  FITS_WRITE,combweightfile,weightmap,whead

  ;---------------------------------------------
  ; SATURATE BAD pixels in the COMBINED IMAGE
  ; DAOPHOT needs to have the bad pixels "saturated",
  ; SExtractor will know which pixels are bad from the "weight" map.
  ;
  ; We could skip the fiximage.pro step but we still need the
  ; individual bpm masks and setting the bad pixels to the background
  ; probably helps in the IMALIGN/IMCOMBINE steps.
  printlog,logf,''
  printlog,logf,'"Saturating" bad pixels in the COMBINED image'
  printlog,logf,''


  badmask = float(weightmap lt 0.5)
  maskdatalevel = max(combim) + 10000       ; set "bad" data level above the highest "good" value
  combim2 = combim*(1.0-badmask) + maskdatalevel*badmask    ; set bad pixels to 100,000
  FITS_WRITE,combfile,combim2,combhead

Endelse ; no scaling of images for combining

; Delete the shifted images
READLINE,outfile,shiftedfiles
FILE_DELETE,shiftedfiles,/allow,/quiet

; Delete mask files
FILE_DELETE,[maskfiles,outmaskfiles],/allow
FILE_DELETE,[maskinfile,maskoutfile,maskshiftsfile,imshiftscript],/allow



;###########################################
; STEP 6: Get PSF for Combined Image
printlog,logf,'----------------------------------------'
printlog,logf,'STEP 6: Getting PSF for Combined Image'
printlog,logf,'----------------------------------------'
; Make .opt files, set saturation just below the mask data level
MKOPT,combfile,satlevel=maskdatalevel-1000
; Using CMN.LST of reference frame if it exists
if file_test(mchbase+'.cmn.lst') then begin
  print,'Using reference image COMMON SOURCE file'
  file_copy,mchbase+'.cmn.lst',mchbase+'_comb.cmn.lst',/over,/allow
endif
; Get the PSF of the combined image
SPAWN,'./getpsf.sh '+file_basename(combfile,'.fits')




;###########################################
; STEP 7: Run allframe prep
;  This iteratively runs SExtractor on the combined image
;  This can take a while.
printlog,logf,'--------------------------------'
printlog,logf,'STEP 7: Running allframe prep'
printlog,logf,'--------------------------------'

ALLFPREP,combfile,als,xoff,yoff,logfile=logfile,error=error,$
         detectprog=detectprog,scriptsdir=scriptsdir,maxiter=finditer,$
         maskfile=combweightfile
if n_elements(error) gt 0 then goto,BOMB



;###########################################
; STEP 8: Running ALLFRAME
printlog,logf,'----------------------------'
printlog,logf,'STEP 8: Running ALLFRAME'
printlog,logf,'----------------------------'

; What we need
; allf.mag     List of coordinates made by allfprep
; allf.mch     List of transformations
; allframe.opt
; obj????.psf
; obj????.als
; obj????.fits

; Delete any temporary ALLFRAME files from possible
; previous runs of allframe.  Otherwise ALLFRAME
; will start from where it left off.
; For each ALLFRAME run there is:
;  mchbasename+'.bck'
;  mchbasename+'.nmg'
;  mchbasename+'.tfr'
; For each file in the MCH file there are:
;  filebasename+'.alf'
;  filebasename+'j.fits'
;  filebasename+'k.fits'
if file_test(mchbase+'.tfr') eq 1 then $
  FILE_COPY,mchbase+'.tfr',mchbase+'.tfr.orig',/over,/allow  ; copy original
FILE_DELETE,base+'.bck',/allow
FILE_DELETE,base+'.nmg',/allow
FILE_DELETE,base+'.tfr',/allow   ; This gets overwritten!!!
FILE_DELETE,base+'j.fits',/allow
FILE_DELETE,base+'k.fits',/allow
FILE_DELETE,base+'.alf',/allow

; Make input file
undefine,cmd
push,cmd,'    '
push,cmd,mchfile               ; mch file
push,cmd,mchbase+'_comb_allf.als'  ; coord file
push,cmd,'    '
;cmdfile = maketemp('temp','.inp')
cmdfile = MKTEMP('temp')
WRITELINE,cmdfile,cmd
;SPAWN,'/net/halo/bin/allframe.2004.fixed < '+cmdfile
;SPAWN,'/net/halo/bin/allframe.2008 < '+cmdfile
SPAWN,'allframe < '+cmdfile

FILE_DELETE,cmdfile,/allow
FILE_DELETE,file_basename(files,'.als')+'j.fits',/allow   ; delete subtracted images



;###########################################
; STEP 9: Combine photometry with MAKEMAG
; This combines the photometry in the N alf files
; and averages chi and sharp
printlog,logf,'--------------------------'
printlog,logf,'STEP 9: Running MAKEMAG'
printlog,logf,'--------------------------'

;FILE_COPY,scriptsdir+'makemag','.',/overwrite
FILE_DELETE,mchbase+'.makemag',/allow

;; Make input file
;magfile = mchbase+'.makemag'
;undefine,cmd
;push,cmd,mchbase+'.tfr'           ; final tfr file
;push,cmd,strtrim(nfiles,2)+',0'   ; nfiles, offset
;push,cmd,magfile                  ; final output file
;push,cmd,'2'                      ; do not renumber
;;cmdfile = maketemp('temp','.inp')
;cmdfile = MKTEMP('temp')
;WRITELINE,cmdfile,cmd
;SPAWN,'./makemag < '+cmdfile
;FILE_DELETE,cmdfile,/allow        ; delete temporary input file

magfile = mchbase+'.makemag'
MAKEMAG,mchbase+'.tfr',magfile


; Prepend the ALF header to the makemag file
line1='' & line2='' & line3=''
openr,unit,/get_lun,mchbase+'.alf'
readf,unit,line1
readf,unit,line2
readf,unit,line3
close,unit
free_lun,unit
head = [line1,line2,line3]
WRITELINE,magfile,head,/prepend



;######################################################
; STEP 10: Adding SExtractor information
printlog,logf,'----------------------------------------'
printlog,logf,'STEP 10: Adding SExtractor information'
printlog,logf,'----------------------------------------'

; combfile_allf.sex can be matched to the makemag file using IDs
; Load the SExtractor file
sexfile = mchbase+'_comb_allf.sex'
if FILE_TEST(sexfile) eq 1 then begin

  fields = ['ID','X','Y','MAG','ERR','FLAG','PROB']
  sex = IMPORTASCII(sexfile,fieldnames=fields,/noprint)
  nsex = n_elements(sex)

  ; Load the MAKEMAG file
  LOADRAW,mchbase+'.makemag',mag,alfhead
  nmag = n_elements(mag)

  ; Match them with IDs
  MATCH,mag.id,sex.id,ind1,ind2,count=nind

  ; Add stellaricity information to mag file
  add_tag,mag,'flag',0L,mag
  add_tag,mag,'prob',0.0,mag
  mag[ind1].flag = sex[ind2].flag
  mag[ind1].prob = sex[ind2].prob

  if nind lt nmag then printlog,logf,'DID NOT MATCH ALL THE STARS!'


  ; Write the final output file
  ;----------------------------
  finalfile = mchbase+'.mag'
  ; How many observations are there
  tags = tag_names(mag)
  ntags = n_elements(tags)
  magind = where(stregex(tags,'^MAG',/boolean) eq 1,nmagind)

  ; Copy the structure to a string array, then print it out
  outarr = strarr(ntags,nmag)
  fmtarr = '('+['I9','F9.3','F9.3',replicate('F9.4',nmagind*2),'F9.4','F9.4','I5','F7.2']+')'
  outfmt='(A9,2A9,'+strtrim(nmagind*2,2)+'A9,2A9,A5,A7)'
  for i=0,ntags-1 do outarr[i,*] = STRING(mag.(i),format=fmtarr[i])
  openw,unit,/get_lun,finalfile
  printf,unit,format=outfmt,outarr
  close,unit
  free_lun,unit

  ; Prepend the ALF header
  WRITELINE,finalfile,[alfhead,' '],/prepend

; DAOPHOT
endif else begin

  ; No SExtractor information, just copy .makemag to .mag
  finalfile = mchbase+'.mag'
  FILE_COPY,mchbase+'.makemag',finalfile,/allow,/over

endelse

printlog,logf,'FINAL ALLFRAME file = ',finalfile

BOMB:

if keyword_set(stp) then stop

end
