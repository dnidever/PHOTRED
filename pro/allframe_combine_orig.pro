;+
;
; ALLFRAME_COMBINE_ORIG
;
; This is the old version of the code that combines images ALLFRAME.
;
; You need to have run daophot, allstar, daomatch and daomaster
; already.  There needs to be a fits, opt, als.opt, ap and als
; file for each file in the MCH file.  There also needs to be an
; associated RAW file for the MCH file.
;
; INPUTS:
;  file            The MCH filename
;  =nocmbimscale   Don't scale the images when combining them.  Not
;                    recommended, but the old way of doing it.  Bright
;                    stars can be missed this way.
;  /combtrim       Trim the combined images to the overlapping region.
;                    This used to be the default, but now the default
;                    is to keep the entire original region.
;  =scriptsdir     The directory that contains all of the necessary scripts.
;  =irafdir        The IRAF home directory.
;  =logfile        A logfile to print to output to.
;  /fake           Run for artificial star tests.
;  /stp            Stop at the end of the program
;
; OUTPUTS:
;  The combined image and mask:
;    FILEBASE.comb.fits       combined image
;    FILEBASE.comb.bpm.fits   bad pixel mask
;    FILEBASE.comb.mask.fits  SExtractor weight map (very similar to bpm)
;    FILEBASE.comb.mch        The transformations of the individual frames
;                               to the combined frame.
;  =maskdatalevel  The "bad" data level above the highest "good" value 
;                    in the combined image.
;  =error          The error message, if there was one, else undefined
;
; USAGE:
;  IDL>allframe_combine_orig,'ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/',finditer=2
;
;
; By D.Nidever   February 2008 
; Automation of steps and scripts by J.Ostheimer and Rachael Beaton
;-


pro allframe_combine_orig,file,stp=stp,scriptsdir=scriptsdir,error=error,logfile=logfile,$
                          irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,$
                          trimcomb=trimcomb,fake=fake,maskdatalevel=maskdatalevel,$
                          xoff=xoff,yoff=yoff

COMMON photred,setup

undefine,error

; Not enough inputs
nfile = n_elements(file)
if (nfile eq 0) then begin
  print,'Syntax - allframe_combine_orig,file,stp=stp,scriptsdir=scriptsdir,satlevel=satlevel,'
  print,'                  ,nocmbimscale=nocmbimscale,error=error,logfile=logfile,'
  print,'                  irafdir=irafdir,trimcomb=trimcomb,usecmn=usecmn,fake=fake'
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
   print,'ALLFRAME_COMBINE_ORIG ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return
endif

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
  printlog,logf,'SCRIPTSDIR NOT INPUT'
  error = 'SCRIPTSDIR NOT INPUT'
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
           'default.sex','default.param','default.nnw','default.conv']
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
endfor ; scripts loop


printlog,logf,'Combining images in ',file

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

; Checking RAW file
rawtest = file_test(mchbase+'.raw')
if rawtest eq 0 then begin
  printlog,logf,mchbase+'.raw NOT FOUND'
  return
endif


;###########################################
; CHECK NECESSARY FILES

; Load the MCH file
LOADMCH,mchfile,files,trans

; Check that the fits, als, opt, and psf files exist
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

endfor

; FAKE, check that we have all the files that we need
if keyword_set(fake) then begin
  ; weights, scale, zero, comb_psf, _shift.mch
  chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_shift.mch']
  bdfiles = where(file_test(chkfiles) eq 0,nbdfiles)
  if nbdfiles gt 0 then begin
    error = 'FAKE.  Some necessary files not found. '+strjoin(chkfiles[bdfiles],' ')
    printlog,logf,error
    return
  endif
endif



;###########################################
; STEP 1: IMALIGN PREP

printlog,logf,''
printlog,logf,'Step A: Getting Weights'
printlog,logf,'-----------------------'

;-----------------------------------
; Computs Weights
if not keyword_set(fake) then begin
  ALLFRAME_GETWEIGHTS,mchfile,weights,scales,sky ;,raw2=raw2
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
  WRITECOL,scalefile,invscales,fmt='(F10.6)'  ; want to scale it UP
  zerofile = mchbase+'.zero'
  WRITECOL,zerofile,-sky,fmt='(F12.4)'  ; want to remove the background, set to 1st frame

; FAKE, use existing ones
endif else begin
  weightfile = mchbase+'.weights'
  scalefile = mchbase+'.scale'
  zerofile = mchbase+'.zero'
  READCOL,weightfile,weights,format='F',/silent
  READCOL,scalefile,invscales,format='F',/silent
  scales = 1.0/invscales
  READCOL,zerofile,sky,format='F',/silent
  sky = -sky
endelse
  

;---------------------------------------
; Get X/Y translations using DAOMASTER
;  NO ROTATION ONLY SHIFTS
;  Need these shifts for IMSHIFT
shiftmch = mchbase+'_shift'
if not keyword_set(fake) then begin
  print,'Measuring X/Y shifts'
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
endif
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
; STEP B: FIX BAD PIXELS
printlog,logf,''
printlog,logf,'Step B: Fixing bad pixels'
printlog,logf,'-------------------------'
FIXIMAGES,'@'+infile,satlevel=satlevel  ;6e4
; This also makes the FILE.mask.fits files for each image

; Find the maximum saturation level
satlevelarr = fltarr(nfiles)
for i=0,nfiles-1 do begin
  ;head = headfits(base[i]+'.fits')
  FITS_READ,base[i]+'.fits',im,head,/no_abort
  saturate = sxpar(head,'SATURATE',count=nsaturate,/silent)
  if nsaturate eq 0 then saturate=max(im)-1000.
  satlevelarr[i] = saturate
endfor
maxsatlevel = max(satlevelarr)


;###########################################
; STEP C: IMALIGN
;  This figures out the X/Y-shifts between the images
;  and creates the shifted images (".shft.fits")
printlog,logf,''
printlog,logf,'Step C: IMALIGN'
printlog,logf,'---------------'
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

;stop

; Trim the images
if keyword_set(trimcomb) then begin

  ; Calculate the trim section
  hd = headfits(reffile)
  xsize = lonarr(nfiles)+sxpar(hd,'NAXIS1',/silent)
  ysize = lonarr(nfiles)+sxpar(hd,'NAXIS2',/silent)
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
  xsize = sxpar(hd,'NAXIS1',/silent)
  ysize = sxpar(hd,'NAXIS2',/silent)
  trimsection = [1,xsize,1,ysize]
  xoff = 0
  yoff = 0
endelse

; Delete the temporary FITS files
FILE_DELETE,tempfits,/allow



;###########################################
; STEP D: MAKE BAD PIXEL/WEIGHT MAP
; in how many images does the pixel need to be bad??
; 1. shift the masks (created by FIXIMAGES.PRO)
;     using the shifts from IRAF_IMALIGN
; 2. trim the masks
; 3. combine the masks
printlog,logf,''
printlog,logf,'Step D: Making Bad Pixel Mask'
printlog,logf,'-----------------------------'

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
    ltv1 = sxpar(head,'LTV1',/silent)  ; 0 if not found
    ltv2 = sxpar(head,'LTV2',/silent)  ; 0 if not found
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
endfor
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
; STEP E: COMBINE IMAGES
printlog,logf,''
printlog,logf,'Step E: IMCOMBINE'
printlog,logf,'-----------------'

; SCALE the images for combining
;-------------------------------
if not keyword_set(nocmbimscale) then begin

  ; Put BPM mask names in file headers
  ;  these will be used by IMCOMBINE
  for i=0,nfiles-1 do begin
    head = headfits(outfiles[i])
    sxaddpar,head,'BPM',outmaskfiles[i]
    modfits,outfiles[i],0,head
  endfor

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
  rdnoise = rdnoise > 0.01   ; must be >=0.01 or it will be 0.00 in the opt file and daophot will crash
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

; Copy the original MCH file to COMB.MCH
FILE_COPY,file,mchdir+'/'+mchbase+'.comb.mch',/allow,/over

; CD back to the original directory
cd,curdir

BOMB:

if keyword_set(stp) then stop

end
