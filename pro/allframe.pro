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
;  =tile          Information on the tiling to use for the combined
;                   image.  If not set then the "original" method is
;                   used.
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
;  /usecmn        Use the common sources file of the reference image.
;  /fake          Run for artificial star tests.
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


pro allframe,file,tile=tile,stp=stp,scriptsdir=scriptsdir,detectprog=detectprog,$
             error=error,logfile=logfile,finditer=finditer0,$
             irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,trimcomb=trimcomb,$
             usecmn=usecmn,fake=fake

COMMON photred,setup

undefine,error

; Not enough inputs
nfile = n_elements(file)
if (nfile eq 0) then begin
  print,'Syntax - allframe,file,tile=tile,scriptsdir=scriptsdir,finditer=finditer,satlevel=satlevel,'
  print,'                  detectprog=detectprog,nocmbimscale=nocmbimscale,error=error,logfile=logfile,'
  print,'                  irafdir=irafdir,trimcomb=trimcomb,usecmn=usecmn,fake=fake,stp=stp'
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
scripts = ['getpsf.sh','allstar.sh','photo.opt','apcor.opt','lstfilter','goodpsf.pro','allframe.opt',$
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

; Combination settings and inputs
if n_elements(tile) eq 0 then cmborig=1
if n_elements(tile) eq 1 then begin
  if size(tile,/type) ne 8 then begin
    printlog,logf,'TILE must be a structure'
    return
  endif
  if tag_exist(tile,'TYPE') eq 0then begin
    printlog,logf,'TILE must have a TYPE column'
    return
  endif
  if tile.type eq 'ORIG' then cmborig=1
endif


printlog,logf,''
printlog,logf,'====================================='
printlog,logf,'RUNNING ALLFRAME on ',file
printlog,logf,'====================================='
printlog,logf,''
printlog,logf,systime(0)

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
fpack = bytarr(nfiles)
for i=0,nfiles-1 do begin
  ;dir = file_dirname(mchfile)
  base = file_basename(files[i],'.als')
  
  ; Checking FITS file
  fitstest = file_test(base+'.fits') or file_test(base+'.fits.fz')
  if fitstest eq 0 then begin
    printlog,logf,base+'.fits/.fits.fz NOT FOUND'
    return
  endif

  ; Uncompress FPACK FITS files if necessary, temporarily
  if file_test(base+'.fits') eq 0 and file_test(base+'.fits.fz') eq 1 then begin
    fpack[i] = 1
    printlog,logf,'Temporarily uncomporessing '+base+'.fits.fz'
    spawn,['funpack',base+'.fits.fz'],/noshell
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

  ; Checking PSF file
  psftest = file_test(base+'.psf')
  if psftest eq 0 then begin
    printlog,logf,base+'.psf NOT FOUND'
    return
  endif

  ; REMOVE ALF if it exists
  if FILE_TEST(base+'.alf') then FILE_DELETE,base+'.alf',/allow,/quiet

endfor

; REMOVE the .mag file if it exists
if FILE_TEST(mchbase+'.mag') then FILE_DELETE,mchbase+'.mag',/allow,/quiet


; FAKE, check that we have all the files that we need
if keyword_set(fake) then begin
  ; weights, scale, zero, comb_psf, _shift.mch
  if keyword_set(cmborig) then $
    chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_shift.mch'] else $
    chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_comb.mch']
  bdfiles = where(file_test(chkfiles) eq 0,nbdfiles)
  if nbdfiles gt 0 then begin
    error = 'FAKE.  Some necessary files not found. '+strjoin(chkfiles[bdfiles],' ')
    printlog,logf,error
    return
  endif
endif


;###########################################
; STEP 1: COMBINE the images

printlog,logf,'--------------------------'
printlog,logf,'STEP 1: Combine the images'
printlog,logf,'--------------------------'
printlog,logf,systime(0)

; Use the original combine code
if keyword_set(cmborig) then begin
  ALLFRAME_COMBINE_ORIG,file,fake=fake,scriptsdir=scriptsdir,error=error,logfile=logfile,$
               irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,trimcomb=trimcomb,$
               maskdatalevel=maskdatalevel,xoff=xoff,yoff=yoff
  combmch = mchbase+'.mch'
; New combine code
endif else begin
  ALLFRAME_COMBINE,file,tile=tile,fake=fake,scriptsdir=scriptsdir,error=error,logfile=logfile,$
               irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,$
               maskdatalevel=maskdatalevel,filestr=filestr
  xoff = 0.0
  yoff = 0.0
  combmch = mchbase+'_comb.mch'
endelse
;  There was an error in combination
if n_elements(error) gt 0 then begin
  printlog,logf,error
  return
endif
combfile = mchbase+'_comb.fits'
combweightfile = mchbase+'_comb.mask.fits'


;###########################################
; STEP 2: Get PSF for Combined Image
printlog,logf,'----------------------------------------'
printlog,logf,'STEP 2: Getting PSF for Combined Image'
printlog,logf,'----------------------------------------'
printlog,logf,systime(0)
if not keyword_set(fake) then begin
  ; Make .opt files, set saturation just below the mask data level
  MKOPT,combfile,satlevel=maskdatalevel-1000
  ; THIS IS NOW DONE IN ALLFRAME_COMBINE/ALLFRAME_COMBINE_ORIG.PRO ABOVE
  ;; Using CMN.LST of reference frame if it exists
  ;if file_test(mchbase+'.cmn.lst') and keyword_set(usecmn) then begin
  ;  print,'Using reference image COMMON SOURCE file'
  ;  file_copy,mchbase+'.cmn.lst',mchbase+'_comb.cmn.lst',/over,/allow
  ;endif
  ; Get the PSF of the combined image
  SPAWN,'./getpsf.sh '+file_basename(combfile,'.fits')

  ; If getpsf failed, lower the spatial PSF variations to linear
  if file_test(file_basename(combfile,'.fits')+'.psf') eq 0 then begin
    printlog,logf,'getpsf.sh failed.  Lowering spatial PSF variations to linear.  Trying again.'
    MKOPT,combfile,satlevel=maskdatalevel-1000,va=1
    SPAWN,'./getpsf.sh '+file_basename(combfile,'.fits')
  stop
  endif

; FAKE, use existing comb.psf file
endif else begin
  printlog,logf,'Using existing '+file_basename(combfile,'.fits')+'.psf file and running ALLSTAR.'
  SPAWN,'./allstar.sh '+file_basename(combfile,'.fits')
  printlog,logf,' '
endelse


;###########################################
; STEP 3: Run allframe prep
;  This iteratively runs SExtractor on the combined image
;  This can take a while.
printlog,logf,'--------------------------------'
printlog,logf,'STEP 3: Running allframe prep'
printlog,logf,'--------------------------------'
printlog,logf,systime(0)
; Make sure we have an allstar.opt file
if file_test('allstar.opt') eq 0 then file_copy,base[0]+'.als.opt','allstar.opt'

ALLFPREP,combfile,als,xoff,yoff,logfile=logfile,error=error,$
         detectprog=detectprog,scriptsdir=scriptsdir,maxiter=finditer,$
         maskfile=combweightfile
if n_elements(error) gt 0 then goto,BOMB


;###########################################
; STEP 4: Running ALLFRAME
printlog,logf,'----------------------------'
printlog,logf,'STEP 4: Running ALLFRAME'
printlog,logf,'----------------------------'
printlog,logf,systime(0)

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
FILE_DELETE,mchbase+'.bck',/allow
FILE_DELETE,base+'.bck',/allow
FILE_DELETE,base+'.nmg',/allow
FILE_DELETE,base+'.tfr',/allow   ; This gets overwritten!!!
FILE_DELETE,base+'j.fits',/allow
FILE_DELETE,base+'k.fits',/allow
FILE_DELETE,base+'.alf',/allow

; Make input file
undefine,cmd
push,cmd,'    '
push,cmd,combmch               ; mch file
push,cmd,mchbase+'_comb_allf.als'  ; coord file
push,cmd,'    '
;cmdfile = maketemp('temp','.inp')
cmdfile = MKTEMP('temp')
WRITELINE,cmdfile,cmd

SPAWN,'allframe < '+cmdfile

FILE_DELETE,cmdfile,/allow
FILE_DELETE,file_basename(files,'.als')+'j.fits',/allow   ; delete subtracted images



;###########################################
; STEP 5: Combine photometry with MAKEMAG
; This combines the photometry in the N alf files
; and averages chi and sharp
printlog,logf,'--------------------------'
printlog,logf,'STEP 5: Running MAKEMAG'
printlog,logf,'--------------------------'
printlog,logf,systime(0)

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
; With the new combined files we are using _comb.mch
;   and _comb.tfr, 10/23/16
; combmch = FILEBASE.mch        ORIG
; combmch = FILEBASE_comb.mch   NEW
; The tfr file will have the same name but with .tfr
MAKEMAG,file_basename(combmch,'.mch')+'.tfr',magfile

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
; STEP 6: Adding SExtractor information
printlog,logf,'----------------------------------------'
printlog,logf,'STEP 6: Adding SExtractor information'
printlog,logf,'----------------------------------------'
printlog,logf,systime(0)

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
printlog,logf,systime(0)

; Delete temporarily funpacked files
bdfpack = where(fpack eq 1,nbdfpack)
if nbdfpack gt 0 then FILE_DELETE,base[bdfpack]+'.fits',/allow

BOMB:

if keyword_set(stp) then stop

end
