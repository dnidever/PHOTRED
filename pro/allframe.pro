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
;  =setupdir      The original base directory which contains photred.setup.
;  =scriptsdir    The directory that contains all of the necessary scripts.
;  =irafdir       The IRAF home directory.
;  =logfile       A logfile to print to output to.
;  /usecmn        Use the common sources file of the reference image.
;  /fake          Run for artificial star tests.
;  =catformat     Catalog format to use: FITS or ASCII.  Default is ASCII.
;  =imager        Imager structure with basic information.
;  =workdir       Use a temporary working directory with this as the base.
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


pro allframe,file,tile=tile,setupdir=setupdir,scriptsdir=scriptsdir,detectprog=detectprog,$
             error=error,logfile=logfile,finditer=finditer0,$
             irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,trimcomb=trimcomb,$
             usecmn=usecmn,fake=fake,catformat=catformat,imager=imager,workdir=workdir,stp=stp

COMMON photred,setup

undefine,error

; Not enough inputs
nfile = n_elements(file)
if (nfile eq 0) or n_elements(setupdir) eq 0 then begin
  print,'Syntax - allframe,file,tile=tile,setupdir=setupdir,scriptsdir=scriptsdir,finditer=finditer,satlevel=satlevel,'
  print,'                  detectprog=detectprog,nocmbimscale=nocmbimscale,error=error,logfile=logfile,'
  print,'                  irafdir=irafdir,trimcomb=trimcomb,usecmn=usecmn,fake=fake,catformat=catformat,'
  print,'                  imager=imager,workdir=workdir,stp=stp'
  error = 'Not enough inputs'
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
   error = !ERROR_STATE.MSG
   PHOTRED_ERRORMSG,logfile=logf
   CATCH, /CANCEL 
   return
endif

;; Get the setup information
nsetup = n_elements(setup)
if nsetup eq 0 then begin
  PHOTRED_LOADSETUP,setup,setupdir=setupdir,count=count
  if count lt 1 then return
endif

; How many FIND iterations
if n_elements(finditer0) eq 0 then finditer=2 else finditer=finditer0
finditer = finditer < 10  ; maximum 10.

; Saturation level
if n_elements(satlevel) eq 0 then satlevel=6e4

; Scaling of the images to be combined
if n_elements(nocmbimscale) eq 0 then nocmbimscale=0

; Getting scripts directory and iraf directory
scriptsdir = READPAR(setup,'SCRIPTSDIR')
irafdir = READPAR(setup,'IRAFDIR')

; Catalog format
if n_elements(catformat) eq 0 then catformat='ASCII'

; No irafdir
if n_elements(scriptsdir) eq 0 then begin
  error = 'SCRIPTSDIR NOT INPUT'
  printlog,logf,error
  return
endif

; Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL
;---------------------------------------------------------------
CHECK_IRAF,iraftest,irafdir=irafdir
if iraftest eq 0 then begin
  error = 'IRAF TEST FAILED.  EXITING'
  print,error
  return
endif

; No scriptsdir
if n_elements(scriptsdir) eq 0 then begin
  error = 'SCRIPTSDIR NOT INPUT'
  printlog,logf,error
  return
endif
; Check if the scripts exist in the current directory
scripts = ['getpsfnofind.sh','allstar.sh','photo.opt','apcor.opt','lstfilter.f','goodpsf.pro','allframe.opt',$
           'default.sex','default.param','default.nnw','default.conv']
nscripts = n_elements(scripts)
; Loop through the scripts
for i=0,nscripts-1 do begin
  info = FILE_INFO(scriptsdir+'/'+scripts[i])
  curinfo = FILE_INFO(scripts[i])
  ; No file
  if info.exists eq 0 or info.size eq 0 then begin
    error = scriptsdir+'/'+scripts[i]+' NOT FOUND or EMPTY'
    printlog,logf,error
    return
  endif
  ; Check if the two files are the same size, if not copy it
  if info.size ne curinfo.size then begin
    FILE_COPY,info.name,curinfo.name,/overwrite
  endif
endfor ; scripts loop
;; Compile lstfilter.f it it wasn't already compiled
;;   compiling it locally allows for different architectures/machines
;;   using the same repository
;if file_test('lstfilter') eq 0 and file_test(scriptsdir+'/lstfilter.f') then begin
  printlog,logfile,'Compiling lstfilter.f'
  file_copy,scriptsdir+'/lstfilter.f','.',/over
  ;; Check which fortran compiler we have
  compiler = ''
  spawn,['which','gfortran'],out,errout,/noshell
  if file_test(strtrim(out[0],2)) eq 1 and errout[0] eq '' then compiler='gfortran'
  spawn,['which','g77'],out,errout,/noshell
  if file_test(strtrim(out[0],2)) eq 1 and errout[0] eq '' then compiler='g77'
  if compiler eq '' then begin
    printlog,logfile,'NO fortran compiler found'
    return
  endif
  ;; Compile
  if file_test('lstfilter') eq 1 then file_delete,'lstfilter'
  spawn,[compiler,'lstfilter.f','-o','lstfilter'],out,errout,/noshell
  if file_test('lstfilter') eq 0 or errout[0] ne '' then begin
    printlog,logfile,'ERROR in compiling lstfilter.f'
    return
  endif
;endif


; Check that the ALLFRAME program exists
SPAWN,['which','allframe'],out,errout,/noshell
allframefile = FILE_SEARCH(out,count=nallframefile)
;allframefile = FILE_SEARCH('/net/halo/bin/allframe.2008',count=nallframefile)
;allframefile = FILE_SEARCH('/net/halo/bin/allframe.2004.fixed',count=nallframefile)
if (nallframefile eq 0) then begin
  ;printlog,logf,'/net/halo/bin/allframe.2004.fixed NOT FOUND'
  ;printlog,logf,'/net/halo/bin/allframe.2008 NOT FOUND'
  error = allframefile+'NOT FOUND'
  printlog,logf,error
  return
endif

; Combination settings and inputs
if n_elements(tile) eq 0 then cmborig=1
if n_elements(tile) eq 1 then begin
  if size(tile,/type) ne 8 then begin
    error = 'TILE must be a structure'
    printlog,logf,error
    return
  endif
  if tag_exist(tile,'TYPE') eq 0then begin
    error = 'TILE must have a TYPE column'
    printlog,logf,error
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
  error = mchfile+' NOT FOUND'
  printlog,logf,error
  return
endif

; Checking RAW file
rawtest = file_test(mchbase+'.raw')
if rawtest eq 0 then begin
  error = mchbase+'.raw NOT FOUND'
  printlog,logf,error
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
    error = base+'.fits/.fits.fz NOT FOUND'
    printlog,logf,error
    return
  endif

  ; Uncompress FPACK FITS files if necessary, temporarily
  if file_test(base+'.fits') eq 0 and file_test(base+'.fits.fz') eq 1 then begin
    fpack[i] = 1
    printlog,logf,'Temporarily uncompressing '+base+'.fits.fz'
    spawn,['funpack',base+'.fits.fz'],/noshell
  endif

  ; Checking OPT file
  opttest = file_test(base+'.opt')
  if opttest eq 0 then begin
    error = base+'.opt NOT FOUND'
    printlog,logf,error
    return
  endif

  ; Checking ALS.OPT file
  alsopttest = file_test(base+'.als.opt')
  if alsopttest eq 0 then begin
    error = base+'.als.opt NOT FOUND'
    printlog,logf,error
    return
  endif

  ; Checking AP file
  aptest = file_test(base+'.ap')
  if aptest eq 0 then begin
    error = base+'.ap NOT FOUND'
    printlog,logf,error
    return
  endif

  ; Checking ALS file
  alstest = file_test(base+'.als')
  if alstest eq 0 then begin
    error = base+'.als NOT FOUND'
    printlog,logf,error
    return
  endif

  ; Checking LOG file
  logtest = file_test(base+'.log')
  if logtest eq 0 then begin
    error = base+'.log NOT FOUND'
    printlog,logf,error
    return
  endif

  ; Checking PSF file
  psftest = file_test(base+'.psf')
  if psftest eq 0 then begin
    error = base+'.psf NOT FOUND'
    printlog,logf,error
    return
  endif

  ; REMOVE ALF if it exists
  if FILE_TEST(base+'.alf') then FILE_DELETE,base+'.alf',/allow,/quiet

endfor

; REMOVE the .mag file if it exists
if FILE_TEST(mchbase+'.mag') then FILE_DELETE,mchbase+'.mag',/allow,/quiet


; FAKE, check that we have all the files that we need
if keyword_set(fake) then begin
  ; In early versions of ALLFRAME _comb.mch was called _shift.mch
  ;  Use new name with link
  if file_test(mchbase+'_comb.mch') eq 0 and file_test(mchbase+'_shift.mch') eq 1 then begin
    print,'Linking '+mchbase+'_comb.mch to '+mchbase+'_shift.mch'
    FILE_LINK,mchbase+'_shift.mch',mchbase+'_comb.mch' 
  endif
  ; weights, scale, zero, comb_psf, _comb.mch
  chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_comb.mch']
  bdfiles = where(file_test(chkfiles) eq 0,nbdfiles)
  if nbdfiles gt 0 then begin
    error = 'FAKE.  Some necessary files not found. '+strjoin(chkfiles[bdfiles],' ')
    printlog,logf,error
    return
  endif
endif


;;------------------------------------
;; Using a temporary working directory
;;------------------------------------
if n_elements(workdir) gt 0 then begin
  ;; Create a temporary directory in WORKDIR
  if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
  tempdir = first_el(MKTEMP('alf',outdir=workdir,/directory))
  FILE_CHMOD,tempdir,/a_execute
  printlog,logf,'Working in temporary directory ',tempdir
  ;; Copy over the files that we need
  ;;  this will copy the contents of symlinks
  FILE_COPY,mchbase+['.mch','.raw','.tfr'],tempdir
  for i=0,nfiles-1 do FILE_COPY,file_basename(files[i],'.als')+'.'+['fits','opt','als.opt','ap','als','log','psf'],tempdir
  ;; Copy resources files and headers if they exist
  for i=0,nfiles-1 do begin
    base1 = file_basename(files[i],'.als')
    if file_test('.'+base1+'.fits') then FILE_COPY,'.'+base1+'.fits',tempdir
    if file_test(base1+'.fits.head') then FILE_COPY,base1+'.fits.head',tempdir
  endfor
  ;; Copy files for FAKE
  if keyword_set(fake) then FILE_COPY,mchbase+['.weights','.scale','.zero','_comb.psf','_comb.mch'],tempdir
  ;; Copy the scripts
  FILE_COPY,scripts,tempdir
  ;; Go there
  CD,tempdir
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
               maskdatalevel=maskdatalevel,filestr=filestr,imager=imager
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
combbase = file_basename(combfile,'.fits')
if not keyword_set(fake) then begin
  ; Make .opt files, set saturation just below the mask data level
  PHOTRED_MKOPT,combfile,va=1,hilimit=maskdatalevel-1000,error=opterror
  if n_elements(opterror) gt 0 then begin
    printlog,logf,opterror
    return
  endif
  ;MKOPT,combfile,satlevel=maskdatalevel-1000
  ; THIS IS NOW DONE IN ALLFRAME_COMBINE/ALLFRAME_COMBINE_ORIG.PRO ABOVE
  ;; Using CMN.LST of reference frame if it exists
  ;if file_test(mchbase+'.cmn.lst') and keyword_set(usecmn) then begin
  ;  print,'Using reference image COMMON SOURCE file'
  ;  file_copy,mchbase+'.cmn.lst',mchbase+'_comb.cmn.lst',/over,/allow
  ;endif
  PHOTRED_GETPSF,combbase,error=error
  if n_elements(error) gt 0 then return

; FAKE, use existing comb.psf file
endif else begin
  printlog,logf,'Using existing '+combbase+'.psf file and running ALLSTAR.'
  SPAWN,'./allstar.sh '+combbase
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

; Sometimes the filenames get too long for allframe,
; use temporary files and symlinks
tbase = (file_basename(MKTEMP('allf',/nodot)))[0]  ; create base, leave so other processes won't take it
tmch = tbase+'.mch'  &  file_delete,tmch,/allow  &  file_link,combmch,tmch
tals = tbase+'.als'  &  file_delete,tals,/allow  &  file_link,mchbase+'_comb_allf.als',tals

; Make input file
undefine,cmd
push,cmd,'    '
;push,cmd,combmch               ; mch file
;push,cmd,mchbase+'_comb_allf.als'  ; coord file
push,cmd,tmch
push,cmd,tals
push,cmd,'    '
;cmdfile = maketemp('temp','.inp')
cmdfile = MKTEMP('temp')
WRITELINE,cmdfile,cmd

SPAWN,'allframe < '+cmdfile

; Rename tfr and nmg files
FILE_MOVE,tbase+'.tfr',combbase+'.tfr',/over,/allow
FILE_MOVE,tbase+'.nmg',combbase+'.nmg',/over,/allow
;FILE_MOVE,tbase+'.tfr',mchbase+'.tfr',/over,/allow
;FILE_MOVE,tbase+'.nmg',mchbase+'.nmg',/over,/allow
FILE_DELETE,cmdfile,/allow
FILE_DELETE,file_basename(files,'.als')+'j.fits',/allow   ; delete subtracted images
FILE_DELETE,[tbase,tmch,tals],/allow   ; delete temporary files and links


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
MAKEMAG,combbase+'.tfr',magfile,error=magerror
if n_elements(magerror) gt 0 then goto,BOMB

; Prepend the ALF header to the makemag file
line1='' & line2='' & line3=''
openr,unit,/get_lun,file_basename(files[0],'.als')+'.alf'
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
sexfile = combbase+'_allf.sex'
if FILE_TEST(sexfile) eq 1 then begin

  ;-------------------------------------
  ; Load sextractor output file
  ; default.param specifies the output columns
  if file_isfits(sexfile) eq 0 then begin
    READLINE,'default.param',fields
    gd = where(strmid(fields,0,1) ne '#' and strtrim(fields,2) ne '',ngd)
    fields = fields[gd]
    sex = IMPORTASCII(sexfile,fieldnames=fields,/noprint)
  endif else begin
    sex = MRDFITS(sexfile,1,/silent)
    if n_tags(sex) eq 1 then sex=MRDFITS(sexfile,2,/silent)
  endelse
  nsex = n_elements(sex)

  ; Load the MAKEMAG file
  LOADRAW,mchbase+'.makemag',mag,alfhead
  nmag = n_elements(mag)

  ; Match them with IDs
  MATCH,mag.id,sex.number,ind1,ind2,count=nind

  ; Add SExtractor information to mag file
  sextags = tag_names(sex)
  ;add_tag,mag,'flag',0L,mag
  ;add_tag,mag,'prob',0.0,mag
  ;mag[ind1].flag = sex[ind2].flag
  ;mag[ind1].prob = sex[ind2].prob
  ;; New columns
  newcols = ['FLAGS','CLASS_STAR','MAG_AUTO','MAGERR_AUTO','BACKGROUND','THRESHOLD','ISOAREA_WORLD',$
             'A_WORLD','B_WORLD','THETA_WORLD','ELLIPTICITY','FWHM_WORLD']
  newname = ['FLAG','PROB','MAG_AUTO','MAGERR_AUTO','BACKGROUND','THRESHOLD','ISOAREA',$
             'ASEMI','BSEMI','THETA','ELLIPTICITY','FWHM']
  for k=0,n_elements(newcols)-1 do begin
    colind = where(sextags eq newcols[k],ncolind)
    if ncolind gt 0 then begin
      add_tag,mag,newname[k],fix('',type=size(sex[0].(colind),/type)),mag
      mag[ind1].(n_tags(mag)-1) = sex[ind2].(colind)
      ;; convert to arcsec
      if newcols[k] eq 'A_WORLD' or newcols[k] eq 'B_WORLD' then mag[ind1].(n_tags(mag)-1) *= 3600
    endif
  endfor   


  if nind lt nmag then printlog,logf,'DID NOT MATCH ALL THE STARS!'


  ; Write the final output file
  ;----------------------------
  finalfile = mchbase+'.mag'
  if catformat eq 'FITS' then begin
    MWRFITS,mag,finalfile,/create,/silent
  endif else begin  ; ASCII
    PRINTSTR,mag,finalfile,/silent
    ;; THIS IS TE OLD WAY OF SAVING THE INFORMATION WITH THE ALS HEADER
    ;; How many observations are there
    ;tags = tag_names(mag)
    ;ntags = n_elements(tags)
    ;magind = where(stregex(tags,'^MAG',/boolean) eq 1,nmagind)
    ;; Copy the structure to a string array, then print it out
    ;outarr = strarr(ntags,nmag)
    ;;fmtarr = '('+['I9','F9.3','F9.3',replicate('F9.4',nmagind*2),'F9.4','F9.4','I5','F7.2']+')'
    ;;outfmt = '(A9,2A9,'+strtrim(nmagind*2,2)+'A9,2A9,A5,A7)'
    ;fmtarr = '('+['I9','F9.3','F9.3',replicate('F9.4',nmagind*2),'F9.4','F9.4','I5','F7.2']+')'
    ;outfmt = '(A9,2A9,'+strtrim(nmagind*2,2)+'A9,2A9,A5,A7)'
    ;for i=0,ntags-1 do outarr[i,*] = STRING(mag.(i),format=fmtarr[i])
    ;openw,unit,/get_lun,finalfile
    ;printf,unit,format=outfmt,outarr
    ;close,unit
    ;free_lun,unit
    ;; Prepend the ALF header
    ;WRITELINE,finalfile,[alfhead,' '],/prepend
  endelse
    
; DAOPHOT
endif else begin

  ; No SExtractor information, just copy .makemag to .mag
  finalfile = mchbase+'.mag'
  if catformat eq 'FITS' then begin
    LOADRAW,mchbase+'.makemag',mag,alfhead
    MWRFITS,mag,finalfile,/create,/silent
  endif else begin  ; ASCII
    FILE_COPY,mchbase+'.makemag',finalfile,/allow,/over
  endelse
    
endelse

printlog,logf,'FINAL ALLFRAME file = ',finalfile
printlog,logf,systime(0)

BOMB:

; Delete temporarily funpacked files
bdfpack = where(fpack eq 1,nbdfpack)
if nbdfpack gt 0 then FILE_DELETE,file_basename(files[bdfpack],'.als')+'.fits',/allow

;;------------------------------------------------
;; Working in temporary directory, copy files back
;;------------------------------------------------
if n_elements(workdir) gt 0 then begin
  printlog,logf,'Copying files back to original directory'
  ;; Make sure all files are writeable
  files0 = file_search(tempdir+'/*',count=nfiles0,/match_initial_dot)
  if nfiles0 gt 0 then FILE_CHMOD,files0,'755'o
  ;; Delete some files
  FILE_DELETE,tempdir+'/'+mchbase+['.mch','.raw','.tfr'],/allow
  for i=0,nfiles-1 do FILE_DELETE,tempdir+'/'+file_basename(files[i],'.als')+'.'+['fits','fits.head','opt','als.opt','ap','als','log','psf'],/allow
  for i=0,nfiles-1 do FILE_DELETE,tempdir+'/.'+file_basename(files[i],'.als')+'.fits',/allow
  if keyword_set(fake) then FILE_DELETE,tempdir+'/'+mchbase+['.weights','.scale','.zero','_comb.psf','_comb.mch'],/allow
  ;; Copy files back
  files = file_search(tempdir+'/*',count=nfiles,/match_initial_dot)  
  FILE_COPY,files,mchdir,/allow,/over
  ;; Delete all temporary files
  FILE_DELETE,files,/allow
  ;; CD back
  CD,curdir
  ;; Delete temporary directory
  ;;  leave the base working directory
  FILE_DELETE,tempdir
endif

if keyword_set(stp) then stop

end
