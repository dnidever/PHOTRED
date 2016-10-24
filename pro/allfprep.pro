;+
;
; ALLFPREP
;
; IDL version of Jamie's allfprep.cl that runs SExtractor
; iteratively on stacked images and then runs ALLSTAR
; to get PSF sources.
;
; The final ALS file will be called file+'_allf.als'
; A file with all of the sources that SExtractor found
; that can be matched to the ALS file using IDs is called
; file+'_allf.sex'
;
; INPUTS:
;  file      Filename of the stacked images
;  xoff      The offset in X between the original and shifted images
;            xorig = xshift + xoff
;  yoff      The offset in Y between the original and shifted images
;            yorig = yshift + yoff
;  =maskfile The name of a mask/weight file to use for the combined image.
;  =maxiter  Maximum number of times to iterate.
;             By default, maxiter=2
;  =scriptsdir  Directory that contains the scripts.
;  =detectprog  The program to use for source detection.  The options
;                 are 'SEXTRACTOR' (the default) or 'DAOPHOT'.
;  =logfile  A logfile to print output to.
;  /stp      Stop at the end of the program
;
; OUTPUTS:
;  als       Final ALLSTAR structure
;  =error    The error, if there was one, else undefined.
;
; USAGE:
;  IDL>allfprep,'ccd1001_comb.fits',als,scriptsdir=scriptsdir,logfile=logfile
;
; By D. Nidever    February 2008 (copied from Jamie's allfprep.cl)
;-

pro allfprep,file,als,xoff,yoff,maxiter=maxiter,scriptsdir=scriptsdir,$
             detectprog=detectprog0,logfile=logfile,stp=stp,error=error,$
             maskfile=maskfile

;scriptsdir = '/net/home/dln5q/daophot/'
undefine,als,error

nfile = n_elements(file)
if nfile eq 0 then begin
  print,'Syntax - allfprep,file,als,xoff,yoff,maxiter=maxiter,scriptsdir=scriptsdir,logfile=logfile,stp=stp'
  return
endif

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1


; Make sure file exists
test = file_test(file)
if test eq 0 then begin
  print,file,' NOT FOUND'
  return
end


; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'ALLFPREP ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return
endif

; What program are we using for source detection?
;------------------------------------------------
if n_elements(detectprog0) gt 0 then detectprog=detectprog0  ; internal variable
if n_elements(detectprog) eq 0 then detectprog='sextractor'  ; the default
detectprog = strlowcase(detectprog)
if detectprog eq 'sex' then detectprog='sextractor'
if detectprog eq 'dao' then detectprog='daophot'
if detectprog ne 'sextractor' and detectprog ne 'daophot' then begin
  print,'DETECTPROG = ',detectprog,' NOT AN OPTION.  Using SEXTRACTOR'
  detectprog = 'sextractor'
endif


; Copy the scripts
;---------------------
; No scriptsdir
if n_elements(scriptsdir) eq 0 then begin
  printlog,logf,'SCRIPTSDIR NOT INPUT'
  error = 'SCRIPTSDIR NOT INPUT'
  return
endif
; Check if the scripts exist in the current directory
scripts = ['default.sex','default.param','default.nnw','default.conv']
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


; Check that the SEXTRACTOR program exists
SPAWN,'which sex',out,errout
sexfile = FILE_SEARCH(out,count=nsexfile)
if (nsexfile eq 0) then begin
  printlog,logf,'SEXTRACTOR PROGRAM NOT AVAILABLE'
  return
endif

; Check that the ALLSTAR program exists
SPAWN,'which allstar',out,errout
alsfile = FILE_SEARCH(out,count=nalsfile)
if (nalsfile eq 0) then begin
  printlog,logf,'ALLSTAR PROGRAM NOT AVAILABLE'
  return
endif

; Check that the DAOPHOT program exists
SPAWN,'which daophot',out,errout
daofile = FILE_SEARCH(out,count=ndaofile)
if (ndaofile eq 0) then begin
  printlog,logf,'DAOPHOT PROGRAM NOT AVAILABLE'
  return
endif



base = FILE_BASENAME(file,'.fits')

; Defaults
nmaxiter = n_elements(maxiter)
if nmaxiter eq 0 then maxiter=2
if n_elements(xoff) eq 0 then xoff=0.0
if n_elements(yoff) eq 0 then yoff=0.0

; Read the .opt file
READLINE,base+'.opt',optlines
optarr = strsplitter(optlines,' ',/extract)
g = where(stregex(optlines,'HI =',/boolean) eq 1,ng)
satlevel = optarr[2,g[0]]
g = where(stregex(optlines,'GA =',/boolean) eq 1,ng)
gain = optarr[2,g[0]]
g = where(stregex(optlines,'FW =',/boolean) eq 1,ng)
fwhm = optarr[2,g[0]]

; Get pixel scale
GETPIXSCALE,file,scale
if scale eq 99.9 then scale=0.5   ; default

; Filenames
origfile = file
subbase = base+'_sub'
subfile = subbase+'.fits'
catfile = subbase+'.cat'
sexfile = base+'_allf.sex'
FILE_COPY,file,subfile,/overwrite

; Make customized SEXTRACTOR file
;--------------------------------
if (detectprog eq 'sextractor') then begin
  sexconfigfile = base+'.sex'
  FILE_COPY,'default.sex',sexconfigfile,/overwrite,/allow
  READLINE,sexconfigfile,sexlines
  sexlines2 = sexlines
  ; CATALOG_NAME
  g = where(stregex(sexlines2,'^CATALOG_NAME',/boolean) eq 1,ng)
  ;catfile = base+'.cat'
  sexlines2[g[0]] = 'CATALOG_NAME    '+catfile+' # name of the output catalog'
  ; SATUR_LEVEL
  g = where(stregex(sexlines2,'^SATUR_LEVEL',/boolean) eq 1,ng)
  sexlines2[g[0]] = 'SATUR_LEVEL     '+satlevel+'         # level (in ADUs) at which arises saturation'
  ; GAIN
  g = where(stregex(sexlines2,'^GAIN',/boolean) eq 1,ng)
  sexlines2[g[0]] = 'GAIN            '+gain+'             # detector gain in e-/ADU.'
  ; PIXEL_SCALE
  g = where(stregex(sexlines2,'^PIXEL_SCALE',/boolean) eq 1,ng)
  sexlines2[g[0]] = 'PIXEL_SCALE     '+strtrim(scale,2)+'             # size of pixel in arcsec (0=use FITS WCS info).'
  ; SEEING_FWHM
  g = where(stregex(sexlines2,'^SEEING_FWHM',/boolean) eq 1,ng)
  fwhmas = float(fwhm)*float(scale)
  sexlines2[g[0]] = 'SEEING_FWHM     '+strtrim(fwhmas,2)+'            # stellar FWHM in arcsec'
  ; MASK/WEIGHT file
  if n_elements(maskfile) ne 0 then begin
    PUSH,sexlines2,'WEIGHT_IMAGE  '+maskfile
    PUSH,sexlines2,'WEIGHT_TYPE   MAP_WEIGHT'
  endif
  ; Write the file
  WRITELINE,sexconfigfile,sexlines2
endif


printlog,logf,'Using '+strupcase(detectprog)+' for Source Detection'

;Loop to find all the objects the first time, bright star problems
flag = 0
nals = 0
count = 1
WHILE (flag eq 0) do begin

  ;print,''
  ;print,'----------------------------'
  ;print,'ITERATION ',strtrim(count,2)
  ;print,'----------------------------'
  ;print,''
  printlog,logf,'--Iteration '+strtrim(count,2)+'--'

  nlastals = nals


  ;######################
  ; Detect New Sources
  ;######################

  ;----------------------------
  ; SExtractor Source Detection
  ;----------------------------
  if (detectprog eq 'sextractor') then begin

    ;-------------------------------------
    ; Initial run of SExtractor -- output sex.cat
    ;!sex allf.fits -c default.sex
    ; Run sextractor on star subtracted file
    printlog,logf,'Running SExtractor'
    FILE_DELETE,catfile,/allow               ; delete sextractor catalog file
    SPAWN,'sex '+subfile+' -c '+sexconfigfile,out,errout
    num = file_lines(catfile)
    printlog,logf,'SExtractor found '+strtrim(num,2)+' sources'

    ;stop

    ;-------------------------------------
    ; Load sextractor output file
    ; default.param specifies the output columns
    fields = ['ID','X','Y','MAG','ERR','FLAGS','STAR']
    sex = IMPORTASCII(catfile,fieldnames=fields,/noprint)
    nsex = n_elements(sex)

    ;-------------------------------------
    ; Make coordinate input file for ALLSTAR

    ; Get ALS header
    line1 ='' & line2=''
    OPENR,unit,/get_lun,base+'.als'
    readf,unit,line1
    readf,unit,line2
    CLOSE,unit
    FREE_LUN,unit
    head = [line1,line2,'']

    ; Concatenate with the ALS file to make a combined
    ; list of sources
    coofile = base+'_all.coo'
    FILE_DELETE,coofile,/allow    ; delete final coordinate file
    nals = n_elements(als) 
    ; ALS file exists, concatenate
    if nals gt 0 then begin

      maxid = max(als.id)         ; The last ALS ID

      nsex = n_elements(sex)
      ; Making ALS structure for new SEX sources
      dum = {ID:0L,X:0.0,Y:0.0,MAG:0.0,ERR:0.0,SKY:0.0,ITER:0.0,CHI:0.0,SHARP:0.0}
      sex2 = replicate(dum,nsex)
      sex2.id = sex.id + maxid    ; offset the IDs, sex IDs start at 1
      sex2.x = sex.x
      sex2.y = sex.y
      sex2.mag = sex.mag
      sex2.err = sex.err
      sex2.sky = median(als.sky)
      sex2.iter = 1
      sex2.chi = 1.0
      sex2.sharp = 0.0

      ; Add to final sextractor file
      ; This will be a file that has all sources SExtractor found
      ; with their IDs properly offset so they can be matched to the
      ; final ALS file
      sex.id = sex.id + maxid     ; offse the IDs
      ;tempfile = maketemp('temp','txt')
      tempfile = MKTEMP('temp')
      WRITECOL,tempfile,sex.id,sex.x,sex.y,sex.mag,sex.err,sex.flags,sex.star,$
               fmt='(I10,F11.3,F11.3,F9.4,F9.4,I4,F6.2)'
      SPAWN,'cat '+tempfile+' >> '+sexfile,out,errout
      FILE_DELETE,tempfile,/allow


      ; Concatenate
      all = [als,sex2]
      nall = n_elements(all)

      ; Write to file
      WRITEALS,coofile,all,alshead

      ;print,'Iter ',strtrim(count,2),':  ',strtrim(nall,2),' stars so far'

    ; First time, no ALS file yet
    endif else begin

      ; Copy SExtractor output to file of all SExtractor sources
      FILE_COPY,catfile,sexfile,/overwrite

      ; Use the sextractor output file
      FILE_COPY,catfile,coofile,/overwrite
      ; Prepend ALS header
      WRITELINE,coofile,head,/prepend

      ;print,'Iter ',strtrim(count,2),':  ',strtrim(nsex,2),' stars so far'

    endelse


  ;-------------------------
  ; DAOPHOT Source Detection
  ;-------------------------
  Endif else begin

    ; Remove the SExtractor file.  If it exists it's an old one.
    FILE_DELETE,sexfile,/allow

    ;-------------------------------------
    ; Run DAOPHOT/FIND and PHOT on star subtracted file
    printlog,logf,'Running DAOPHOT'

    ; Copy the OPT file from _comb.opt to _comb_sub.opt
    if FILE_TEST(subbase+'.opt') eq 0 then $
      FILE_COPY,base+'.opt',subbase+'.opt',/over,/allow

    ; Copy the .opt file daophot.opt 
    if FILE_TEST('daophot.opt') eq 0 then $
      FILE_COPY,subbase+'.opt','daophot.opt',/over,/allow
 
    ; Make temporary photo.opt file
    tphotofile = MKTEMP('photo')
    undefine,photlines
    push,photlines,'A1 = '+STRING(fwhm*3.0,format='(F7.4)')
    push,photlines,'IS = 45.0000'
    push,photlines,'OS = 50.0000'
    WRITELINE,tphotofile,photlines
 
    ; Make a temporary script to run FIND
    undefine,lines
    push,lines,'#!/bin/sh'
    ;push,lines,'daophot="/net/astro/bin/daophot"'
    push,lines,'export image=${1}'
    push,lines,'rm ${image}.temp.log      >& /dev/null'
    push,lines,'rm ${image}.temp.coo      >& /dev/null'
    push,lines,'rm ${image}.temp.ap       >& /dev/null'
    push,lines,'daophot << END_DAOPHOT >> ${image}.temp.log'
    push,lines,'OPTIONS'
    push,lines,'${image}.opt'
    push,lines,''
    push,lines,'ATTACH ${image}.fits'
    push,lines,'FIND'
    push,lines,'1,1'
    push,lines,'${image}.temp.coo'
    push,lines,'y'
    push,lines,'PHOTOMETRY'
    push,lines,FILE_BASENAME(tphotofile)
    push,lines,''
    push,lines,'${image}.temp.coo'
    push,lines,'${image}.temp.ap'
    push,lines,'EXIT'
    push,lines,'END_DAOPHOT'
    ;tempfile = maketemp('dao','.sh')
    tempfile = MKTEMP('dao')    ; absolute path
    WRITELINE,tempfile,lines
    FILE_CHMOD,tempfile,'755'o
 
    ; Run the program
    SPAWN,tempfile+' '+subbase,out,errout
    FILE_DELETE,tempfile    ; delete the temporary script

    ;-------------------------------------
    ; Load DAOPHOT FIND and PHOT output files
    LOADCOO,subbase+'.temp.coo',coo,coohead
    LOADAPER,subbase+'.temp.ap',aper,aperhead
    FILE_DELETE,[subbase+'.temp.coo',subbase+'.temp.ap'],/allow,/quiet
    ncoo = n_elements(coo)
    printlog,logf,'DAOPHOT found '+strtrim(ncoo,2)+' sources' 


    ; Get ALS header 
    line1 ='' & line2=''
    OPENR,unit,/get_lun,base+'.als'
    readf,unit,line1
    readf,unit,line2
    CLOSE,unit
    FREE_LUN,unit
    head = [line1,line2,'']    

    ; Concatenate with the ALS file to make a combined
    ; list of sources
    coofile = base+'_all.coo'
    FILE_DELETE,coofile,/allow    ; delete final coordinate file
    nals = n_elements(als)
    ; ALS file exists, concatenate
    if (nals gt 0) then begin
    
      maxid = max(als.id)         ; The last ALS ID

      ; Make a fake ALS structure
      ;--------------------------
      ndao = n_elements(coo)
      ; Making ALS structure for new DAOPHOT sources
      dum = {ID:0L,X:0.0,Y:0.0,MAG:0.0,ERR:0.0,SKY:0.0,ITER:0.0,CHI:0.0,SHARP:0.0}
      dao = replicate(dum,ndao)
      dao.id = coo.id + maxid    ; offset the IDs, sex IDs start at 1
      dao.x = coo.x
      dao.y = coo.y
      dao.mag = aper.mag[0]
      dao.err = aper.err[0]
      dao.sky = aper.sky
      dao.iter = 1 
      dao.chi = 1.0
      dao.sharp = 0.0

      ; Concatenate
      all = [als,dao]
      nall = n_elements(all)
    
      ; Write to file
      WRITEALS,coofile,all,alshead
  
      ;print,'Iter ',strtrim(count,2),':  ',strtrim(nall,2),' stars so far'

    ; First time, no ALS file yet
    endif else begin


      ; Make a fake ALS structure
      ;--------------------------
      ndao = n_elements(coo)
      ; Making ALS structure for new DAOPHOT sources
      dum = {ID:0L,X:0.0,Y:0.0,MAG:0.0,ERR:0.0,SKY:0.0,ITER:0.0,CHI:0.0,SHARP:0.0}
      dao = replicate(dum,ndao)
      dao.id = coo.id
      dao.x = coo.x
      dao.y = coo.y
      dao.mag = aper.mag[0]
      dao.err = aper.err[0]
      dao.sky = aper.sky
      dao.iter = 1 
      dao.chi = 1.0
      dao.sharp = 0.0

      ; Use the DAOPHOT FIND/PHOT data
      WRITEALS,coofile,dao,head
  
      ;print,'Iter ',strtrim(count,2),':  ',strtrim(nsex,2),' stars so far'
    
    endelse

    ;stop

  Endelse  ; DAOPHOT

  ;stop

  ;-----------------------------------------------------
  ; Run ALLSTAR on all sources found so far
  ; on original frame
  printlog,logf,'Running ALLSTAR'
  FILE_DELETE,subfile,/allow          ; delete subfile
  FILE_DELETE,subbase+'.als',/allow   ; delete als output file

  ; Make input file
  undefine,cmd
  push,cmd,'    '
  push,cmd,base+'.fits'       ; image file
  push,cmd,base+'.psf'        ; psf file
  push,cmd,coofile            ; coordinate file
  push,cmd,subbase+'.als'     ; new als file
  push,cmd,subbase+'.fits'    ; subtracted fits file
  ;cmdfile = maketemp('temp','.inp')
  cmdfile = MKTEMP('temp')
  WRITELINE,cmdfile,cmd
  SPAWN,'allstar < '+cmdfile,out,errout

  FILE_DELETE,cmdfile,/allow

  ; Load ALS file
  LOADALS,subbase+'.als',als,alshead
  nals = n_elements(als)
  printlog,logf,'ALLSTAR found '+strtrim(nals,2)+' sources'

  ; How many new stars
  nnew = nals-nlastals
  printlog,logf,strtrim(nnew,2)+' new stars found'
  if nnew lt 10 then flag=1
  if nnew lt round(nals*0.01) then flag=1  ; more than 1% of total
  if count ge maxiter then flag=1

  ; Increment counter
  count++

  ;stop


ENDWHILE

; The X/Y pixel coordinates should already be in the reference image
; coordinate system

; Need to offset the final X/Y pixel coordinates for XOFF/YOFF
printlog,logf,'Applying offsets: Xoff='+strtrim(xoff,2)+'  Yoff='+strtrim(yoff,2)
LOADALS,subbase+'.als',als,alshead
als.x = als.x + xoff
als.y = als.y + yoff
WRITEALS,base+'_allf.als',als,alshead
printlog,logf,'Final ALS file = '+base+'_allf.als'
if FILE_TEST(sexfile) eq 1 then $
  printlog,logf,'Final SExtractor file = '+sexfile

if keyword_set(stp) then stop

end
