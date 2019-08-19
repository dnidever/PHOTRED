;+
;
; PHOTRED_WCS
;
; This gets the WCS for the images.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The WCS of the images is updated
;
; By D.Nidever  Feb 2008
;-

pro photred_wcs,redo=redo,stp=stp,testing=testing,std=std

COMMON photred,setup

print,''
print,'###################'
print,'RUNNING PHOTRED_WCS'
print,'###################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'WCS'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','writeline','readlist','printlog','fits_read','sxpar',$
         'sexig2ten','queryvizier','undefine','writecol','readcol','photred_loadsetup',$
         'photred_updatelists','push','undefine','printlog','photred_mkopt','mkopt','imfwhm',$
         'wcsfit','matchstars','wcsfit_imacs','mk_imacsmask','hdr2wcstnx','wcstnx_xy2rd',$
         'wcstnx_rd2xy','parsetnx','wcstnx2hdr','wcstnxcor','rd2xieta','xieta2rd',$
         'head_xyad','head_adxy','getpixscale','mktemp','loadcoo','loadaper','photred_getgain',$
         'photred_getrdnoise','roi_cut','srcmatch','array_indices2','first_el','mad',$
         'maxloc','minloc','mpfit','photred_getinput','range','readpar','scale','scale_vector',$
         'slope','stringize','strsplitter','add_tag','loadinput','randomize','rotsphcen','rotsph',$
         'touchzero','importascii','signs','strmult','strtrim0','combine_structs','stress',$
         'strep']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; LOAD THE SETUP FILE if not passed
;--------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count,std=std
  if count lt 1 then return
endif


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument
telescope = strlowcase(READPAR(setup,'TELESCOPE'))
instrument = strlowcase(READPAR(setup,'INSTRUMENT'))

; Reference Catalog Name
wcsrefname = READPAR(setup,'WCSREFNAME')
if wcsrefname ne 'USNO-B1' and wcsrefname ne '2MASS-PSC' and wcsrefname ne 'UCAC4' and $
   wcsrefname ne 'GAIA/GAIA' then wcsrefname='GAIA/GAIA'

; WCSFIT Search distance
searchdist = READPAR(setup,'SEARCHDIST')
if searchdist ne '0' and searchdist ne '' and searchdist ne '-1' then $
   searchdist=float(searchdist) else undefine,searchdist
if n_elements(searchdist) gt 0 then if searchdist lt 0. then undefine,searchdist

; Maximum WCS RMS value
wcsrmslim = READPAR(setup,'WCSRMSLIM')
if wcsrmslim ne '0' and wcsrmslim ne '' and wcsrmslim ne '-1' then $
   wcsrmslim=float(wcsrmslim) else undefine,wcsrmslim
if n_elements(wcsrmslim) gt 0 then if wcsrmslim lt 0. then undefine,wcsrmslim

; WCS catalog err limit
wcscaterrlim = READPAR(setup,'WCSCATERRLIM')
if wcscaterrlim ne '0' and wcscaterrlim ne '' and wcscaterrlim ne '-1' then $
   wcscaterrlim=float(wcscaterrlim) else undefine,wcscaterrlim

; Pixel scale, this is only used for non-standard setups
wcspixscale = READPAR(setup,'PIXSCALE')
if wcspixscale ne '0' and wcspixscale ne '' and wcspixscale ne '-1' then $
   wcspixscale=float(wcspixscale) else undefine,wcspixscale
if n_elements(wcspixscale) gt 0 then if wcspixscale lt 0. then undefine,wcspixscale

; WCSLEFT, this is only used for non-standard setups
wcsleft = READPAR(setup,'WCSLEFT')
gdleft = first_el(where(strpos(['N','S','E','W'],strupcase(wcsleft)) ne -1,ngdleft))
if ngdleft eq 0 then undefine,wcsleft
;if wcsleft eq '0' or wcsleft eq '' or wcsleft eq '-1' then undefine,wcsleft

; WCSUP, this is only used for non-standard setups
wcsup = READPAR(setup,'WCSUP')
gdup = first_el(where(strpos(['N','S','E','W'],strupcase(wcsup)) ne -1,ngdup))
if ngdup eq 0 then undefine,wcsup
;if wcsup eq '0' or wcsup eq '' or wcsup eq '-1' then undefine,wcsup

; SKIPWCS
skipwcs = READPAR(setup,'SKIPWCS')
skipwcs = strtrim(skipwcs,2)
if skipwcs eq '1' then skipwcs=1 else undefine,skipwcs

; Hyperthread?
hyperthread = READPAR(setup,'hyperthread')
if hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if strtrim(hyperthread,2) eq '0' then hyperthread=0
if n_elements(hyperthread) eq 0 then hyperthread=0

; Getting NMULTI
nmulti = READPAR(setup,'NMULTI')
nmulti = long(nmulti)

; Use NMULTI_WCS if set
nmultiwcs = READPAR(setup,'NMULTI_WCS')
if nmultiwcs ne '0' and nmultiwcs ne '' and nmultiwcs ne '-1' then nmulti=long(nmultiwcs)
nmulti = nmulti > 1  ; must be >=1

; SKIPCHECK, skip all of the detailed file checking
skipcheck = READPAR(setup,'SKIPCHECK',count=nskipcheck)
if nskipcheck eq 0 then undefine,skipcheck

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif

;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS or FITS.FZ files
; OUTLIST        FITS or FITS.FZ files
; SUCCESSLIST    FITS or FITS.FZ files


; Get input
;-----------
precursor = 'SPLIT'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext=['fits','fits.fz'])
ninputlines = lists.ninputlines

; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines



; Skipping WCS
;-------------
if keyword_set(skipwcs) then begin
  printlog,logfile,''
  printlog,logfile,'SKIPWCS=1  SKIPPING WCS FITTING. Putting all input files into the outlist'
  printlog,logfile,''

  ; All files are put into the outlist
  outlist = inputlines
  ; UPDATE the Lists
  PHOTRED_UPDATELISTS,lists,outlist=outlist,setupdir=curdir,/silent

  goto,FINISH

endif

;; GOOD IDEA BUT THIS IS WHAT THE INPUT/SUCCESS LISTS ARE FOR!!
;; Check for previous success
;if not keyword_set(redo) then begin
;  printlog,logfile,'Checking for previous successes'
;  ; Loop through the input files
;  done = intarr(ninputlines)
;  for i=0,ninputlines-1 do begin
;    longfile = inputlines[i]
;    file = FILE_BASENAME(longfile)
;    filedir = FILE_DIRNAME(longfile)
;    if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then fpack=1 else fpack=0
;    if fpack eq 1 then base=FILE_BASENAME(file,'.fits.fz') else $
;      base = FILE_BASENAME(file,'.fits')
;    if file_test(longfile) eq 0 then begin
;      done[i] = 1
;      goto,BOMB1
;    endif
;    ; Successful
;    if fpack eq 1 then head=PHOTRED_READFILE(longfile,exten=1,/header) else $
;       head = PHOTRED_READFILE(longfile,exten=0,/header)
;    ctype1 = strtrim(SXPAR(head,'CTYPE1',count=nctype1,/silent),2)
;    dum = where(stregex(head,'WCSFIT: RMS',/boolean) eq 1,nwcsfit)
;    if (nctype1 gt 0 and ctype1 ne '0' and nwcsfit gt 0) then done[i]=1
;    BOMB1:
;  endfor
;  gd = where(done eq 0,ngd,comp=bd,ncomp=nbd)
;  if ngd eq 0 then begin
;    printlog,logfile,'All files previously successful.  None to process.'
;    ; All files are put into the outlist
;    outlist = inputlines
;    ; UPDATE the Lists
;    PHOTRED_UPDATELISTS,lists,outlist=outlist,setupdir=curdir,/silent
;    goto,FINISH
;  endif
;  printlog,logfile,strtrim(nbd,2),' previous successes'
;  printlog,logfile,strtrim(ngd,2),' files still to run'
;  inputlines = inputlines[gd]
;  ninputlines = ngd
;endif


;####################################################
;#  PROCESSING THE FILES
;####################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''

; Initialzing some arrays
UNDEFINE,outlist,successlist,failurelist

; Figure out which files to make OPT files for
undefine,cmd,cmddir,cmdlongfile,cmdfile

; Loop through the input files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then fpack=1 else fpack=0
  if fpack eq 1 then base=FILE_BASENAME(file,'.fits.fz') else $
    base = FILE_BASENAME(file,'.fits')
  printlog,logfile,strtrim(i+1,2),' ',longfile

  ; Check that the file exists
  test = FILE_TEST(longfile)
  if (test eq 0) then begin
    printlog,logfile,longfile,' NOT FOUND'
    PUSH,failurelist,longfile
    goto,BOMB
  endif


  CD,filedir

  ;; Check on the files
  if not keyword_set(skipcheck) then begin

    ; Load the header
    head = PHOTRED_READFILE(longfile,/header)
    object = SXPAR(head,'OBJECT',/silent)

    ; Does this have multiple extensions
    rfile = filedir+'/.'+file
    if file_test(rfile) eq 0 then begin
      FITS_OPEN,file,fcb,message=message0
      nextend = fcb.nextend
      FITS_CLOSE,fcb
      if nextend gt 0 then mef=1 else mef=0
    endif else mef=0   ; resource files are not MEF

    ; We only do SPLIT files, NOT MEF files
    if (mef eq 1) then begin
      printlog,'This is a MEF file.  Need SPLIT files'
      PUSH,failurelist,longfile
      goto,BOMB
    endif

    ; Make sure that |BITPIX| > 16
    ; Only works for single chip images
    bitpix = long(SXPAR(head,'BITPIX',/silent))
    if (bitpix eq 8 or bitpix eq 16) and (mef eq 0) and (fpack eq 0) then begin
      printlog,logfile,'BIXPIX = ',strtrim(bitpix,2),'.  Making image FLOAT'

      ; Read in the image
      ; FITS_READ,file,im,head,/no_abort,message=message
      ;im = MRDFITS(file,0,head,status=status,/silent)
      im = PHOTRED_READFILE(file,head,error=error)

      ; Make sure BZERO=0
      bzero = sxpar(head,'BZERO',count=nbzero,/silent)
      if nbzero gt 0 then sxaddpar,head,'BZERO',0.0

      ; Write the FLOAT image
      if n_elements(error) eq 0 and size(im,/type) lt 4 then $
        FITS_WRITE_RESOURCE,file,float(im),head

      ; There was a problem reading the image
      if n_elements(error) gt 0 then begin
        printlog,logfile,'PROBLEM READING IN ',file,' ',error
        PUSH,failurelist,longfile
        goto,BOMB
      endif
    endif
  endif  ; not skipcheck


  ;---------------------------------------------------------------
  ; GOING THROUGH THE DIFFERENT INSTRUMENTS
  ;---------------------------------------------------------------



  ; MOSAIC IMAGES
  ;---------------
  if (instrument eq 'mosaic') then begin

    ; SPLIT MOSAIC frames have a default TNX WCS
    ; The distortion terms are normally okay
    printlog,logfile,''
    printlog,logfile,file,' is a SPLIT MOSAIC IMAGE'
    printlog,logfile,''

    ; Add to the PBS command list
    cmd1 = "WCSFIT,'"+file+"',up='E',left='S',refname='"+wcsrefname+"'"
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
    if n_elements(wcscaterrlim) gt 0 then cmd1+=",caterrlim="+strtrim(wcscaterrlim,2)
    if keyword_set(redo) then cmd1+=',/redo'
    PUSH,cmd,cmd1
    PUSH,cmddir,filedir
    PUSH,cmdlongfile,longfile
  Endif ; MOSAIC


  ; IMACS IMAGES
  ;-------------
  if (instrument eq 'imacs') then begin
    printlog,logfile,''
    printlog,logfile,file,' is an IMACS IMAGE'
    printlog,logfile,''

    ; Add to the PBS command list
    cmd1 = "WCSFIT_IMACS,'"+file+"',refname='"+wcsrefname+"'"
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
    if n_elements(wcscaterrlim) gt 0 then cmd1+=",caterrlim="+strtrim(wcscaterrlim,2)
    if keyword_set(redo) then cmd1+=',/redo'
    PUSH,cmd,cmd1
    PUSH,cmddir,filedir
    PUSH,cmdlongfile,longfile
  Endif ; IMACS


  ; LBC IMAGES
  ;------------
  if (instrument eq 'lbc') then begin
    printlog,logfile,''
    printlog,logfile,file,' is a SPLIT LBC IMAGE'
    printlog,logfile,''

    printlog,logfile,''
    printlog,logfile,'*************************************'
    printlog,logfile,'LBC NOT PROPERLY SET UP YET!!!'
    printlog,logfile,'STILL NEED TO DEAL WITH DISTORTIONS!!!'
    printlog,logfile,'*************************************'
    printlog,logfile,''

    ; Getting chip number, 1-4
    chip = first_el(strsplit(base,'_',/extract),/last)
    chip = long(chip)

    ; CHIP 4 files are missing CDELT1/CDELT2
    head = PHOTRED_READFILE(file,/header)
    SXADDPAR,head,'CDELT1',6.2220000e-05
    SXADDPAR,head,'CDELT2',-6.2220000e-05
    if file_test(rfile) eq 1 then begin
      info = file_info(file)
      FITS_WRITE_RESOURCE,file,0,head
    endif else begin
      MODFITS,file,0,head
    endelse
    printlog,logfile,'Adding CDELT1/CDELT2 parameters for ',file

    ; LBC images
    uparr = ['N','N','N','E']
    leftarr = ['E','E','E','S']
    up = uparr[chip-1]
    left = leftarr[chip-1]
    pixscale = 0.224

    ; Add to the PBS command list
    undefine,error
    cmd1 = "WCSFIT,'"+file+"',up="+up+",left="+left+",pixscale="+strtrim(pixscale,2)+",refname='"+wcsrefname+"'"
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
    if n_elements(wcscaterrlim) gt 0 then cmd1+=",caterrlim="+strtrim(wcscaterrlim,2)
    if keyword_set(redo) then cmd1+=',/redo'
    PUSH,cmd,cmd1
    PUSH,cmddir,filedir
    PUSH,cmdlongfile,longfile
  Endif ; LBC


  ; SWOPE IMAGES
  ;-------------
  if (telescope eq 'swope') then begin

    printlog,logfile,''
    printlog,logfile,file,' is a SWOPE IMAGE'
    printlog,logfile,''

    ; Swope images only have RA/DEC by default
    ; Swope defaults
    up = 'E'
    left = 'N'
    if n_elements(wcsup) gt 0 then begin
      printlog,logfile,'Using input UP=',wcsup
      up = wcsup
    endif
    if n_elements(wcsleft) gt 0 then begin
      printlog,logfile,'Using input LEFT=',wcsleft
      left = wcsleft
    endif
    if n_elements(wcspixscale) gt 0 then begin
      printlog,logfile,'Using input PIXSCALE=',strtrim(wcspixscale,2),' "/pixel'
      pixscale = wcspixscale
    endif else begin
      pixscale = 0.6955
    endelse

    ; Add to the PBS command list
    cmd1 = "WCSFIT,'"+file+"',up='"+up+"',left='"+left+"',pixscale="+strtrim(pixscale,2)+",refname='"+wcsrefname+"'"
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
    if n_elements(wcscaterrlim) gt 0 then cmd1+=",caterrlim="+strtrim(wcscaterrlim,2)
    if keyword_set(redo) then cmd1+=',/redo'
    PUSH,cmd,cmd1
    PUSH,cmddir,filedir
    PUSH,cmdlongfile,longfile
  Endif ; SWOPE


  ; Unknown image type
  ;-------------------
  if (instrument ne 'mosaic' and instrument ne 'imacs' and instrument ne 'lbc' and $
      telescope ne 'swope') then begin

    ;printlog,logfile,''
    ;printlog,logfile,file,' is NOT a MOSAIC/IMACS/LBC/SWOPE image.  Trying anyway.'
    ;printlog,logfile,''

    if n_elements(wcspixscale) gt 0 then pixscale=wcspixscale

    ; Pixscale NOT input
    if n_elements(pixscale) eq 0 then $
    GETPIXSCALE,file,pixscale
    if (pixscale lt 0.0) then begin
      printlog,'ERROR - Cannot get pixel scale for ',file
      PUSH,failurelist,longfile
      goto,BOMB
    endif

    ; Add to the PBS command list
    cmd1 = "WCSFIT,'"+file+"',pixscale="+strtrim(pixscale,2)+",refname='"+wcsrefname+"'"
    if n_elements(wcsleft) gt 0 then cmd1+=",left='"+wcsleft+"'"
    if n_elements(wcsup) gt 0 then cmd1+=",up='"+wcsup+"'"
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
    if n_elements(wcscaterrlim) gt 0 then cmd1+=",caterrlim="+strtrim(wcscaterrlim,2)
    if keyword_set(redo) then cmd1+=',/redo'
    PUSH,cmd,cmd1
    PUSH,cmddir,filedir
    PUSH,cmdlongfile,longfile
  Endif

  BOMB:

  ; Back to original directory
  CD,curdir

ENDFOR

; Run WCSFIT with PBS_DEAMON
;----------------------------
printlog,logfile,'Running WCSFIT on images'
printlog,logfile,systime(0)
ncmd = n_elements(cmd)
printlog,logfile,strtrim(ncmd,2),' files to run WCSFIT on'

; no update
;print,'NO UPDATE'
;cmd += ',/noupdate'

if ncmd gt 0 then begin
  cmd = "cd,'"+cmddir+"' & "+cmd  ; go to the directory
  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='wfit',hyperthread=hyperthread,/idle,waittime=1,scriptsdir=scriptsdir
endif

; Check for success/failures
for i=0,ncmd-1 do begin
  ; Successful
  if strmid(cmdlongfile[i],6,7,/reverse_offset) eq 'fits.fz' then head=PHOTRED_READFILE(cmdlongfile[i],exten=1,/header) else $
     head = PHOTRED_READFILE(cmdlongfile[i],exten=0,/header)
  ctype1 = strtrim(SXPAR(head,'CTYPE1',count=nctype1,/silent),2)
  dum = where(stregex(head,'WCSFIT: RMS',/boolean) eq 1,nwcsfit)
  if (nctype1 gt 0 and ctype1 ne '0' and nwcsfit gt 0) then begin
    PUSH,successlist,cmdlongfile[i]
    PUSH,outlist,cmdlongfile[i]
  
  ; Failure
  endif else begin
    printlog,logfile,'PROBLEMS FITTING WCS FOR ',file_basename(cmdlongfile[i])
    PUSH,failurelist,cmdlongfile[i]
  endelse
endfor

FINISH:

;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


printlog,logfile,'PHOTRED_WCS Finished  ',systime(0)

if keyword_set(stp) then stop

end
