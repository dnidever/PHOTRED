pro photred_wcs,redo=redo,stp=stp,testing=testing,std=std

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
if wcsrefname ne 'USNO-B1' and wcsrefname ne '2MASS-PSC' then undefine,wcsrefname

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




;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS files
; OUTLIST        FITS files
; SUCCESSLIST    FITS files


; Get input
;-----------
;precursor = ['SPLIT','RENAME']
precursor = 'SPLIT'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='fits')
ninputlines = lists.ninputlines

;if (ninputlines eq 0) then begin
;  precursor = 'RENAME'
;  lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='fits')
;  ninputlines = lists.ninputlines
;endif

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
  PHOTRED_UPDATELISTS,lists,outlist=outlist,/silent

  goto,FINISH

endif



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
  base = FILE_BASENAME(file,'.fits')

  ; Check that the file exists
  test = FILE_TEST(longfile)
  if (test eq 0) then begin
    printlog,logfile,longfile,' NOT FOUND'
    PUSH,failurelist,longfile
    goto,BOMB
  endif


  CD,filedir

  ; Load the header
  head = HEADFITS(longfile)
  object = SXPAR(head,'OBJECT')

  ; Does this have multiple extensions
  UNDEFINE,im,hh
  FITS_READ,file,im,hh,exten=1,message=message,/no_abort
  if n_elements(im) gt 0 and strtrim(message,2) eq '' then mef=1 else mef=0

  ; We only do SPLIT files, NOT MEF files
  if (mef eq 1) then begin
    printlog,'This is a MEF file.  Need SPLIT files'
    PUSH,failurelist,longfile
    goto,BOMB
  endif



  ; Make sure that |BITPIX| > 16
  ; Only works for single chip images
  bitpix = long(SXPAR(head,'BITPIX'))
  if (bitpix eq 8 or bitpix eq 16) and (mef eq 0) then begin
    printlog,logfile,'BIXPIX = ',strtrim(bitpix,2),'.  Making image FLOAT'

    ; Read in the image
    FITS_READ,file,im,head,/no_abort,message=message

    ; Make sure BZERO=0
    bzero = sxpar(head,'BZERO',count=nbzero)
    if nbzero gt 0 then sxaddpar,head,'BZERO',0.0

    ; Write the FLOAT image
    if (message[0] eq '') then $
    FITS_WRITE,file,float(im),head

    ; There was a problem reading the image
    if (message[0] ne '') then begin
      printlog,logfile,'PROBLEM READING IN ',file
      PUSH,failurelist,longfile
      goto,BOMB
    endif

  endif



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
    cmd1 = "WCSFIT,'"+file+"',up='E',left='S',refname="+wcsrefname
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
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
    cmd1 = "WCSFIT_IMACS,'"+file+"',refname="+wcsrefname
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
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
    head = headfits(file)
    SXADDPAR,head,'CDELT1',6.2220000e-05
    SXADDPAR,head,'CDELT2',-6.2220000e-05
    MODFITS,file,0,head
    printlog,logfile,'Adding CDELT1/CDELT2 parameters for ',file

    ; LBC images
    uparr = ['N','N','N','E']
    leftarr = ['E','E','E','S']
    up = uparr[chip-1]
    left = leftarr[chip-1]
    pixscale = 0.224

    ; Add to the PBS command list
    undefine,error
    cmd1 = "WCSFIT,'"+file+"',up="+up+",left="+left+",pixscale="+strtrim(pixscale,2)+",refname="+wcsrefname
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
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
    cmd1 = "WCSFIT,'"+file+"',up="+up+",left="+left+",pixscale="+strtrim(pixscale,2)+",refname="+wcsrefname
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
    if keyword_set(redo) then cmd1+=',/redo'
    PUSH,cmd,cmd1
    PUSH,cmddir,filedir
    PUSH,cmdlongfile,longfile
  Endif ; SWOPE


  ; Unknown image type
  ;-------------------
  if (instrument ne 'mosaic' and instrument ne 'imacs' and instrument ne 'lbc' and $
      telescope ne 'swope') then begin

    printlog,logfile,''
    printlog,logfile,file,' is NOT a MOSAIC/IMACS/LBC/SWOPE image.  Trying anyway.'
    printlog,logfile,''

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
    cmd1 = "WCSFIT,'"+file+"',pixscale="+strtrim(pixscale,2)+",refname="+wcsrefname
    if n_elements(wcsleft) gt 0 then cmd1+=",left="+wcsleft
    if n_elements(wcsup) gt 0 then cmd1+=",up="+wcsup
    if n_elements(searchdist) gt 0 then cmd1+=",searchdist="+strtrim(searchdist,2)
    if n_elements(wcsrmslim) gt 0 then cmd1+=",rmslim="+strtrim(wcsrmslim,2)
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
ncmd = n_elements(cmd)
printlog,logfile,strtrim(ncmd,2),' files to run WCSFIT on'

if ncmd gt 0 then begin
  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='wfit',hyperthread=hyperthread,/idle,waittime=5
endif


; Check for success/failures
for i=0,ncmd-1 do begin
  ; Successful
  head = HEADFITS(cmdlongfile[i])
  ctype1 = strtrim(SXPAR(head,'CTYPE1'),2)
  if (n_elements(error) eq 0 and ctype1 ne '0') then begin
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
                    failurelist=failurelist


printlog,logfile,'PHOTRED_WCS Finished  ',systime(0)

if keyword_set(stp) then stop

end
