;+
;
; PHOTRED_ASTROM
;
; This gets coordinates from the FITS WCS and puts them
; in the photometry files.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The final calibrated photometry with accurate astrometry
;
; By D.Nidever  Mar 2008
;-

pro photred_astrom,redo=redo,stp=stp

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_ASTROM'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'ASTROM'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','head_xyad','hdr2wcstnx','parsetnx','wcstnx_xy2rd',$
         'wcstnxcor','xieta2rd','add_tag','printstr','photred_getinput','photred_updatelists',$
         'photred_loadsetup','push','undefine','printlog','strsplitter','combine_structs','touchzero',$
         'writeline','first_el','mktemp','stress','strep','loadmch']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1
; Minimum number of detections
ndetmin = READPAR(setup,'NDETMIN')
if (ndetmin eq '-1' or ndetmin eq '0' or ndetmin eq '') then undefine,ndetmin else ndetmin=long(ndetmin)

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')

; Catalog format to use
catformat = READPAR(setup,'catformat')
if catformat eq '0' or catformat eq '' or catformat eq '-1' then catformat='ASCII'
if catformat ne 'ASCII' and catformat ne 'FITS' then catformat='ASCII'


;###################
; GETTING INPUTLIST
;###################
; INLIST         PHOT files
; OUTLIST        AST files
; SUCCESSLIST    PHOT files

; Get input
;-----------
; Precursors: rename > wcs > split > daophot > match > allframe >
; apcor > *astrom* > calib > combine > deredden > save
precursor = ['ALLFRAME','MATCH']
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo)
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines



;##################################################
;#  PROCESSING THE FILES
;##################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,systime(0)

undefine,outlist,successlist,failurelist

; Loop through the input PHOT files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)

  printlog,logfile,''
  printlog,logfile,'=========================================='
  printlog,logfile,'Getting coordinates for '+file
  printlog,logfile,'=========================================='


  CD,filedir

  ending = first_el(strsplit(longfile,'.',/extract),/last)

  ; Wrong input file
  if (ending ne 'mch' and ending ne 'mag') then begin
    printlog,logfile,file+' ERROR - Input files must end in ".mag" or ".mch"'
    PUSH,failurelist,longfile
    goto,BOMB
  endif

  ; ALLFRAME output
  ;------------------
  If (ending eq 'mag') then begin
    base = FILE_BASENAME(file,'.mag')
    mchfile = base+'.mch'
    magfile = base+'.mag'
    ; Switched to the stacked/comb file for the WCS, 10/23/16
    ;  in the original combination procedure these two files ahd
    ;  the identical WCS, in the new combination procedure they
    ;  are very different but the "reference frame" is the 
    ;  combined frame.
    ;fitsfile = base+'.fits'
    fitsfile = base+'_comb.fits'
    photfile = magfile
  Endif ; 'mag' file

  ; DAOPHOT output
  ;------------------
  If (ending eq 'mch') then begin
    base = FILE_BASENAME(file,'.mch')
    mchfile = file
    rawfile = base+'.raw'
    fitsfile = base+'.fits'
    if file_test(fitsfile) eq 0 then fitsfile=base+'.fits.fz'
    photfile = rawfile
  Endif

  ; Check that the MAG/RAW, MCH and FITS files exist
  ;-------------------------------------------------
  mchtest = FILE_TEST(mchfile)
  if mchtest eq 1 then nmchlines = FILE_LINES(mchfile) else nmchlines=0
  phottest = FILE_TEST(photfile)
  if phottest eq 1 then nphotlines = FILE_LINES(photfile) else nphotlines=0
  fitstest = FILE_TEST(fitsfile)
  if (nmchlines eq 0) or (nphotlines eq 0) or (fitstest eq 0) then begin
    PUSH,failurelist,longfile
    if mchtest eq 0 then printlog,logfile,mchfile+' NOT FOUND'
    if mchtest eq 1 and nmchlines eq 0 then printlog,logfile,mchfile+' HAS 0 LINES'
    if phottest eq 0 then printlog,logfile,photfile+' NOT FOUND'
    if phottest eq 1 and nphotlines eq 0 then printlog,logfile,photfile+' HAS 0 LINES'
    if fitstest eq 0 then printlog,logfile,fitsfile+' NOT FOUND'
    goto,BOMB
  endif


  ; Load the MCH file
  ;------------------
  LOADMCH,mchfile,alsfiles
  nalsfiles = n_elements(alsfiles)    
  numobs = nalsfiles

  ; Load the photometry file
  ;-------------------------  
  phot0 = PHOTRED_READFILE(photfile)
  nphot = n_elements(phot0)
  schema = phot0[0]
  struct_assign,{dum:''},schema
  schema = CREATE_STRUCT(schema,'RA',0.0d0,'DEC',0.0d0)
  phot = REPLICATE(schema,nphot)
  struct_assign,phot0,phot
  undefine,phot0
  printlog,logfile,'Nstars = '+strtrim(nphot,2)


  ; Load the FITS header
  if strmid(fitsfile,6,7,/reverse_offset) eq 'fits.fz' then begin
    head = PHOTRED_READFILE(fitsfile,exten=1,/header)
    ; Fix the NAXIS1/2 values in the header
    sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
    sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
  endif else begin
    head = PHOTRED_READFILE(fitsfile,/header)
  endelse

  ; Checking that the header has a WCS
  EXTAST,head,astr
  nastr = n_elements(astr)
  if (nastr eq 0) then begin
    PUSH,failurelist,longfile
    printlog,logfile,fitsfile,' DOES NOT HAVE A WCS'
    goto,BOMB
  endif

  ; Check that the structure has X/Y
  tags = TAG_NAMES(phot)
  xgd = where(tags eq 'X',nxgd)
  ygd = where(tags eq 'Y',nygd)
  if (nxgd eq 0) or (nygd eq 0) then begin
    PUSH,failurelist,longfile
    printlog,logfile,file,' DOES NOT HAVE X/Y COORDINATES'
    goto,BOMB
  endif


  ; Converting to IDL X/Y convention, starting at (0,0)
  ; DAOPHOT has X/Y start at (1,1)
  x = phot.x - 1.0
  y = phot.y - 1.0


  ; Get RA/DEC coordinates for X/Y
  HEAD_XYAD,head,x,y,ra,dec,/degree

  ;; Add the RA/DEC tags
  ;ragd = where(tags eq 'RA',nragd)
  ;if (nragd) eq 0 then ADD_TAG,phot,'RA',0.0d0,phot
  ;decgd = where(tags eq 'DEC',ndecgd)
  ;if (ndecgd) eq 0 then ADD_TAG,phot,'DEC',0.0d0,phot

  ; Put RA/DEC into the structure
  phot.ra = ra
  phot.dec = dec


  ; Apply minimum number of detections, NDETMIN
  if n_elements(ndetmin) gt 0 then begin
    printlog,logfile,'Applying NDETMIN = '+strtrim(ndetmin,2)
    ; ID X Y MAG1 MAG1ERR MAG2 MAG2ERR MAG3 MAG3ERR MAG4 MAG4ERR
    magind = where(stregex(tags,'^MAG',/boolean) eq 1 and stregex(tags,'ERR$',/boolean) eq 0,nmagind)
    printlog,logfile,strtrim(nmagind,2)+' magnitude columns'
    ndet = lonarr(n_elements(phot))
    for i=0,nmagind-1 do ndet += (phot.(magind[i]) lt 50)
    gddet = where(ndet ge ndetmin,ngddet)
    if ngddet gt 0 then begin
      printlog,logfile,'Keeping '+strtrim(ngddet,2)+' of '+strtrim(n_elements(phot),2)+' sources'
    endif else begin
      printlog,logfile,'NO SOURCES PASS.  Keeping the first source.'      
      gddet = 0
      ngddet = 1
    endelse
    phot = phot[gddet]
  endif
  
  
  ; Output the structure to the AST file
  astfile = base+'.ast'
  printlog,logfile,'File with RA/DEC coordinates is: ',astfile

  if catformat eq 'FITS' then begin
    MWRFITS,phot,astfile,/create,/silent
  endif else begin   ; ASCII
    PRINTSTR,phot,astfile,/silent
  endelse
    
  ; Check that the file AST file is there
  asttest = FILE_TEST(astfile)
  if (asttest eq 1) then begin
    PUSH,outlist,filedir+'/'+astfile
    PUSH,successlist,longfile
  endif else begin
    PUSH,failurelist,longfile
    printlog,logfile,astfile,' NOT FOUND'
  endelse


  BOMB:

  CD,curdir

  ;##########################################
  ;#  UPDATING LIST FILES
  ;##########################################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir,/silent

  ;stop

ENDFOR


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


printlog,logfile,'PHOTRED_ASTROM Finished  ',systime(0)

if keyword_set(stp) then stop




end
