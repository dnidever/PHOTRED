pro photred_astrom,redo=redo,stp=stp

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

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')




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
    fitsfile = base+'.fits'
    photfile = magfile
  Endif ; 'mag' file

  ; DAOPHOT output
  ;------------------
  If (ending eq 'mch') then begin
    base = FILE_BASENAME(file,'.mch')
    mchfile = file
    rawfile = base+'.raw'
    fitsfile = base+'.fits'
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
  ; This is copied from LOADRAW.PRO
  ;  We should really use LOADRAW.PRO, but I'm not sure
  ;  if it works for .mag files.

  ; Figure out the number of columns
  line1='' & line2='' & line3='' & line4='' & line5=''
  openr,unit,/get_lun,photfile
  readf,unit,line1
  readf,unit,line2
  readf,unit,line3

  ; First line for the first star
  readf,unit,line4
  arr4 = strsplit(line4,' ',/extract)
  narr4 = n_elements(arr4)
  ncol = narr4

  ; Check for continuation lines
  endflag = 0
  nstarline = 1
  WHILE (endflag ne 1) do begin

    line4 = ''
    readf,unit,line4

    ; If there are too many frames/columns then DAOMASTER puts
    ; these on separate lines and lead with ~27 blank spaces
    ;  Not sure if this is needed for MAG files as well

    ; This is a continuation line
    if strtrim(strmid(line4,0,15),2) eq '' then begin
      arr4 = strsplit(line4,' ',/extract)
      narr4 = n_elements(arr4)
      ncol += narr4
      nstarline++
    endif else endflag=1
  ENDWHILE
  close,unit
  free_lun,unit

  ; Stars in this file
  numstar = (FILE_LINES(photfile)-3L )/nstarline

  nextra = ncol - 2*numobs    ; nextra includes, ID, X, Y, CHI, SHARP, etc.

  ; mastable is where everything is stored, id, x, y, unsolved magnitudes, chi, sharp
  mastable = fltarr(numstar,2*numobs+nextra)

  ; Reading in the magnitude file
  get_lun,unit
  openr, unit, photfile

  line=''
  readf,unit,line
  readf,unit,line
  readf,unit,line

  ; Loop through the stars
  for j=0.,numstar-1 do begin
    instr=''
    instr1=''
    inline = fltarr(2*numobs+nextra)

    ; Loop through the lines per star
    for k=0l,nstarline-1 do begin
      readf, unit, instr1
      ;instr += instr1
      ; remove extra 25 spaces at the beginning of extra/wrap lines
      if k eq 0 then instr+=instr1 else instr+=strmid(instr1,25)
    endfor
    ; We need to use the formatted read because sometimes there are
    ; NO spaces between the numbers in the columns.
    ; ***BUT if the daomaster format changes then this will be MESSED UP!!!!***
    ; MAKEMAG.PRO uses a slightly different format than daomaster
    ;  more space for larger integers
    ;  fmt='(A1,I8,2F9.3,'+strtrim(nfiles*2+2,2)+'F9.4)'
    if ending eq 'mag' then fmt='(I9,2F9.3,'+strtrim(ncol-3,2)+'F9.4)' else $  ; makemag output
      fmt='(I7,2F9.3,'+strtrim(ncol-3,2)+'F9.4)'   ; daomaster output
    reads,instr,inline,format=fmt
    mastable[j,0:2*numobs+nextra-1] = inline[0:2*numobs+nextra-1]
  endfor

  close, unit
  free_lun,unit

  ; Now convert to PHOT structure
  dum = {id:0L,x:0.0,y:0.0}
  for m=0,numobs-1 do begin
    dum = CREATE_STRUCT(dum,'MAG'+strtrim(m+1,2),0.0,'MAG'+strtrim(m+1,2)+'ERR',0.0)
  end
  dum = CREATE_STRUCT(dum,'CHI',0.0,'SHARP',0.0)
  if ending eq 'mag' then dum = CREATE_STRUCT(dum,'FLAG',0L,'PROB',0.0)
  dum = CREATE_STRUCT(dum,'RA',0.0d0,'DEC',0.0d0)
  phot = REPLICATE(dum,numstar)
  ncopy = 2*numobs+nextra
  for j=0,ncopy-1 do phot.(j)=reform(mastable[*,j])
  
  nphot = n_elements(phot)
  printlog,logfile,'Nstars = '+strtrim(nphot,2)


  ; Load the FITS header
  head = HEADFITS(fitsfile)

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

  ; Output the structure to the AST file
  astfile = base+'.ast'
  printlog,logfile,'File with RA/DEC coordinates is: ',astfile
  PRINTSTR,phot,astfile

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
                    failurelist=failurelist,/silent

  ;stop

END


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'PHOTRED_ASTROM Finished  ',systime(0)

if keyword_set(stp) then stop




end
