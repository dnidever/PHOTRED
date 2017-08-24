;+
;
; STDRED_ASTROM
;
; This gets the coordinates for stars using the WCS
; in the FITS header.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The TOT files with RA/DEC coordinates
;
; By D.Nidever  May 2008
;-

pro stdred_astrom,redo=redo,stp=stp

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING STDRED_ASTROM'
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
printlog,logfile,'Starting STDRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','head_xyad','hdr2wcstnx','parsetnx','wcstnx_xy2rd',$
         'wcstnxcor','xieta2rd','add_tag','printstr','photred_getinput','photred_updatelists',$
         'photred_loadsetup','push','undefine','printlog','strsplitter','combine_structs','touchzero',$
         'writeline','first_el','mktemp','stress','strep']
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
  PHOTRED_LOADSETUP,setup,count=count,/std
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
precursor = 'DAOGROW'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='tot')
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
  base = FILE_BASENAME(file,'.tot')
  filedir = FILE_DIRNAME(longfile)
  fitsfile = base+'.fits'
  if file_test(fitsfile) eq 0 then fitsfile=base+'.fits.fz'

  printlog,logfile,''
  printlog,logfile,'=========================================='
  printlog,logfile,'Getting coordinates for ',file
  printlog,logfile,'=========================================='


  CD,filedir

  ; Test that the TOT and FITS files exists
  test = FILE_TEST(file)
  if test eq 1 then nlines = FILE_LINES(file)
  fitstest = FILE_TEST(fitsfile)
  if (test eq 0) or (nlines eq 0) or (fitstest eq 0) then begin
    PUSH,failurelist,longfile
    if test eq 0 then printlog,logfile,file,' NOT FOUND'
    if test eq 1 and nlines eq 0 then printlog,logfile,file,' HAS 0 LINES'
    if fitstest eq 0 then printlog,logfile,fitsfile,' NOT FOUND'
    goto,BOMB
  endif

  ; Loading the TOT file
  ;---------------------
  ; MAGFAP is the aperture magnitude for the final aperture for this star
  ; APCORR is the aperture correction for MAGFAP
  ; FINALAP is the number for the final aperture for this star
  ;            i.e. 5 for the 5th aperture
  ; MAG = MAGFAP + APCORR
  fieldnames = ['ID','X','Y','MAG','ERR','SKY','MAGFAP','APCORR','FINALAP']
  fieldtypes = [3,4,4,4,4,4,4,4,3]
  phot = IMPORTASCII(file,fieldnames=fieldnames,fieldtypes=fieldtypes,skip=3,/noprint)
  nphot = n_elements(phot)
  ;if size(phot,/type) ne 8 then begin
  ;  PUSH,failurelist,longfile
  ;  printlog,logfile,fitsfile,' DOES NOT HAVE A WCS'
  ;  goto,BOMB
  ;endif
  printlog,logfile,'Nstars = ',strtrim(nphot,2)

  ; Load the header
  if strmid(fitsfile,6,7,/reverse_offset) eq 'fits.fz' then begin
    head = HEADFITS(fitsfile,exten=1) else $
    ; Fix the NAXIS1/2 values in the header
    sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
    sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
  endif else begin
    head = HEADFITS(fitsfile)
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

  ; Add the RA/DEC tags
  ragd = where(tags eq 'RA',nragd)
  if (nragd) eq 0 then ADD_TAG,phot,'RA',0.0d0,phot
  decgd = where(tags eq 'DEC',ndecgd)
  if (ndecgd) eq 0 then ADD_TAG,phot,'DEC',0.0d0,phot

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

ENDFOR


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'STDRED_ASTROM Finished  ',systime(0)

if keyword_set(stp) then stop


end
