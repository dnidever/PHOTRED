;+
;
; PHOTRED_SPLIT
;
; This splits up multi-extension files (MEF) into individual
; chips/ampls using IRAF MSCRED's MSCSPLIT command.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The split MEF files.
;
; By D.Nidever  Feb 2008
;-

pro photred_split,redo=redo,stp=stp,std=std

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_SPLIT'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'SPLIT'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','photred_loadsetup','photred_getinput','photred_updatelists',$
         'maketemp','rndint','writeline','printlog','undefine','push','mktemp','strsplitter',$
         'touchzero','first_el','check_iraf']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; LOAD THE SETUP FILE if not passed
;-------------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count,std=std
  if count lt 1 then return
endif


; Get the IRAF directory from the setup file
;-------------------------------------------
irafdir = READPAR(setup,'IRAFDIR')
if irafdir eq '0' or irafdir eq '-1' then irafdir = '~/iraf/'
irafdir = FILE_SEARCH(irafdir,/fully_qualify,count=nirafdir)
if nirafdir eq 0 then begin
  print,'NO IRAF DIRECTORY'
  return
endif

; Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL
;---------------------------------------------------------------
CHECK_IRAF,iraftest,irafdir=irafdir
if iraftest eq 0 then begin
  print,'IRAF TEST FAILED.  EXITING'
  return
endif


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')

; Are we keeping the MEF files after splitting
dokeepmef = READPAR(setup,'KEEPMEF')
keepmef = 1
if dokeepmef eq '0' or dokeepmef eq '-1' then keepmef = 0




;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS or FITS.FZ files
; OUTLIST        FITS or FITS.FZ files
; SUCCESSLIST    FITS or FITS.FZ files


; Get input
;-----------
precursor = 'RENAME'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext=['fits','fits.fz'])
ninputlines = lists.ninputlines

; No files to process
;---------------------
if ninputlines lt 1 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines



;######################################
;#  PROCESSING THE FILES
;######################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''
printlog,logfile,systime(0)

; SPLIT MEF frames and pass through single chip files

; Initializing arrays
UNDEFINE,outlist,successlist,failurelist

; Get the image information
FOR i=0,ninputlines-1 do begin

  file = inputlines[i]
  filedir = FILE_DIRNAME(file)
  name = FILE_BASENAME(file)
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then base=FILE_BASENAME(file,'.fits.fz') else $
    base = FILE_BASENAME(file,'.fits')

  CD,filedir

  ; Does this file exist?
  test = FILE_TEST(file)
  if (test eq 0) then begin
    printlog,logfile,file,' NOT FOUND'
    PUSH,failurelist,file
    goto,BOMB
  endif

  ; Does this have multiple extensions
  FITS_OPEN,file,fcb,message=message0
  next = fcb.nextend
  FITS_CLOSE,fcb
  if next gt 0 then mef=1 else mef=0

  ;----------------------
  ; MULTI-EXTENSION FILE
  ;----------------------
  If (mef eq 1) then begin

    printlog,logfile,''
    printlog,logfile,file,' is a MULTI-EXTENSION file with ',strtrim(next,2),' extensions. SPLITTING'

    ; If the split files already exist delete them
    if (next ge 10) then begin
      strnext = strtrim(next,2)
      last = strmid(strnext,1,1)

      oldsplitfiles1 = FILE_SEARCH(base+'_[0-9].fits',/fully,count=noldsplitfiles1)
      if noldsplitfiles1 gt 0 then push,oldsplitfiles,oldsplitfiles1
      oldsplitfiles2 = FILE_SEARCH(base+'_1[0-'+last+'].fits',/fully,count=noldsplitfiles2)
      if noldsplitfiles2 gt 0 then push,oldsplitfiles,oldsplitfiles2
      noldsplitfiles = n_elements(oldsplitfiles)

    endif else begin
      oldsplitfiles = FILE_SEARCH(base+'_[0-'+strtrim(next,2)+'].fits',/fully,count=noldsplitfiles)
    endelse

    if noldsplitfiles gt 0 then begin
      printlog,logfile,'Deleting previous split files for ',name
      FILE_DELETE,oldsplitfiles,/allow,/quiet
    endif


    ; Split with MSCSPLIT
    ;--------------------
 
    ; Making IRAF script
    undefine,iraflines
    push,iraflines,'print("")' ; adorta: FIRST LINE WILL BE IGNORED!!
    push,iraflines,'cd '+filedir
    push,iraflines,'mscred'
    push,iraflines,'mscsplit(input="'+name+'",output="",mefext=".fits",delete=no,verbose=yes)'
    push,iraflines,'logout'
    ;cmdfile = MAKETEMP('temp','.cl')
    cmdfile = MKTEMP('temp')   ; absolute path
    WRITELINE,cmdfile,iraflines

    CD,irafdir

    ; Running IRAF
    undefine,out
    ;SPAWN,'cl < '+filedir+'/'+cmdfile,out,errout
    SPAWN,'cl < '+cmdfile,out,errout

    ; Return to original directory
    CD,filedir

    ; The output
    printlog,logfile,out

    ; Erasing the temporary file
    FILE_DELETE,cmdfile,/allow,/quiet


    ; How many were split?
    ;--------------------
    ; How many split files are there, don't include _0.fits
    if (next ge 10) then begin
      strnext = strtrim(next,2)
      last = strmid(strnext,1,1)

      undefine,splitfiles    
      splitfiles1 = FILE_SEARCH(base+'_[1-9].fits',/fully,count=nsplitfiles1)
      if nsplitfiles1 gt 0 then PUSH,splitfiles,splitfiles1
      ;splitfiles2 = FILE_SEARCH(base+'_1[0-'+last+'].fits',/fully,count=nsplitfiles2)
      maxten = long(next)/10
      splitfiles2 = FILE_SEARCH(base+'_[1-'+strtrim(maxten,2)+'][0-9].fits',/fully,count=nsplitfiles2)
      if nsplitfiles2 gt 0 then PUSH,splitfiles,splitfiles2
      nsplitfiles = n_elements(splitfiles)
    
    endif else begin
      splitfiles = FILE_SEARCH(base+'_[1-'+strtrim(next,2)+'].fits',/fully,count=nsplitfiles)
    endelse
    printlog,logfile,strtrim(nsplitfiles,2),' split files found for ',file


    ; SUCCESS!
    if (nsplitfiles eq next) then begin
      PUSH,outlist,splitfiles
      PUSH,successlist,file

    ; FAILURE!
    ; Do NOT add the _0.fits file to outputarr
    endif else begin
      PUSH,failurelist,file
    endelse


    ; Erase the MEF file
    ;-------------------
    if (keepmef eq 0) and (nsplitfiles eq next) then begin
      printlog,logfile,'Deleting the MEF file ',file
      FILE_DELETE,file,/allow,/quiet
    endif


  ;--------------------
  ; SINGLE CHIP FILE
  ;--------------------
  Endif else begin

    zero = STREGEX(file,'_0.fits',/boolean)

    ; Just pass it on
    if (zero eq 0) then begin
      PUSH,outlist,file
      PUSH,successlist,file
      printlog,logfile,file,' is a SINGLE chip file'
  
    ; FAILURE!
    ; Do NOT add the _0.fits files to outputarr
    endif else begin
      PUSH,failurelist,file
    endelse


  Endelse

  BOMB:

  CD,curdir

  ;#####################
  ; UPDATE the Lists
  ;#####################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent

  ;stop


ENDFOR


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'PHOTRED_SPLIT Finished  ',systime(0)

if keyword_set(stp) then stop

end
