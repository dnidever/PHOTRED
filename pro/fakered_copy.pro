pro fakered_copy,redo=redo,stp=stp

;+
;
; FAKERED_COPY
;
; This copies all of the field files to initialize a mock
; for the artificial star tests.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The field files are copied to new mock directories.
;
; By D.Nidever  May 2016
;-

COMMON photred,setup

print,''
print,'######################'
print,'RUNNING FAKERED_COPY'
print,'######################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'COPY'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
inputfile = 'logs/'+thisprog+'.inlist'
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting FAKERED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','printlog','undefine','strsplitter','mktemp',$
         'check_iraf','push','maketemp','rndint','writeline','photred_getinput','photred_getgain',$
         'photred_updatelists','photred_getexptime','photred_getuttime','photred_getfilter','first_el',$
         'photred_loadsetup','photred_getairmass','photred_getdate','badpar','airmass',$
         'photred_getrdnoise','touchzero','sexig2ten']
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
  PHOTRED_LOADSETUP,setup,count=count,/fake
  if count lt 1 then return
endif

; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Separate field directories
spefielddir = READPAR(setup,'SEPFIELDDIR')
if sepfielddir eq '0' or sepfielddir eq '-1' or sepfielddir eq '' then undefine,sepfielddir else sepfielddir=1

; Number of separate mocks
nmocks = READPAR(setup,'NMOCKS')
if nmocks eq '0' or nmocks eq '-1' or nmocks eq '' then nmocks=1 else nmocks=long(nmocks)>1

; Number of artificial stars per mock
nmockstars = READPAR(setup,'NMOCKSTARS')
if nmockstars eq '0' or nmockstars eq '-1' or nmockstars eq '' then nmockstars=50000L else nmockstars=long(nmockstars)



; COPY
; fits, opt, als.opt, psf, mch
;  weights, scale, zero, comb_psf, _shift.mch
;  chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_shift.mch']




;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS files, it get files from directory
; OUTLIST        FITS files
; SUCCESSLIST    FITS files
undefine,outlist,successlist,failurelist

; Add all fits files to the INLIST
fitsfiles = FILE_SEARCH('*.fits',count=nfitsfiles,/fully)

; Remove files that start with "F[1-9]"
; since they were probably renamed already
; This APPENDS files to the input file
; Also don't include "check.fits"
if (nfitsfiles gt 0) then begin
  gd = where(stregex(FILE_BASENAME(fitsfiles),'^F[1-9]',/boolean) eq 0 and $
             (FILE_BASENAME(fitsfiles) ne 'check.fits'),ngd)
  if ngd gt 0 then $
  WRITELINE,inputfile,fitsfiles[gd],/append
endif

; Get input
;-----------
lists = PHOTRED_GETINPUT(thisprog,redo=redo)
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines



;##########################################################
;#  PROCESSING THE FILES
;##########################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''


; Making object fieldname list
;-------------------------------

; Object frames
gd = where(calibarr eq 0 and stdarr eq 0,ngd)

; What is the name of the "fields" file?
fieldsfile = 'fields'
if keyword_set(testing) then fieldsfile='fields.testing'

; Some objects found
if (ngd gt 0) then begin

  fieldnames = fieldarr[gd]

  ; Get unique field names
  ; CASE-INSENSITIVE
  Ufieldnames = strupcase(fieldnames)
  ui = uniq(Ufieldnames,sort(Ufieldnames))
  ui = ui[sort(ui)]
  fields = fieldnames[ui]
  nfields = n_elements(fields)

  printlog,logfile,strtrim(nfields,2),' fields found'
  printlog,logfile,fields

  ; Load old "fields" file
  undefine,oldshortnames,oldfields
  if FILE_TEST('fields') eq 1 then begin
    READCOL,'fields',oldshortnames,oldfields,format='A,A',/silent
  endif

  ; Combine them
  undefine,newfields
  push,newfields,oldfields
  push,newfields,fields
  Unewfields = strupcase(newfields)   ; case-insensitive
  ui = uniq(Unewfields,sort(Unewfields))
  ui = ui[sort(ui)]
  newfields = newfields[ui]
  nnewfields = n_elements(newfields)

  ; Make new "fields" file
  printlog,logfile
  printlog,logfile,'Printing fieldnames to "'+fieldsfile+'" file'
  shortnames = 'F'+strtrim(lindgen(nnewfields)+1,2)
  out = shortnames+'   '+newfields
  WRITELINE,fieldsfile,out

  fields = newfields
  nfields = n_elements(fields)

  ; Make separate field directories
  if keyword_set(sepfielddir) then $
    for i=0,nfields-1 do FILE_MKDIR,shortnames[i]

; No object frames
endif else begin

  printlog,logfile,'No OBJECT frames found'
  undefine,fields
  nfields = 0

endelse

; Initialzing some arrays
newfilearr = strarr(ninputlines)


printlog,logfile,'-------------------------------------------------------------------------------------------------------------'
printlog,logfile,'     TYPE          FILENAME               NEW FILENAME            OBJECT         FILTER  EXPTIME   UT TIME'
printlog,logfile,'-------------------------------------------------------------------------------------------------------------'

; Rename the files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  base = FILE_BASENAME(file,'.fits')

  ; Load the header
  head = HEADFITS(longfile)

  object = SXPAR(head,'OBJECT',/silent)
  exptime = PHOTRED_GETEXPTIME(longfile)
  filter = PHOTRED_GETFILTER(longfile)
  ut = PHOTRED_GETUTTIME(longfile)


  ; This is a calibration frame
  ;------------------------------
  if (calibarr[i] eq 1 or stdarr[i] eq 1) then begin

    ; Make the calib directory if it doesn't exist
    if FILE_TEST('calib',/directory) eq 0 then FILE_MKDIR,'calib'

    newfile = 'calib/'+file
    newfilearr[i] = newfile

    if calibarr[i] eq 1 then com='Calib'
    if stdarr[i] eq 1 then com='Standard'
    fmt = '(A10,A20,A25,A25,A5,A8,A16)'
    printlog,logfile,com,file,newfile,object,filter,exptime,ut,format=fmt

    ; Moving files
    if not keyword_set(testing) then begin
      FILE_MOVE,file,newfile,/over

      ; Did we have success
      test = FILE_TEST(newfile)
      ;SUCCESS!
      if test eq 1 then begin
        PUSH,successlist,longfile

      ; FAILURE!
      endif else begin
        PUSH,failurelist,longfile
      endelse
      ;if test eq 1 then successarr[i] = 1
    endif


  ; Object frame
  ;----------------
  endif else begin

    ; Get field name and "short" name, i.e. F3.ccd1001.fits
    ifield = first_el(strsplit(object,' ',/extract))
    ind = where(strupcase(fields) eq strupcase(ifield),nind)  ; case-insensitive
    shfield = shortnames[ind[0]]

    sep = '-'   ; '.'
    newfile = shfield+sep+file
    if keyword_set(sepfielddir) then newfile=shfield+'/'+newfile  ; separate field directory
    newfilearr[i] = newfile


    com = 'Object'
    fmt = '(A10,A20,A25,A25,A5,F8.2,A17)'
    printlog,logfile,com,file,newfile,object,filter,exptime,ut,format=fmt


    ; Moving
    if not keyword_set(testing) then begin
      FILE_MOVE,file,newfile,/over

      ; Did we have success
      test = FILE_TEST(newfile)
      ; SUCCESS!
      if test eq 1 then begin
        newfile2 = FILE_SEARCH(newfile,/fully)        ; fully qualified name
        PUSH,outlist,newfile2
        PUSH,successlist,longfile

      ; FAILURE!
      endif else begin
        PUSH,failurelist,longfile
      endelse
    endif  ; not testing

  endelse


  ;#####################
  ; UPDATE the Lists
  ;#####################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent
  ;stop

ENDFOR
printlog,logfile,'-------------------------------------------------------------------------------------------------------------'


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile
printlog,logfile,'PLEASE CHECK THAT ALL OF THE FILES WERE NAMED PROPERLY'
printlog,logfile,'IF NOT, RENAME THEM BY HAND AND UPDATE THE OUTLIST APPROPRIATELY'
printlog,logfile

; There were header problem. RETALL
if (headerproblem eq 1) then begin
  printlog,logfile,''
  printlog,logfile,'There were HEADER problems. RETURNING.'
  printlog,logfile,''
  RETALL
endif

printlog,logfile,'FAKERED_COPY Finished  ',systime(0)

if keyword_set(stp) then stop

end
