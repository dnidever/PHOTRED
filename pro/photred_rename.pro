;+
;
; PHOTRED_RENAME
;
; This renames files so that it includes their field information.
; For example, ccd1001.fits gets renamed to F1.ccd1001.fits
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The object fits files are renamed with the field information
;  prepended.
;
; By D.Nidever  Feb 2008
;-

pro photred_rename,redo=redo,stp=stp,testing=testing,$
     continue_on_error=continue_on_error

COMMON photred,setup

print,''
print,'######################'
print,'RUNNING PHOTRED_RENAME'
print,'######################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'RENAME'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
inputfile = 'logs/'+thisprog+'.inlist'
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


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
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif

; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument, observatory
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')
observatory = READPAR(setup,'OBSERVATORY')
if observatory eq '0' or observatory eq '-1' or observatory eq '' then undefine,observatory
if strlowcase(telescope) eq 'blanco' then observatory='ctio'
if strlowcase(telescope) eq 'swope' then observatory='lco'
if strlowcase(telescope) eq 'magellan' then observatory='lco'
if strlowcase(telescope) eq 'lbt' then observatory='mgio'

; Separate field directories
sepfielddir = READPAR(setup,'SEPFIELDDIR')
if sepfielddir eq '0' or sepfielddir eq '-1' or sepfielddir eq '' then undefine,sepfielddir else sepfielddir=1


;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS files, it get files from directory
; OUTLIST        FITS files
; SUCCESSLIST    FITS files
undefine,outlist,successlist,failurelist

; Add all fits files to the INLIST
fitsfiles = FILE_SEARCH(['*.fits','*.fits.fz'],count=nfitsfiles,/fully)

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

; Remove files that start with "F[1-9]" or end in _0.fits
; *This is a final check on the entire INLIST*
; *and is not redundant with the code above*
bd1 = where(stregex(FILE_BASENAME(inputlines),'^F[1-9]',/boolean) eq 1,nbd1)
if nbd1 gt 0 then begin
  printlog,logfile,'Removed ',strtrim(nbd1,2),' files that already start with F#'
  PUSH,bad,bd1
endif
bd2 = where(stregex(FILE_BASENAME(inputlines),'_0.fits$',/boolean) eq 1,nbd2)
if nbd2 gt 0 then begin
  printlog,logfile,'Removed ',strtrim(nbd2,2),' files that end in _0.fits'
  PUSH,bad,bd2
endif
nbad = n_elements(bad)
if nbad gt 0 then begin
  ui = uniq(bad,sort(bad))   ; want unique ones
  ui = ui[sort(ui)]
  bad = bad[ui]
  nbad = n_elements(bad)

  ; Update the lists
  PUSH,failurelist,inputlines[bad]   ; add to failurelist
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent
endif

nleft = ninputlines - nbad
if (nleft gt 0) then begin
  if nbad gt 0 then REMOVE,bad,inputlines
  ninputlines = n_elements(inputlines)
endif else begin
  undefine,inputlines
  printlog,logfile,'NO FILES TO PROCESS'
  return
endelse



;##########################################################
;#  PROCESSING THE FILES
;##########################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''
printlog,systime(0)

; Initializing arrays
fieldarr = strarr(ninputlines)
calibarr = intarr(ninputlines)
stdarr = intarr(ninputlines)


; Check that the FITS header information can be properly interpreted
;-------------------------------------------------------------------
printlog,logfile,'CHECKING HEADER KEYWORDS'
headerproblem = 0
for i=0,ninputlines-1 do begin

  file = inputlines[i]
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
    base = FILE_BASENAME(file,'.fits.fz')
    head = HEADFITS(file,exten=1)
  endif else begin
    base = FILE_BASENAME(file,'.fits')
    head = HEADFITS(file)
  endelse
  com=''

  ; Checking GAIN
  gain = PHOTRED_GETGAIN(file)
  if gain le 0.0 then com=com+' GAIN ERROR,'

  ; Checking READNOISE
  rdnoise = PHOTRED_GETRDNOISE(file)
  if rdnoise le 0.0 then com=com+' READNOISE ERROR,'

  ; Checking UT-TIME
  uttime = PHOTRED_GETUTTIME(file)
  if (uttime eq '') then com=com+' UT-TIME ERROR,'

  ; Checking FILTER
  filter = SXPAR(head,'FILTER',count=nfilter,/silent)
  if (nfilter eq 0) then com=com+' FILTER ERROR,'

  ; Checking EXPTIME
  exptime = SXPAR(head,'EXPTIME',count=nexptime,/silent)
  if (nexptime eq 0) then com=com+' EXPTIME ERROR,'

  ; Checking RA
  ra = SXPAR(head,'RA',count=nra,/silent)
  if nra eq 0 then ra = SXPAR(head,'CRVAL1',count=nra,/silent)
  if (nra eq 0) then com=com+' RA ERROR,'

  ; Checking DEC
  dec = SXPAR(head,'DEC',count=ndec,/silent)
  if ndec eq 0 then dec = SXPAR(head,'CRVAL2',count=ndec,/silent)
  if (ndec eq 0) then com=com+' DEC ERROR,'

  ; Checking DATE
  date = PHOTRED_GETDATE(file)
  if (date eq '') then com=com+' DATE ERROR,'

  ; Checking AIRMASS
  airmass = PHOTRED_GETAIRMASS(file,obs=observatory,/update)
  if (airmass lt 0.9) then com=com+' AIRMASS ERROR'

  ; There were header problems
  if com ne '' then begin
    printlog,logfile,base,com
    if not keyword_set(continue_on_error) then PUSH,failurelist,file
    headerproblem = 1
    if not keyword_set(continue_on_error) then testing = 1
  endif

endfor

; UPDATE the Lists
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent


; There were HEADER problems
if (headerproblem eq 1) then begin
  printlog,logfile,''
  printlog,logfile,'HEADER KEYWORD problems.'
  if not keyword_set(continue_on_error) then retall
  ;printlog,logfile,'HEADER problems.  TESTING ONLY'
  ;printlog,logfile,''
endif else begin
  printlog,logfile,'HEADER KEYWORDS OKAY'
endelse


; What is the "Object" for each frame
;-------------------------------------
FOR i=0,ninputlines-1 do begin

  file = inputlines[i]
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
    base = FILE_BASENAME(file,'.fits.fz')
    head = HEADFITS(file,exten=1)
  endif else begin
    base = FILE_BASENAME(file,'.fits')
    head = HEADFITS(file)
  endelse

  object = SXPAR(head,'OBJECT',/silent)
  field = first_el(strsplit(object,' ',/extract))
  fieldarr[i] = field

  ; zero in name?
  zeroname = stregex(base,'zero',/fold_case,/boolean)
  ; bias in name?
  biasname = stregex(base,'bias',/fold_case,/boolean)
  ; flat in name?
  flatname = stregex(base,'flat',/fold_case,/boolean)

  ; zero in object string
  zeroobj = stregex(object,'zero',/fold_case,/boolean)
  ; bias in object string
  biasobj = stregex(object,'bias',/fold_case,/boolean)
  ; flat in object string
  flatobj = stregex(object,'flat',/fold_case,/boolean)
  ; twilight
  twiobj = stregex(object,'twil',/fold_case,/boolean)
  ; sky
  skyobj = stregex(object,'sky',/fold_case,/boolean)
  ; pointing
  pointobj = stregex(object,'pointing',/fold_case,/boolean)
  ; focus
  focusobj = stregex(object,'focus',/fold_case,/boolean)
  ; test
  testobj = stregex(object,'test',/fold_case,/boolean)


  ; standard star fields
  object2 = strcompress(object,/remove_all)  ; remove all whitespace
  sa98 = stregex(object2,'sa98',/fold_case,/boolean) OR $
         stregex(object2,'sa-98',/fold_case,/boolean)
  sa110 = stregex(object2,'sa110',/fold_case,/boolean) OR $
          stregex(object2,'sax-110',/fold_case,/boolean)
  sa114 = stregex(object2,'sa114',/fold_case,/boolean) OR $
          stregex(object2,'sa-114',/fold_case,/boolean)
  n3680 = stregex(object2,'n3680',/fold_case,/boolean) OR $
          stregex(object2,'ngc3680',/fold_case,/boolean)

  ; Check for standards against standard file
  stdfile=readpar(setup,'STDFILE')
  if (FILE_TEST(stdfile) eq 1) then begin
    ;readcol,stdfile,stdname,stdra,stdec,f='a,a,a',/silent
    readcol,stdfile,stdname,f='a',/silent
    stdcheck = max(stregex(stdname,object2,/fold_case,/boolean))
  endif else stdcheck=0
  ; This is a non-object frame
  if (zeroname eq 1) or (biasname eq 1) or (flatname eq 1) or $
     (zeroobj eq 1) or (biasobj eq 1) or (flatobj eq 1) or $
     (twiobj eq 1) or (skyobj eq 1) or (pointobj eq 1) or $
     (focusobj eq 1) or (testobj eq 1) then calibarr[i] = 1

  ; Standard star field
  if (sa98 eq 1) or (sa110 eq 1) or (sa114 eq 1) or $
     (n3680 eq 1) or (stdcheck eq 1) then stdarr[i] = 1

ENDFOR

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

  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
    base = FILE_BASENAME(file,'.fits.fz')
    head = HEADFITS(file,exten=1)
  endif else begin
    base = FILE_BASENAME(file,'.fits')
    head = HEADFITS(file)
  endelse

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

printlog,logfile,'PHOTRED_RENAME Finished  ',systime(0)

if keyword_set(stp) then stop

end
