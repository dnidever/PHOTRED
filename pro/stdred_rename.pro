;+
;
; STDRED_RENAME
;
; This renames files so that it includes their filter and night
; information.  For example, ccd1001.fits gets renamed to
; Mn1-ccd1001.fits.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /testing  Do NOT rename anything.  Just testing.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The fits filenames are prepended with filter and night
;  information.
;
; By D.Nidever  May 2008
;-

pro stdred_rename,redo=redo,stp=stp,testing=testing

COMMON photred,setup

print,''
print,'######################'
print,'RUNNING STDRED_RENAME'
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
printlog,logfile,'Starting STDRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','printlog','undefine','strsplitter',$
         'check_iraf','push','maketemp','rndint','writeline','photred_getinput',$
         'photred_updatelists','photred_getexptime','photred_getuttime','photred_getfilter','first_el',$
         'photred_loadsetup','photred_getairmass','photred_getdate','badpar','airmass','mktemp',$
         'photred_getgain','photred_getrdnoise','touchzero','sexig2ten']
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
  PHOTRED_LOADSETUP,setup,count=count,/std
  if count lt 1 then return
endif



; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument, observatory
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')
observatory = READPAR(setup,'OBSERVATORY')
if strlowcase(telescope) eq 'blanco' then observatory='ctio'
if strlowcase(telescope) eq 'swope' then observatory='lco'
if strlowcase(telescope) eq 'magellan' then observatory='lco'
if strlowcase(telescope) eq 'lbt' then observatory='mgio'
if observatory eq '0' or observatory eq '-1' or observatory eq '' then begin
  printlog,logfile,'OBSERVATORY not input'
  return
endif
; Get observatory structure
OBSERVATORY,observatory,obs_struct
if obs_struct.name eq '' then begin
  printlog,logfile,'Observatory >>',observatory,'<< NOT FOUND'
  return
endif


;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS or FITS.FZ files, it get files from directory
; OUTLIST        FITS or FITS.FZ files
; SUCCESSLIST    FITS or FITS.FZ files
undefine,outlist,successlist,failurelist

; Add all fits files to the INLIST
; This APPENDS files to the input file
fitsfiles = FILE_SEARCH(['*.fits','*.fits.fz'],count=nfitsfiles,/fully)
WRITELINE,inputfile,fitsfiles,/append

; Get input
;-----------
;lists = PHOTRED_GETINPUT(thisprog,redo=redo)
; Do NOT allow /redo or the FITS will get renamed TWICE
lists = PHOTRED_GETINPUT(thisprog)
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
bd = where(stregex(FILE_BASENAME(inputlines),'_0.fits$',/boolean) eq 1,nbd)
if nbd gt 0 then begin
  printlog,logfile,'Removed ',strtrim(nbd,2),' files that end in _0.fits'

  ; Update the lists
  PUSH,failurelist,inputlines[bd]   ; add to failurelist
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent
endif

nleft = ninputlines - nbd
if (nleft gt 0) then begin
  if nbd gt 0 then REMOVE,bd,inputlines
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
  if gain lt 0.0 then com=com+' GAIN ERROR,'

  ; Checking READNOISE
  rdnoise = PHOTRED_GETRDNOISE(file)
  if rdnoise lt 0.0 then com=com+' READNOISE ERROR,'

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
    PUSH,failurelist,file
    headerproblem = 1
    testing = 1
  endif

endfor

; UPDATE the Lists
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent

; There were HEADER problems
if (headerproblem eq 1) then begin
  printlog,logfile,''
  printlog,logfile,'HEADER KEYWORD problems.'
  retall
  ;printlog,logfile,'HEADER problems.  TESTING ONLY'
  ;printlog,logfile,''
endif else begin
  printlog,logfile,'HEADER KEYWORDS OKAY'
endelse



; Getting NIGHT and FILTER information
;-------------------------------------
jdarr = fltarr(ninputlines)
mjdarr = lonarr(ninputlines)
filtarr = strarr(ninputlines)
objarr = strarr(ninputlines)
calibarr = intarr(ninputlines)
FOR i=0,ninputlines-1 do begin

  file = inputlines[i]
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
    base = FILE_BASENAME(file,'.fits.fz')
    head = HEADFITS(file,exten=1)  ; load the header
  endif else begin
    base = FILE_BASENAME(file,'.fits')
    head = HEADFITS(file)  ; load the header
  endelse

  ; Getting object information
  object = SXPAR(head,'OBJECT',/silent)
  objarr[i] = object

  ; Getting filter information
  filt = PHOTRED_GETFILTER(file)
  filtarr[i] = filt

  ; Getting NIGHT information
  time = PHOTRED_GETUTTIME(file)
  timearr = strsplit(time,':',/extract)
  hour = long(timearr[0])
  min = long(timearr[1])
  sec = float(timearr[2])
  date = PHOTRED_GETDATE(file)
  datearr = strsplit(date,'-',/extract)
  year = long(datearr[0])
  month = long(datearr[1])
  day = long(datearr[2])

  jd = JULDAY(month,day,year,hour,min,sec)
  jdarr[i] = jd

  ; Get MJD night number
  mjdarr[i] = PHOTRED_GETMJD(file,observatory)

  ; Is this a calibration image (bias, flat, etc.)?
  ;------------------------------------------------
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

  ; This is a non-object frame
  if (zeroname eq 1) or (biasname eq 1) or (flatname eq 1) or $
     (zeroobj eq 1) or (biasobj eq 1) or (flatobj eq 1) or $
     (twiobj eq 1) or (skyobj eq 1) or (pointobj eq 1) or $
     (focusobj eq 1) or (testobj eq 1) then calibarr[i] = 1

ENDFOR


; Getting non-calibration frames
gd = where(calibarr eq 0,ngd)

; Some object frames
if (ngd gt 0) then begin

  ; Getting unique FILTERS
  ;-----------------------
  ui = uniq(filtarr[gd],sort(filtarr[gd]))
  filters = filtarr[gd[ui]]
  nfilters = n_elements(filters)

  ; Getting unique NIGHTS
  ;----------------------
  ;; Correct to LOCAL time
  ;timeoff = obs_struct.tz
  ;jdlocal = jdarr - timeoff/24.0
  ;minjd = floor(min(jdlocal))
  ;; JD starts at NOON, so 0.5 is midnight
  ;jdlocal2 = jdlocal-minjd
  ;night = floor(jdlocal2)+1
  night = mjdarr-min(mjdarr)+1

  ; unique nights
  ui = uniq(night[gd],sort(night[gd])) ; only care about object frames
  nights = night[gd[ui]]
  nnights = n_elements(nights)
  
  seqnight = intarr(ninputlines)-1      ; nights in a sequence, 1, 2, 3, etc...

  ; Make sequential night numbers
  for i=0,ngd-1 do begin
    inight = night[gd[i]]
    nightind = where(nights eq inight,nnightind)
    seqnight[gd[i]] = nightind[0]+1
  endfor

endif



; Initializing some arrays
newfilearr = strarr(ninputlines)


printlog,logfile,'-------------------------------------------------------------------------------------------------------------'
printlog,logfile,'     TYPE          FILENAME               NEW FILENAME            OBJECT         FILTER  EXPTIME     UT TIME'
printlog,logfile,'-------------------------------------------------------------------------------------------------------------'

; Rename the files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
    base = FILE_BASENAME(file,'.fits.fz')
    head = HEADFITS(longfile,exten=1)  ; load the header
  endif else begin
    base = FILE_BASENAME(file,'.fits')
    head = HEADFITS(longfile)          ; load the header
  endelse

  object = SXPAR(head,'OBJECT',/silent)
  exptime = PHOTRED_GETEXPTIME(longfile)
  filter = PHOTRED_GETFILTER(longfile)
  ut = PHOTRED_GETUTTIME(longfile)


  ; This is a calibration frame
  ;------------------------------
  if (calibarr[i] eq 1) then begin

    ; Make the calib directory if it doesn't exist
    if FILE_TEST('calib',/directory) eq 0 then FILE_MKDIR,'calib'

    newfile = 'calib/'+file
    newfilearr[i] = newfile

    com='Calib'
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

    ; Make prefix, filtershortname+'n'+nightnumber+'-'+fitsname
    prefix = filtarr[i]+'n'+strtrim(seqnight[i],2)
    sep = '-'
    newfile = prefix+sep+file
    newfilearr[i] = newfile


    com = 'Object'
    ;fmt = '(A10,A20,A25,A25,A5,A8,A16)'
    fmt = '(A10,A20,A25,A25,A5,F8.2,A16)'
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

printlog,logfile,'STDRED_RENAME Finished  ',systime(0)

if keyword_set(stp) then stop

end
