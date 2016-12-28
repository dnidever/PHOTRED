;+
;
; PHOTRED_COMBINE
;
; This program combines the data from all the chips
; See comb_mscphot.pro
;
; INPUTS:
;  /redo   Redo files that were already done.
;  /force  Force combination even if not all amps are there.
;  /stp   Stop at the end of the program.
;
; OUTPUTS:
;  The combined photometry files for each pointing
;
; By D.Nidever  Mar 2008
;-

pro photred_combine,redo=redo,stp=stp,force=force,posonly=cmbposonly

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_COMBINE'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'COMBINE'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','add_tag','photmatch','phot_overlap','srcmatch',$
         'photred_getinput','photred_updatelists','photred_loadsetup','push','undefine','printlog',$
         'printstr','mktemp','first_el','odd','strsplitter','combine_structs','touchzero','writeline',$
         'minloc','range','stress','strep']
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


; LOAD information from the "photred.setup" file
;-----------------------------------------------
; REDO
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1
; TELESCOPE
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(strtrim(telescope,2))
if (telescope eq '0' or telescope eq '' or telescope eq '-1') then begin
  printlog,logfile,'NO TELESCOPE FOUND.  Please add to >>photred.setup<< file'
  return
endif
; INSTRUMENT
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(strtrim(instrument,2))
if (instrument eq '0' or instrument eq '' or instrument eq '-1') then begin
  printlog,logfile,'NO INSTRUMENT FOUND.  Please add to >>photred.setup<< file'
  return
endif
; FORCE combination
cmbforce = READPAR(setup,'CMBFORCE')
cmbforce = strupcase(strtrim(cmbforce,2))
if keyword_set(force) or (cmbforce ne '-1' and cmbforce ne '0') then force=1
; CMBPOSONLY, only use astrometry for matching
if n_elements(cmbposonly) eq 0 then begin
  cmbposonly = READPAR(setup,'CMBPOSONLY')
  if cmbposonly eq '0' or cmbposonly eq '' or cmbposonly eq '-1' then cmbposonly=0
  if cmbposonly eq '1' then cmbposonly=1 else cmbposonly=0
endif

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif



; LOAD THE "imagers" FILE
;----------------------------
printlog,logfile,'Loading imager information'
imagerstest = FILE_TEST(scriptsdir+'/imagers')
if (imagerstest eq 0) then begin
  printlog,logfile,'NO >>imagers<< file in '+scriptsdir+'  PLEASE CREATE ONE!'
  return
endif
; The columns need to be: Telescope, Instrument, Naps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(scriptsdir+'/imagers',fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#')
imagers.telescope = strupcase(strtrim(imagers.telescope,2))
imagers.instrument = strupcase(strtrim(imagers.instrument,2))
imagers.observatory = strupcase(strtrim(imagers.observatory,2))
singleind = where(imagers.namps eq 1,nsingle)
if nsingle gt 0 then imagers[singleind].separator = ''
if (n_tags(imagers) eq 0) then begin
  printlog,logfile,'NO imagers in '+scriptsdir+'/imagers'
  return
endif



; What IMAGER are we using??
;---------------------------
ind_imager = where(imagers.telescope eq strupcase(telescope) and imagers.instrument eq strupcase(instrument),nind_imager)
if nind_imager eq 0 then begin
  printlog,logfile,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
  return
endif
thisimager = imagers[ind_imager[0]]
; print out imager info
printlog,logfile,''
printlog,logfile,'USING IMAGER:'
printlog,logfile,'Telescope = '+thisimager.telescope
printlog,logfile,'Instrument = '+thisimager.instrument
printlog,logfile,'Namps = '+strtrim(thisimager.namps,2)
printlog,logfile,"Separator = '"+thisimager.separator+"'"
printlog,logfile,''



;###################
; GETTING INPUTLIST
;###################
; INLIST         PHOT files
; OUTLIST        PHOT files
; SUCCESSLIST    PHOT files

; Get input
;-----------
precursor = 'CALIB'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='phot')
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

; Load the "fields" file
;-------------------
testfields = FILE_TEST('fields')
if (testfields eq 0) then begin
  printlog,logfile,'>>fields<< FILE NOT FOUND'
  return
endif
READCOL,'fields',shortfields,fields,format='A,A',/silent
nfields = n_elements(fields)
if (nfields eq 0) then begin
  printlog,logfile,'NO FIELDS in >>fields<< FILE'
  return
endif

; Getting the unique short FIELD names from the filenames
;--------------------------------------------------------
shortfieldarr = strarr(ninputlines)
For i=0,ninputlines-1 do begin
  sep = '-'  ; '.'
  file = FILE_BASENAME(inputlines[i])
  shortfield = first_el(strsplit(file,sep,/extract))
  shortfieldarr[i] = shortfield
End

ui = uniq(shortfieldarr,sort(shortfieldarr))
ui = ui[sort(ui)]
sfields = shortfieldarr[ui]
nsfields = n_elements(sfields)
printlog,logfile,strtrim(nsfields,2)+' unique fields found'


;--------------------------------
; LOOP through the unique FIELDS
;--------------------------------
FOR i=0,nsfields-1 do begin

  ishortfield = sfields[i]

  printlog,logfile,''
  printlog,logfile,'==================================='
  printlog,logfile,'Combining photometry files for '+ishortfield
  printlog,logfile,'==================================='
  printlog,logfile,''
  printlog,logfile,systime(0)

  gd = where(shortfieldarr eq ishortfield,ngd)
  fieldlines = inputlines[gd]
  nfieldlines = n_elements(fieldlines)
  base = FILE_BASENAME(fieldlines,'.phot')
  basedir = FILE_DIRNAME(fieldlines[0])
  printlog,logfile,strtrim(nfieldlines,2)+' File(s) for FIELD='+ishortfield

  ; Getting the FIELD name
  gdfield = where(shortfields eq ishortfield,ngdfield)
  if (ngdfield gt 0) then begin
    ifield = fields[gdfield[0]]
  endif else begin
    printlog,logfile,'NO MATCH FOR '+ishortfield+' in the >>fields<< file'
    PUSH,failurelist,inputlines[gd]
    goto,BOMB1
  endelse



  ;###########################
  ; MULTI-AMP IMAGERS
  ;###########################
  If (thisimager.namps gt 1) then begin


    ; Check that we get ALL files for this GROUP
    ;-------------------------------------------
    ; Some previous successes to check
    nsuccess = lists.nsuccesslines
    if (nsuccess gt 0) then begin
      printlog,logfile,''
      printlog,logfile,'Some previous successes.  Making sure we have all files for this field'

      successbase = FILE_BASENAME(lists.successlines,'.phot')
      matchind = where(stregex(successbase,'^'+ishortfield+'-',/boolean),nmatchind)

      ; Found some matches
      if (nmatchind gt 0) then begin

        ; Check if these are already in the INLIST
        undefine,ind1,ind2,num_alreadyinlist
        MATCH,successbase[matchind],base,ind1,ind2,count=num_alreadyinlist
        num_notinlist = nmatchind - num_alreadyinlist

        ; Some not in INLIST yet
        if (num_notinlist gt 0) then begin      
      
          printlog,logfile,'Found '+strtrim(num_notinlist,2)+' previously successful file(s) for this group NOT YET in the '+$
                           'INLIST.  Adding.'
          indtoadd = matchind
          if num_alreadyinlist gt 0 then REMOVE,ind1,indtoadd
          PUSH,base,successbase[indtoadd]
          PUSH,fieldlines,lists.successlines[indtoadd]

          ; Setting REDO=1 so the files in the success list will be redone.
          if not keyword_set(redo) then begin
            printlog,logfile,'Setting REDO=1'
            redo = 1
          endif

        end  ; some not in inlist yet
      end  ; some files from this group in success file

      ; Make sure they are unique
      ui = uniq(fieldlines,sort(fieldlines))
      ui = ui[sort(ui)]
      fieldlines = fieldlines[ui]
      base = base[ui]
      nfieldlines = n_elements(fieldlines)

      printlog,logfile,''

    endif ; some successes

    ; Not enough files for this imager
    if (nfieldlines ne thisimager.namps) and not keyword_set(force) then begin
      printlog,logfile,'Only '+strtrim(nfieldlines,2)+' for '+ifield+'.  Need '+strtrim(thisimager.namps,2)+' for '+$
                       thisimager.telescope+'+'+thisimager.instrument
      PUSH,failurelist,fieldlines
      goto,BOMB1
    endif


    printlog,logfile,'COMBINING '+strtrim(nfieldlines,2)+' files for '+ifield

    ; Getting the main name
    firstname = FILE_BASENAME(fieldlines[0],'.phot')
    basedir = FILE_DIRNAME(fieldlines[0])

    ; Getting basename
    ending = first_el(strsplit(firstname,thisimager.separator,/extract),/last)
    endlen = strlen(ending)
    len = strlen(firstname)
    basename = strmid(firstname,0,len-endlen-1)
 
    ;-------------------------------------------------
    ; LOOP through all the PHOT files for this field
    ;-------------------------------------------------
    ncombined = 0
    undefine,str,all
    undefine,fieldnames0,fieldtypes0
    For j=0,nfieldlines-1 do begin

      file = fieldlines[j]
      filebase = FILE_BASENAME(file,'.phot')

      ; Check that the PHOT file exists
      test = FILE_TEST(file)
      if test eq 1 then photlines = FILE_LINES(file) else photlines=0
      if (photlines eq 0) then begin
        PUSH,failurelist,file
        if test eq 0 then printlog,logfile,file+' NOT FOUND'
        if test eq 1 and nlines eq 0 then printlog,logfile,file+' HAS 0 LINES'
        goto,BOMB2
      endif

      ; Get the fieldnames and fieldtypes
      ; We need ID to be STRING not LONG
      tempfile = MKTEMP('temp')
      SPAWN,'head '+file+' >> '+tempfile,out,errout
      temp = IMPORTASCII(tempfile,/header,/noprint)
      FILE_DELETE,tempfile,/allow,/quiet
      fieldnames = TAG_NAMES(temp)
      nfieldnames = n_elements(fieldnames)
      fieldtypes = lonarr(nfieldnames)
      for k=0,nfieldnames-1 do fieldtypes[k] = SIZE(temp[0].(k),/type)
      idind = where(fieldnames eq 'ID',nidind)
      fieldtypes[idind[0]] = 7
      ; FORCE it to be the same for ALL chip files, otherwise we sometimes
      ;  get "type mismatch" errors with double/floats
      if j eq 0 then begin
        fieldnames0 = fieldnames
        fieldtypes0 = fieldtypes
        file0 = file
      endif
        
      ; Check that the fieldnames are the same
      if n_elements(fieldnames0) ne n_elements(fieldnames) then begin
        PUSH,failurelist,file
        printlog,logfile,''
        printlog,logfile,file+' format does NOT agree with '+file0
        goto,BOMB2
      endif
      if total(strcmp(fieldnames0,fieldnames)) ne n_elements(fieldnames0) then begin
        PUSH,failurelist,file
        printlog,logfile,''
        printlog,logfile,file+' format does NOT agree with '+file0
        goto,BOMB2
      endif

      ; Load the PHOT file
      str = IMPORTASCII(file,fieldnames=fieldnames0,fieldtypes=fieldtypes0,skip=1,/noprint)
      nstr = n_elements(str)


      ; Print out the file info
      printlog,logfile,''
      printlog,logfile,'ADDING '+filebase+' Nstars='+strtrim(nstr,2)
      printlog,logfile,''


      ; Adding EXT tag for amp number
      tags = TAG_NAMES(str)
      extgd = where(tags eq 'EXT',nextgd)
      if nextgd eq 0 then ADD_TAG,str,'EXT',0,str

      ; Getting the amplifier number (extension)
      ;-----------------------------------------
      ext = first_el(strsplit(filebase,thisimager.separator,/extract),/last)
      str.ext = ext

      ; Updating the IDs
      ; FIELD_EXT.IDNUMBER, i.e. 190L182a_5.17366
      ;---------------------------------------------
      ;id2 = ifield+thisimager.separator+ext+'.'+strtrim(str.id,2)
      id2 = ifield+'_'+ext+'.'+strtrim(str.id,2)
      str.id = id2


      ; PROBABLY CUT THIS SECTION OUT!!!!

     ; ; Adding Field XB/YB coordinates, MOSAIC only
     ; ;--------------------------------------------
     ; if (instrument eq 'MOSAIC') then begin
     ;
     ;   ; Add the XB/YB fields
     ;   xbgd = where(tags eq 'XB',nxbgd)
     ;   if nxbgd eq 0 then ADD_TAG,str,'XB',0.0,str
     ;   ybgd = where(tags eq 'YB',nybgd)
     ;   if nybgd eq 0 then ADD_TAG,str,'YB',0.0,str
     ;   
     ;   ; 16 extensions, amplifiers
     ;   if (nfieldlines gt 10) then begin
     ;
     ;     ; X/Y are amplifier coordinates
     ;
     ;     ; Adding XCH/YCH fields
     ;     xchgd = where(tags eq 'XCH',nxchgd)
     ;     if nxchgd eq 0 then ADD_TAG,str,'XCH',0.0,str
     ;     ychgd = where(tags eq 'YCH',nychgd)
     ;     if nychgd eq 0 then ADD_TAG,str,'YCH',0.0,str
     ;
     ;     ; The original X/Y are amp coordinates
     ;     str.xch = str.x
     ;     str.ych = str.y
     ;     if (ext mod 2) eq 0 then str.x=str.x+1024L   ; correcting the X values
     ;
     ;     ; Correcting the X/Y coordinates
     ;     xoff = ((ext-1) mod 8)*1024.
     ;     yoff = ((ext-1)/8)*4096.
     ;
     ;     str.xb = str.x + xoff
     ;     str.yb = str.y + yoff
     ;
     ;   endif   ; 16 extensions
     ;
     ;   ; 8 extensions, chips
     ;   if (nfieldlines le 8) then begin
     ;
     ;     ; X/Y are chip coordinates
     ;
     ;     ; Correcting the X/Y coordinates
     ;     xoff = ((ext-1) mod 4)*2048.
     ;     yoff = ((ext-1)/4)*4096.
     ;
     ;     str.xb = str.x + xoff
     ;     str.yb = str.y + yoff
     ;   endif  ; 8 extensions
     ;
     ; endif  ; mosaic, adding XB/YB


      ; Check that the structure has RA/DEC
      ;------------------------------------
      gdra = where(tags eq 'RA',ngdra)
      gddec = where(tags eq 'DEC',ngddec)
      if (ngdra eq 0) or (ngddec eq 0) then begin
        PUSH,failurelist,file
        printlog,logfile,file+' DOES NOT HAVE RA/DEC FIELDS'
        goto,BOMB2
      endif


      ; COMBINING
      ;-----------
      ; Adding together with PHOT_OVERLAP.PRO to deal with overlaps
      ; at the edges of chips/amplifiers
      ; This uses RA/DEC so if these aren't in the structures it won't
      ; "combine" stars properly.
      if (n_elements(all) gt 0) then begin

        PHOT_OVERLAP,all,str,outstr,silent=silent,posonly=cmbposonly
        noutstr = n_elements(outstr)

        ; There was a problem
        if (noutstr eq 0) then begin
          printlog,logfile,'Problem in PHOT_OVERLAP for '+file
          PUSH,failurelist,file
          goto,BOMB2
        endif        

        all = outstr
        ncombined++

      ; First time, initialize the "all" structure
      endif else begin
        all = str
        ncombined++
      endelse

      BOMB2:

    Endfor  ; loop through the individual PHOT files

    ; At least one file failed
    ; Field and files are successful only if *ALL* succeeded
    ;--------------------------------------------------------
    if (ncombined ne nfieldlines) then begin
      printlog,logfile,ifield,' FAILED because ',strtrim(ngd-ncombined,2),' Files FAILED'
      PUSH,failurelist,fieldlines
      goto,BOMB1
    endif

    ; Are there any stars?
    ;----------------------
    nall = n_elements(all)
    if (nall eq 0) then begin
      printlog,logfile,'NO STARS for FIELD='+ifield
      PUSH,failurelist,fieldlines
      goto,BOMB1
    endif

    ; We were successful!
    ;-----------------------
    printlog,logfile,''
    printlog,logfile,'Combined '+strtrim(ncombined,2)+' files'
    printlog,logfile,'Nstars = '+strtrim(nall,2)+' for FIELD='+ifield

    ; Output the data
    ;----------------
    outname = basedir+'/'+basename+'.cmb'
    printlog,logfile,'OUTPUTTING data to '+outname

    PRINTSTR,all,outname

    ; Make sure that it exists, and add to OUTLIST
    outtest = FILE_TEST(outname)
    if outtest eq 1 then outlines = FILE_LINES(outname) else outlines=0
    if (outlines gt 0) then begin
      PUSH,outlist,outname                           ; add final file to outlist
      PUSH,successlist,fieldlines                ; add all individual files to successlist
    endif else begin
      if outtest eq 0 then printlog,logfile,outname+' NOT FOUND'
      if outtest eq 1 and outlines eq 0 then printlog,logfile,outname+' HAS 0 LINES'
      PUSH,failurelist,fieldlines
    endelse


  ;###########################
  ; SINGLE AMP IMAGERS
  ;###########################
  Endif else begin

    ; More than one file
    if nfieldlines gt 1 then begin
      printlog,logfile,strtrim(nfieldlines,2)+' files found.  Only ONE possible with this imager'
      PUSH,failurelist,fieldlines
      goto,BOMB1
    endif

    printlog,logfile,'Updating IDs for '+ifield

    ; Change the IDs
    ;---------------
    file = fieldlines[0]
    filebase = FILE_BASENAME(file,'.phot')

    ; Check that the PHOT file exists
    test = FILE_TEST(file)
    if test eq 1 then photlines = FILE_LINES(file) else photlines=0
    if (photlines eq 0) then begin
      PUSH,failurelist,file
      if test eq 0 then printlog,logfile,file+' NOT FOUND'
      if test eq 1 and photlines eq 0 then printlog,logfile,file+' HAS 0 LINES'
      goto,BOMB1
    endif

    ; Get the fieldnames and fieldtypes
    ; We need ID to be STRING not LONG
    tempfile = MKTEMP('temp')
    SPAWN,'head '+file+' >> '+tempfile,out,errout
    temp = IMPORTASCII(tempfile,/header,/noprint)
    FILE_DELETE,tempfile,/allow,/quiet
    fieldnames = TAG_NAMES(temp)
    nfieldnames = n_elements(fieldnames)
    fieldtypes = lonarr(nfieldnames)
    for k=0,nfieldnames-1 do fieldtypes[k] = SIZE(temp[0].(k),/type)
    idind = where(fieldnames eq 'ID',nidind)
    fieldtypes[idind[0]] = 7

    ; Load the PHOT file
    str = IMPORTASCII(file,fieldnames=fieldnames,fieldtypes=fieldtypes,skip=1,/noprint)
    nstr = n_elements(str)

    ; Print out the file info
    printlog,logfile,filebase,' Nstars='+strtrim(nstr,2)

    ; Updating the IDs
    ; FIELD.IDNUMBER, i.e. 190L182a.17366
    ;---------------------------------------------
    id2 = ifield+'.'+strtrim(str.id,2)
    str.id = id2

    ; Output the data
    ; We will update the original file.
    ;-----------------------------------
    ;printlog,logfile,'Moving original file to '+file+'.orig'
    ;FILE_MOVE,file,file+'.orig',/overwrite,/allow_same       ; save the original
    outname = basedir+'/'+base+'.cmb'
    printlog,logfile,'OUTPUTTING data to '+outname
    PRINTSTR,str,outname

    ; Make sure that it exists, and add to OUTLIST
    outtest = FILE_TEST(outname)
    if outtest eq 1 then outlines = FILE_LINES(outname) else outlines=0
    if (outlines gt 0) then begin
      PUSH,outlist,outname
      PUSH,successlist,file
    endif else begin
      if outtest eq 0 then printlog,logfile,outname+' NOT FOUND'
      if outtest eq 1 and outlines eq 0 then printlog,logfile,outname+' HAS 0 LINES'
      ;printlog,logfile,'Moving original file back to '+file
      ;FILE_MOVE,file+'.orig',file,/overwrite,/allow_same
      PUSH,failurelist,file
    endelse


  Endelse  ; single file

  BOMB1:



  ;##########################################
  ;#  UPDATING LIST FILES
  ;##########################################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent

ENDFOR


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'PHOTRED_COMBINE Finished  '+systime(0)

if keyword_set(stp) then stop


 end
