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
; Catalog format to use
catformat = READPAR(setup,'catformat')
if catformat eq '0' or catformat eq '' or catformat eq '-1' then catformat='ASCII'
if catformat ne 'ASCII' and catformat ne 'FITS' then catformat='ASCII'
; MCHUSETILES
mchusetiles = READPAR(setup,'MCHUSETILES')
if mchusetiles eq '0' or mchusetiles eq '' or mchusetiles eq '-1' then undefine,mchusetiles
tilesep = '+'

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
Endfor

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
      endif  ; some not in inlist yet
    endif  ; some files from this group in success file
    ; Make sure they are unique
    ui = uniq(fieldlines,sort(fieldlines))
    ui = ui[sort(ui)]
    fieldlines = fieldlines[ui]
    base = base[ui]
    nfieldlines = n_elements(fieldlines)
    printlog,logfile,''
  endif         ; some successes
  
  ;; Not enough files for this imager
  if (thisimager.namps gt 1) and (nfieldlines ne thisimager.namps) and not keyword_set(force) then begin
    printlog,logfile,'Only '+strtrim(nfieldlines,2)+' for '+ifield+'.  Need '+strtrim(thisimager.namps,2)+' for '+$
                     thisimager.telescope+'+'+thisimager.instrument
    PUSH,failurelist,fieldlines
    goto,BOMB1
  endif

  printlog,logfile,'COMBINING '+strtrim(nfieldlines,2)+' files for '+ifield

  ;; Getting the main name
  firstname = FILE_BASENAME(fieldlines[0],'.phot')
  basedir = FILE_DIRNAME(fieldlines[0])
  if keyword_set(mchusetiles) then basedir=FILE_DIRNAME(basedir)  ;; F1/F1-T10/F1-00507801+T10.phot
  ;; Getting basename
  basename = firstname   ;; namps=1 and not using tiles
  if keyword_set(mchusetiles) then basename = (strsplit(firstname,'\'+tilesep+'T',/extract))[0]
  if not keyword_set(mchusetiles) and (thisimager.namps gt 1) then begin
    ending = first_el(strsplit(firstname,thisimager.separator,/extract),/last)
    endlen = strlen(ending)
    len = strlen(firstname)
    basename = strmid(firstname,0,len-endlen-1)
  endif

  ;; Get file information and find unique exposures
  undefine,allfiles
  for j=0,nfieldlines-1 do begin
    dir1 = file_dirname(fieldlines[j])
    base1 = file_basename(fieldlines[j],'.phot')
    LOADMCH,dir1+'/'+base1+'.mch',indivfiles
    push,allfiles,dir1+'/'+file_basename(indivfiles,'.als')+'.fits'
  endfor
  bd = where(file_test(allfiles) eq 0,nbd)
  if nbd gt 0 then allfiles[bd]+='.fz'
  ;; Get exposure names
  ;; single-amp F1-12340044.als
  ;; multi-amp  F1-12340044_10.als
  arr1 = strsplitter(file_basename(allfiles),'-',/extract)
  allexpnum = reform(arr1[1,*])         ;; single-amp
  if thisimager.namps gt 1 then begin   ;; multi-amp
    arr2 = strsplitter(reform(arr1[1,*]),thisimager.separator,/extract)
    allexpnum = reform(arr2[0,*])
  endif
  uiexp = uniq(allexpnum,sort(allexpnum))
  uexpnum = allexpnum[uiexp]
  nexp = n_elements(uexpnum)
  ;; Get exposure information
  expstr = replicate({expnum:'',filter:'',dateobs:'',mjd:0.0d0,exptime:0.0,cenra:0.0d0,cendec:0.0d0},nexp)
  PHOTRED_GATHERFILEINFO,allfiles[uiexp],filestr
  struct_assign,filestr,expstr
  expstr.expnum = uexpnum
  for j=0,nexp-1 do expstr[j].mjd = date2jd(expstr[j].dateobs,/mjd)
  si = sort(expstr.mjd)  ; put in chronological order
  expstr = expstr[si]
  printlog,logfile,strtrim(nexp,2)+' exposures'
  
  ;-------------------------------------------------
  ; LOOP through all the PHOT files for this field
  ;-------------------------------------------------
  ncombined = 0
  undefine,str,all
  undefine,fieldnames0,fieldtypes0
  For j=0,nfieldlines-1 do begin
    file = fieldlines[j]
    filebase = FILE_BASENAME(file,'.phot')
    filedir = FILE_DIRNAME(file)

    ; Check that the PHOT file exists
    test = FILE_TEST(file)
    if test eq 1 then photlines = FILE_LINES(file) else photlines=0
    if (photlines eq 0) then begin
      PUSH,failurelist,file
      if test eq 0 then printlog,logfile,file+' NOT FOUND'
      if test eq 1 and nlines eq 0 then printlog,logfile,file+' HAS 0 LINES'
      goto,BOMB2
    endif

    ;; Load the PHOT file
    str = PHOTRED_READFILE(file,count=nstr)

    ;; Get the column names and types
    fieldnames = tag_names(str)
    fieldtypes = lonarr(n_elements(fieldnames))
    for k=0,n_elements(fieldnames)-1 do fieldtypes[k]=size(str[0].(k),/type)
    ;; ID must be a string
    idind = where(fieldnames eq 'ID',nidind)
    if fieldtypes[idind] ne 7 then begin
      fieldtypes[idind] = 7
      str0 = str[0] & struct_assign,{dum:''},str0
      schema = create_struct(fieldnames[0],fix(str0.(0),type=fieldtypes[0]))
      for k=1,n_elements(fieldnames)-1 do schema=create_struct(schema,fieldnames[k],fix(str0.(k),type=fieldtypes[k]))
      str_orig = str
      str = replicate(schema,nstr)
      struct_assign,str_orig,str
      str.id = strtrim(str.id,2)
      undefine,str_orig
    endif

    ;; Print out the file info
    printlog,logfile,''
    printlog,logfile,'ADDING '+filebase+' Nstars='+strtrim(nstr,2)
    printlog,logfile,''

    ;; Adding EXT tag for amp number
    tags = TAG_NAMES(str)
    if (thisimager.namps gt 1) or keyword_set(mchusetiles) then begin
      extgd = where(tags eq 'EXT',nextgd)
      if nextgd eq 0 then ADD_TAG,str,'EXT',0,str
    endif

    ;; Single-amp and NOT using tiles
    ;;-------------------------------
    If (not keyword_set(mchusetiles)) and (thisimager.namps eq 1) then begin
      ;; Updating the IDs
      ;; FIELD.IDNUMBER, i.e. 190L182a.17366
      ;;---------------------------------------------
      id2 = ifield+'.'+strtrim(str.id,2)
      str.id = id2
    Endif
       
    ;; Multi-amps and NOT using tiles
    ;;-------------------------------
    If (not keyword_set(mchusetiles)) and (thisimager.namps gt 1) then begin
      ; FORCE the format to be the same for ALL chip files, otherwise we sometimes
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
      ; Getting the amplifier number (extension)
      ;-----------------------------------------
      ext = photred_getchipnum(filebase,thisimager)
      str.ext = ext
      ; Updating the IDs
      ; FIELD_EXT.IDNUMBER, i.e. 190L182a_5.17366
      ;---------------------------------------------
      id2 = ifield+'_'+ext+'.'+strtrim(str.id,2)
      str.id = id2
    Endif
      
    ;; Using TILES, need to reformat photometry structure
    If keyword_set(mchusetiles) then begin
      ;; The chip/amp information is not in the name
      ;; each image could have a different chip
      ;; F1-00507801+T2, no chip extension
      ;; Use the TILE Number as the EXT
      ext = first_el(strsplit(filebase,tilesep+'T',/extract),/last)
      str.ext = ext
      ; Updating the IDs
      ; FIELD_TILE.IDNUMBER, i.e. 190L182a_5.17366
      ;---------------------------------------------
      id2 = ifield+'_'+ext+'.'+strtrim(str.id,2)
      str.id = id2
    Endif
      
    ;; Reformat the photometry structure
    ;; Put the photometry in exposure columns
    ;;---------------------------------------
    ;; Get information on the individual images
    LOADMCH,filedir+'/'+filebase+'.mch',alsfiles
    fitsfiles = filedir+'/'+file_basename(alsfiles,'.als')+'.fits'
    bdfits = where(file_test(fitsfiles) eq 0,nbdfits)
    if nbdfits gt 0 then fitsfiles[bdfits]+='.fz'
    PHOTRED_GATHERFILEINFO,fitsfiles,filestr
    ;; Get exposure names
    add_tag,filestr,'expnum','',filestr
    arr1 = strsplitter(file_basename(filestr.file),'-',/extract)
    allexpnum = reform(arr1[1,*])         ;; single-amp
    if thisimager.namps gt 1 then begin   ;; multi-amp
      arr2 = strsplitter(reform(arr1[1,*]),thisimager.separator,/extract)
      allexpnum = reform(arr2[0,*])
    endif
    filestr.expnum = allexpnum
    ;; Reformat the photometry structure
    str_orig = str & undefine,str
    PHOTRED_COMBINE_REFORMATPHOT,str_orig,filestr,expstr,str

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
      ;; There was a problem
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

  ;; At least one file failed
  ;; Field and files are successful only if *ALL* succeeded
  ;;--------------------------------------------------------
  if (ncombined ne nfieldlines) then begin
    printlog,logfile,ifield,' FAILED because ',strtrim(ngd-ncombined,2),' Files FAILED'
    PUSH,failurelist,fieldlines
    goto,BOMB1
  endif

  ;; Are there any stars?
  ;;----------------------
  nall = n_elements(all)
  if (nall eq 0) then begin
    printlog,logfile,'NO STARS for FIELD='+ifield
    PUSH,failurelist,fieldlines
    goto,BOMB1
  endif
  
  ;; We were successful!
  ;;-----------------------
  printlog,logfile,''
  printlog,logfile,'Combined '+strtrim(ncombined,2)+' files'
  printlog,logfile,'Nstars = '+strtrim(nall,2)+' for FIELD='+ifield

  ;; Output the data
  ;;----------------
  outname = basedir+'/'+basename+'.cmb'
  printlog,logfile,'OUTPUTTING data to '+outname
  if catformat eq 'ASCII' then begin
    PRINTSTR,all,outname
    PRINTSTR,expstr,outname+'.meta'
  endif else begin
    MWRFITS,all,outname,/create
    MWRFITS,expstr,outname,/silent   ;; add meta-data on exposure information 
  endelse

  ;; Make sure that it exists, and add to OUTLIST
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

  BOMB1:

  ;##########################################
  ;#  UPDATING LIST FILES
  ;##########################################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,setupdir=curdir,/silent
ENDFOR


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


printlog,logfile,'PHOTRED_COMBINE Finished  '+systime(0)

if keyword_set(stp) then stop

end
