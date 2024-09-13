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
; HYPERTHREAD
hyperthread = READPAR(setup,'hyperthread')
if hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if strtrim(hyperthread,2) eq '0' then hyperthread=0
if n_elements(hyperthread) eq 0 then hyperthread=0
; NMULTI
if n_elements(nmulti) eq 0 then begin
  nmulti = READPAR(setup,'NMULTI')
  nmulti = long(nmulti)
  ; Use NMULTI_WCS if set
  nmultiwcs = READPAR(setup,'NMULTI_WCS')
  if nmultiwcs ne '0' and nmultiwcs ne '' and nmultiwcs ne '-1' then nmulti=long(nmultiwcs)
endif
nmulti = nmulti > 1  ; must be >=1
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
; SEPCHIPDIR
sepchipdir = READPAR(setup,'sepchipdir')
if sepchipdir eq '0' or sepchipdir eq '' or sepchipdir eq '-1' then sepchipdir=0

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
imagerfile = scriptsdir+'/imagers'
; The columns need to be: Telescope, Instrument, Naps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(imagerfile,fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#')
imagers.telescope = strupcase(strtrim(imagers.telescope,2))
imagers.instrument = strupcase(strtrim(imagers.instrument,2))
imagers.observatory = strupcase(strtrim(imagers.observatory,2))
singleind = where(imagers.namps eq 1,nsingle)
if nsingle gt 0 then imagers[singleind].separator = ''
if (n_tags(imagers) eq 0) then begin
  printlog,logfile,'NO imagers in '+imagerfile
  return
endif

filterfile = scriptsdir+'/filters'
if file_test(filterfile) eq 0 then begin
  print,'no >>filters<< file in ',scriptsdir
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
undefine,cmd,cmddir

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
  printlog,logfile,'Checking files for field '+ishortfield
  printlog,logfile,'-----------------------------'

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
    continue
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
    continue
  endif

  sfieldlines = "['"+strjoin(fieldlines,"','")+"']"
  cmd1 = "PHOTRED_COMBINE_SINGLE,"+sfieldlines
  cmd1 += ",'"+imagerfile+"','"+ifield+"'"
  cmd1 += ",telescope='"+telescope+"',instrument='"+instrument+"'"
  cmd1 += ",catformat='"+catformat+"'"
  if keyword_set(force) then cmd1 += ',/force'
  if keyword_set(posonly) then cmd1 += ',/posonly'
  if keyword_set(mchusetiles) then cmd1 += ',/mchusetiles'
  if keyword_set(sepchipdir) gt 0 then cmd1 += ',/sepchipdir'

  PUSH,cmd,cmd1
  PUSH,cmddir,curdir+'/'+ishortfield

  ;; make sure there's a "filters" file in each field directory
  if file_test(curdir+'/'+ishortfield+'/filters') eq 0 then $
    file_copy,filterfile,curdir+'/'+ishortfield,/over,/allow

ENDFOR
ncmd = n_elements(cmd)

if ncmd gt 0 then begin
  cmd = "cd,'"+cmddir+"' & "+cmd ; go to the directory
  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='cmb',hyperthread=hyperthread,/idle,$
             waittime=1,/cdtodir,scriptsdir=scriptsdir
endif

;; Check that the outfiles exist
for i=0,nsfields-1 do begin
  ishortfield = sfields[i]
  gd = where(shortfieldarr eq ishortfield,ngd)
  fieldlines = inputlines[gd]
  ;; Getting the main name for the output file
  firstname = FILE_BASENAME(fieldlines[0],'.phot')
  basedir = FILE_DIRNAME(fieldlines[0])
  if keyword_set(mchusetiles) then basedir=FILE_DIRNAME(basedir)  ;; F1/F1-T10/F1-00507801+T10.phot
  if not keyword_set(mchusetiles) and keyword_set(sepchipdir) then basedir=FILE_DIRNAME(basedir)  ;; F1/chip01/
  ;; Getting basename
  basename = firstname   ;; namps=1 and not using tiles
  if keyword_set(mchusetiles) then basename = (strsplit(firstname,'\'+tilesep+'T',/extract))[0]
  if not keyword_set(mchusetiles) and (thisimager.namps gt 1) then begin
    ending = first_el(strsplit(firstname,thisimager.separator,/extract),/last)
    endlen = strlen(ending)
    len = strlen(firstname)
    basename = strmid(firstname,0,len-endlen-1)
  endif
  outname = ishortfield+'/'+basename+'.cmb'
  ;; Make sure that it exists, and add to OUTLIST
  outtest = FILE_TEST(outname)
  if outtest eq 1 then outlines = FILE_LINES(outname) else outlines=0
  if (outlines gt 0) then begin
    PUSH,outlist,outname                       ; add final file to outlist
    PUSH,successlist,fieldlines                ; add all individual files to successlist
  endif else begin
    if outtest eq 0 then printlog,logfile,outname+' NOT FOUND'
    if outtest eq 1 and outlines eq 0 then printlog,logfile,outname+' HAS 0 LINES'
    PUSH,failurelist,fieldlines
  endelse
endfor

;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


printlog,logfile,'PHOTRED_COMBINE Finished  '+systime(0)

if keyword_set(stp) then stop

end
