;+
;
; FAKERED_ADDSTAR
;
; This adds artificial stars to the mock images.
;
; INPUTS:
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The artificial stars are added to the mock images.
;
; By D.Nidever  March 2017
;   Antonio Dorta,  major update   June/July 2017
;-

function getField, fn, remove_path=remove_path
  ; Given a file, returns the FIELD. Filename is expected with
  ; format FX-NNNNNNN_YY.ext (FX will be returned)
  if keyword_set(remove_path) then fn = file_basename(fn)
  tmp = strsplit(fn, '-', /EXTRACT)   ; Field located in the beginning before "-"
  return, tmp[0]
end

function getChip, fn, remove_path=remove_path
  ; Given a file, returns the CHIP. Filename is expected with
  ; format FX-NNNNNNN_YY.ext (YY will be returned)
  if keyword_set(remove_path) then fn = file_basename(fn)
  tmp = strsplit(fn, '_', /EXTRACT)   ; Chip located between "_" and extension
  if n_elements(tmp) ne 2 then return,''
  tmp = strsplit(tmp[1], '.', /EXTRACT)   ; Remove extension
  return, tmp[0]
end


pro fakered_addstar,redo=redo,stp=stp,waittime=waittime

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING FAKERED_ADDSTAR'
print,'########################'
print,''

if not keyword_set(waittime) then waittime = 15

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'ADDSTAR'
logfile = 'logs'+PATH_SEP()+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting FAKERED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readpar', 'runshellcmd', 'printlog', 'undefine', 'check_python', 'check_pyhtcondor', 'push', $
         'photred_updatelists', 'first_el', 'photred_loadsetup']
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

;                          VARIABLE          NAME IN SETUP       SETUP  DEFAULT VALUE   LOGFILE       
hyperthread     = getparam(hyperthread     , 'hyperthread'     , setup, '0'           , logfile, /bool)
nmulti          = getparam(nmulti          , 'nmulti'          , setup, '1'           , logfile)
scriptsdir      = getparam(scriptsdir      , 'scriptsdir'      , setup, ''            , logfile)
sepfielddir     = getparam(sepfielddir     , 'sepfielddir'     , setup, '0'           , logfile, /bool)
sepchipdir      = getparam(sepchipdir      , 'sepchipdir'      , setup, '0'           , logfile, /bool)
starssplit      = getparam(starssplit      , 'starssplit'      , setup, '0'           , logfile, /bool)
starsfile       = getparam(starsfile       , 'starsfile'       , setup, ''            , logfile)
starscols       = getparam(starscols       , 'starscols'       , setup, ''            , logfile)
datadir         = getparam(datadir         , 'datadir'         , setup, '.'           , logfile)
datatransfer    = getparam(datatransfer    , 'datatransfer'    , setup, 'skip'        , logfile)
chipsfile       = getparam(chipsfile       , 'chipssfile'      , setup, '*chips.fits' , logfile)
maxccdsize      = getparam(maxccdsize      , 'maxccdsize'      , setup, '2048,4096'   , logfile)
magext          = getparam(magext          , 'magext'          , setup, ''            , logfile)
radcent         = getparam(radcent         , 'radcent'         , setup, '*'           , logfile)
dimfield        = getparam(dimfield        , 'dimfield'        , setup, '*'           , logfile)
distance        = getparam(distance        , 'distance'        , setup, '0'           , logfile)
pythonbin       = getparam(pythonbin       , 'pythonbin'       , setup, 'python'      , logfile)
htcondor        = getparam(htcondor        , 'htcondor'        , setup, '0'           , logfile)
htcondor_cmd    = getparam(htcondor_cmd    , 'htcondor_cmd'    , setup, ''            , logfile)
starsshuffle    = getparam(starsshuffle    , 'starsshuffle'    , setup, '0'           , logfile, /bool)
maxmocks        = getparam(maxmocks        , 'maxmocks'        , setup, '0'           , logfile)
refinemch       = getparam(refinemch       , 'refinemch'       , setup, '0'           , logfile)

; SOME EXTRA VALIDATIONS

; nmulti should have a numeric type (NOT string)
nmulti = long(nmulti)

; Check if scriptsdir exists
if scriptsdir eq '' or file_test(scriptsdir,/directory) ne 1 then begin
  printlog,logfile,'ERROR: SCRIPTSDIR "',scriptsdir,' NOT FOUND'
  error='SCRIPTDIRS "'+scriptsdir+' NOT FOUND'
  return
endif

; Check if starsfile exists
if starsfile eq '' or file_test(starsfile,/read) ne 1 then begin
  printlog,logfile,'ERROR: STARSFILE "',starsfile,' NOT FOUND'
  error='STARSFILE "'+starsfile+' NOT FOUND'
  return
endif

; Check if starscols is NOT empty
if starscols eq '' then begin
  printlog,logfile,'ERROR: STARSCOLS NOT DEFINED'
  error='STARSCOLS "'+starsfile+' NOT DEFINED'
  return
endif

; Check whether datatransfer is one of the valid methods
datatransfer_valid = ['skip', 'copy', 'move', 'link', 'inlist']
if where(datatransfer_valid eq datatransfer) eq -1 then begin 
  printlog,logfile,''
  printlog,logfile,'WARNING!! "' + datatransfer + '" is not a valid transfer option'
  printlog,logfile,'Please, check your setup file. Skipping file transferring...'
  printlog,logfile,''
  datatransfer = 'skip'
endif

; Check if the scripts exist in the current directory
scripts = ['fakered_transfer.sh']
nscripts = n_elements(scripts)
; Loop through the scripts
for i=0,nscripts-1 do begin
  fullpath = scriptsdir+PATH_SEP()+scripts[i]

  if file_test(fullpath,/read) ne 1 then begin
    printlog,logfile,fullpath,' NOT FOUND or EMPTY'
    error = fullpath+' NOT FOUND or EMPTY'
    return
  endif
endfor ; scripts loop


;---------------------------------------------------------------
; Run CHECK_PYTHON.PRO to make sure that you can run PYTHON from IDL
;---------------------------------------------------------------
CHECK_PYTHON,pythontest,pythonbin=pythonbin
if pythontest eq 0 then begin
  print,'PYTHON TEST FAILED.  EXITING'
  return
endif

; Check that the DAOPHOT programs exist
SPAWN,'which daophot',out,errout
daophotfile = FILE_SEARCH(out,count=ndaophotfile)
if (ndaophotfile eq 0) then begin
  print,'DAOPHOT PROGRAM NOT AVAILABLE'
  return
endif

if refinemch gt 0 then begin
  ; Check that the DAOMSATER programs exist
  SPAWN,'which daomaster',out,errout
  daomasterfile = FILE_SEARCH(out,count=ndaomasterfile)
  if (ndaomasterfile eq 0) then begin
    print,'DAOMASTER PROGRAM NOT AVAILABLE'
    return
  endif
endif



;##########################################################
;#  TRANSFERRING INPUT FILES (if needed)
;##########################################################
printlog,logfile,''
printlog,logfile,systime(0)

; Check if the user wants to SKIP the transfer (files are already there!)
if datatransfer ne 'skip' then begin

  printlog,logfile,''
  printlog,logfile,'--------------------------'
  printlog,logfile,' TRANSFERRING INPUT FILES '
  printlog,logfile,'--------------------------'
  printlog,logfile,''

  cmd = scriptsdir + PATH_SEP()+ 'fakered_transfer.sh'
  params = " '" + datatransfer + "' '" + datadir + "' '" + strcompress(sepfielddir) + "' '" + strcompress(sepchipdir) $
               + "' '" + starsfile + "' '" + strcompress(starssplit)  + "' '" + chipsfile + "'"
  msg = 'Finishing this stage since input files were not properly transferred'
  if runshellcmd(cmd+params,msg=msg,/printcmd,/printout) ne 0 then return
endif



;----------------------------
; MANAGE LOGS
;----------------------------
; INLIST         MCH files with format F*-*_*.mch
; OUTLIST        FITS files that were created
; SUCCESSLIST    FITS files
undefine,inlist,outlist,successlist,failurelist
 
; Add all mch files with format FX-NNNNN-YY.mch to the INLIST
;cmd='find `pwd` -regextype sed -regex ".*F[0-9]*-[0-9]*_[0-9]*\.mch"' 
;ec=runshellcmd(cmd,output=inlist,/quiet)

READLIST,'logs'+PATH_SEP()+thisprog+'.inlist',mchinputlines,/exist,/unique,/fully,count=nmchinputlines,logfile=logfile
if n_elements(mchinputlines) eq 0 then mchinputlines = []

;--------------------------------


;##########################################################
;#  PROCESSING THE FILES
;##########################################################
printlog,logfile,''
printlog,logfile,'------------------------'
printlog,logfile,' PROCESSING INPUT FILES '
printlog,logfile,'------------------------'
printlog,logfile,''


CD,current= curdir
undefine,cmd,cmddir
reldir      = ""
base_chipsfile  = curdir+PATH_SEP()+string(chipsfile)
base_script     = scriptsdir + PATH_SEP()+"fakered_addstar.py"
; Check if we have extra directories (for fields and/or chips)
if sepfielddir eq 1 then reldir += ".." + PATH_SEP()
if sepchipdir  eq 1 then reldir += ".." + PATH_SEP()

printlog,logfile,"Processing " + strcompress(nmchinputlines) + " files from INPUT list"
for i=0,nmchinputlines-1 do begin
  base_dir = file_dirname(mchinputlines[i])
  ; Make commands for addstar python program
  ; Syntax: python python_script.py  <1:mch_file> <2:chips_file> <3:maxmocks> <4:starscols> <5:starsshufle> 
  ;                                  <6:magext>   <7:maxccdsize> <8:dimfield> <9:radcent> <10:distance> <11:refinemch>
  cmd1  = pythonbin                 + " "   + base_script            + " '"    $
        + mchinputlines[i]          + "' '" + base_chipsfile         + "' '"   $
        + strcompress(maxmocks)     + "' '" + strcompress(starscols) + "' '"   $
        + strcompress(starsshuffle) + "' '" + strcompress(magext)    + "' '"   $
        + strcompress(maxccdsize)   + "' '" + strcompress(radcent)   + "' '"   $
        + strcompress(dimfield)     + "' '" + strcompress(distance)  + "' '"   $
        + strcompress(refinemch)    + "'"
    ; Add commands and directories
  PUSH,cmd,cmd1
  PUSH,cmddir,base_dir
endfor


if htcondor ne 0 then begin
  ; HTCondor is ENABLED!!

  ; Build the HTCondor submit commands (Use |=| for assignations and |;| to separate commands
  ; Add user's HTCondor commands and needed environment variables (PATH, LD_LIBRARY_PATH and HOME)
  htcondor_cmd = "#shared|;|" + htcondor_cmd + "|;|environment|=|PATH="+getenv('PATH')+";LD_LIBRARY_PATH="+ $
                  getenv('LD_LIBRARY_PATH')+";HOME="+getenv('HOME')

  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='py',hyperthread=hyperthread,waittime=waittime,    $
           htcondor=htcondor_cmd, scriptsdir=scriptsdir, pythonbin=pythonbin

endif else  $
  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='py',hyperthread=hyperthread,waittime=waittime,    $
           htcondor='', scriptsdir=scriptsdir, pythonbin=pythonbin,/cdtodir



;---------------------------------------
; UPDATE the Lists
;---------------------------------------

successlist=[]
failurelist=[]
; We will process the inlist to see which input files have produced outputs (.fits files)
; and then move it to SUCCESS or failed and then move to FAILURE
for i=0,nmchinputlines-1 do begin
  ; For each file in inlist, separate dirname(path) and basename(filename:fn)
  path    = file_dirname(mchinputlines[i])
  fn_mch  = file_basename(mchinputlines[i])
  ; ------------------------
  ; From basename get ONLY FIELD and CHIP (bear in mind that if sepfielddir and sepchipdir are false,
  ; all files will be together, so we need to get field and chip from each input file
  field = getField(fn_mch)
  chip  = getChip(fn_mch)
  ; ------------------------
  ; Get ALL .add files for all Mocks of that field and chip in that path
  fn_adds = path + PATH_SEP() + field + 'M*-*_'+ chip + '.add'
  check_adds = file_search(fn_adds, /FULLY)
  if n_elements(check_adds) lt 1 then begin
    push, failurelist, mchinputlines[i]          ; If there are no .add files, it is a failure
    continue
  endif

  ; Now check that for EACH .add should be a .fits file (if not, it is a failure)
  check_ok = 1
  for j=0,n_elements(check_adds)-1 do begin
    fn_fits = file_basename(check_adds[j], "add") + "fits"   ; change .add to .fits
    check_fits = path + PATH_SEP() + fn_fits
    info = file_info(check_fits)                             ; get info about file
    if info.exists eq 0 or info.size eq 0 then begin
      push, failurelist, mchinputlines[i]                    ; a fits file is missing or it has size=0, FAILURE!
      check_ok = 0
      break
    endif
  endfor
  if check_ok eq 1 then push, successlist, mchinputlines[i]  ; ALL fits files are there, SUCCESS!
endfor

; OUTLIST is all new .mch files for all mocks
outlist = file_search("", "F*M*-*_??.mch", /FULLY)

; CREATE LISTS
lists = photred_getinput(thisprog)
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,  $
                    failurelist=failurelist,/silent


; ---------------------------------------
;  GATHER AND MERGE PARTIAL FILES
; ---------------------------------------

; Some information like Aperture Correction (apcor) and Transformation Equations 
; are created partially for each chip and stored in partial files. 
; Now we will gather all partial info and merge it in one only file that will
; be used in next stage(s)
partial_files = [{fn_full:'apcor.lst',   fn_part:'*partial_apcor*.lst'}, $
                 {fn_full:'calib.trans', fn_part:'*partial_calib*.trans'}]

for i=0,n_elements(partial_files)-1 do begin
  ; If the destination file exists, save a BAK
  if file_test(partial_files[i].fn_full) then $    
    file_move, partial_files[i].fn_full, partial_files[i].fn_full+'.bak', /OVERWRITE

  ; Search for all partial files, read them and write all content in a single file
  partial_fn = file_search("", partial_files[i].fn_part, /FULLY)
  for j=0,n_elements(partial_fn)-1 do begin
    readline,  partial_fn[j], partial_data
    writeline, partial_files[i].fn_full, partial_data, /APPEND
  endfor
endfor

; ---------------------------------------
;  BUILD THE NEW fields FILE
; ---------------------------------------
; We will add the new "fake" fields (fieldmock: FXM1, FXM2, ... FXMN) 
; to the fields file, since it is required by later stages (like COMBINE)
; The long field name of the new fieldmocks will be the same as the old field
; Example (with 2 mocks per field):
;    OLD fields file                NEW field file
;      F1    Field10                  F1    Field10
;      F2    Field43                  F1M1  Field10
;                                     F1M2  Field10
;                                     F2    Field43
;                                     F2M1  Field43
;                                     F2M2  Field43
FIELDS_FILE='fields'
if file_test(FIELDS_FILE) then begin      
  newshortfields = []
  newfields = []
  ; Read old field file
  READCOL,FIELDS_FILE,oldshortfields,oldfields,format='A,A',/silent
  ; Backup old file (-> fields.bak)
  file_move, FIELDS_FILE, FIELDS_FILE+".bak", /OVERWRITE
  ; Process the outlist to get the new shortfields (they include mocks)
  for i=0,n_elements(outlist)-1 do begin
    fn = file_basename(outlist[i])
    ; Get the fieldmock (FXMN) and the field (FN) of each outlist file
    prefix = strsplit(fn, "-", /EXTRACT)
    fieldmock = prefix[0]
    prefix = strsplit(fieldmock, 'M', /EXTRACT)
    field = prefix[0]
   
    ; Check if the field (FX) is already included in the new fields 
    newidx = where(newshortfields eq field)
    if newidx lt 0 then begin
      ; It is not inclueded, add it!!
      oldidx = where(oldshortfields eq field)
      if oldidx ge 0 then begin
        ; Add both shortfield and field
        PUSH,newshortfields,field
        PUSH,newfields,oldfields[oldidx]
      endif else printlog,logfile,"WARNING: NO field " + field + " found in original field file!!"
    endif
       
    ; Check if the fieldmock (FXMN) is already included in the new fields 
    newidx = where(newshortfields eq fieldmock)
    if newidx lt 0 then begin
      ; It is not inclueded, add it!!
      oldidx = where(oldshortfields eq field)
      if oldidx ge 0 then begin
        ; Add both shortfield and field
        PUSH,newshortfields,fieldmock
        PUSH,newfields,oldfields[oldidx]
      endif
    endif
  endfor

  ; Write to file...
  out = newshortfields+'   '+newfields
  WRITELINE,FIELDS_FILE,out
endif else printlog,logfile,"WARNING: NO fields field found!!"
  
  

; --------------------------------------
; OLD CODE TO GENERATE LIST (using shell scripts) 
; --------------------------------------
; Update the list with all NEW mch files create for each mock (F*M*....mch)
;cmd = 'find `pwd` -name "F*M*-*_??.mch" > logs/ALLFRAME.inlist'
;ec=runshellcmd(cmd,msg='ALLFRAME.inlist might not have been properly created!')

; Make a copy of old apcor.lst and create a new file with all info from partial files created by Python
;ec=runshellcmd('mv -f apcor.lst apcor.lst.bak', /quiet)
;cmd = 'find . -type f -name "*partial_apcor*.lst" -exec cat {} >> apcor.lst \;'
;ec=runshellcmd(cmd,msg='Info about Aperture Correction in apcor.lst file needed in CALIB stage might NOT be valid!!')

; Make a copy of old calib.trans and create a new file with all info from partial files created by Python
;ec=runshellcmd('mv -f calib.trans calib.trans.bak', /quiet)
;cmd = 'find . -type f -name "*partial_calib*.trans" -exec cat {} >> calib.trans \;'
;ec=runshellcmd(cmd,msg='Info about Transformation Equations in calib.trans file needed in CALIB stage might NOT be valid!!')

printlog,logfile,'FAKERED_ADDSTAR Finished  ',systime(0)

if keyword_set(stp) then stop

end
