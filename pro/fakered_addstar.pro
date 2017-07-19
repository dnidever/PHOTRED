;+
;
; FAKERED_ADDSTAR
;
; This adds artificial stars to the mock images.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The artificial stars are added to the mock images.
;
; By D.Nidever  March 2017
;   Antonio Dorta,  major update   June/July 2017
;-

pro fakered_addstar,redo=redo,stp=stp

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING FAKERED_ADDSTAR'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'ADDSTAR'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
inputfile = 'logs/'+thisprog+'.inlist'
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting FAKERED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readpar', 'printlog', 'undefine', 'check_python', 'check_pyhtcondor', 'push', $
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

; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Separate field directories
;                       VARIABLE       NAME IN SETUP    SETUP  DEFAULT VALUE   LOGFILE       
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
if datatransfer ne 'skip' and datatransfer ne 'copy' and datatransfer ne 'move' and datatransfer ne 'link' then datatransfer='skip'

; Check if the scripts exist in the current directory
scripts = ['fakered_transfer.sh']
nscripts = n_elements(scripts)
; Loop through the scripts
for i=0,nscripts-1 do begin
  fullpath = scriptsdir+'/'+scripts[i]

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

; Check that the DAOMSATER programs exist
SPAWN,'which daomaster',out,errout
daomasterfile = FILE_SEARCH(out,count=ndaomasterfile)
if (ndaomasterfile eq 0) then begin
  print,'DAOMASTER PROGRAM NOT AVAILABLE'
  return
endif


;##########################################################
;#  PROCESSING THE FILES
;##########################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''


;----------------------
; RUNNING THE COMMANDS
;----------------------
printlog,logfile,''
printlog,logfile,systime(0)

if datatransfer ne 'skip' then begin
  cmd = scriptsdir + '/fakered_transfer.sh'
  params = " '" + datatransfer + "' '" + datadir + "' '" + strcompress(sepfielddir) + "' '" + strcompress(sepchipdir) $
               + "' '" + starsfile + "' '" + strcompress(starssplit)  + "' '" + chipsfile + "'"

  print,cmd+params
  SPAWN,cmd+params,out,errout,EXIT_STATUS=exitcode

  if exitcode ne 0 then begin
    ; Check errors (exitcode is NOT 0) 
    printlog,logfile,'' 
    printlog,logfile,"#############"
    printlog,logfile,"# ERROR!!   #"
    printlog,logfile,"#############"
    printlog,logfile,"There was an error when executing script", cmd
    printlog,logfile,"Error was:", errout
    printlog,logfile,'' 
    return
  endif
  print,out
endif

; READ FIELDS AND CHIPS
readline, 'addstar_fields', fields, count=nfields
readline, 'addstar_chips' , chips,  count=nchips


;--------------------------------
; PROCESS ALL FIELDS AND CHIPS
;--------------------------------

undefine,cmd,cmddir
htcondor_cmd = ""
reldir      = ""
base_chip   = ""
base_script = ""
base_dir    = "."  
workdir     = ""
CD,current= curdir

; Check if we have extra directories (for fields and/or chips)
if sepfielddir eq 1 then reldir += "../"
if sepchipdir  eq 1 then reldir += "../"

for fld=0,nfields-1 do begin
  for chp=0,nchips-1 do begin
   
    ; Build workdir according to extra directories 
    workdir_tmp = ""
    if sepfielddir eq 1 then workdir_tmp += "/"+fields[fld]
    if sepchipdir  eq 1 then workdir_tmp += "/chip"+chips[chp]

    if htcondor eq '0' or htcondor ne 'transfer_files' then begin
    	; Using NO HTCondor OR shared directory: change directory to workdir
      base_dir = curdir+workdir_tmp
      base_chip = reldir
      base_script = scriptsdir + "/"
   endif else begin
  		; Get and store working dir needed by HTCondor when transferring files
      workdir += "," + workdir_tmp
    endelse
 
    ; Make commands for addstar python program
    cmd1  = pythonbin +  " "  + base_script + "fakered_addstar.py '"  $
          + string(fields[fld])      + "' '" + string(chips[chp])      + "' '"       $
          + base_chip                        + string(chipsfile)       + "' '"       $
          + strcompress(starscols)   + "' '" + strcompress(magext)     + "' '"       $
          + strcompress(maxccdsize)  + "' '" + strcompress(radcent)    + "' '"       $
          + strcompress(dimfield)    + "' '" + strcompress(distance)   + "'"

    ; Add commands and directories
    PUSH,cmd,cmd1
    PUSH,cmddir,base_dir

  endfor
endfor


if htcondor ne '0' then begin
  ; HTCondor is ENABLED!!

  ; Check whether we are going to use shared directory or transfer files
  if htcondor eq 'transfer_files' then htcondor_cmd = "#transfer_files" $
  else htcondor_cmd = "#shared"

  ; Add common commands: Rank (choose first the fastest machines) and 
  ; Set basic environment variables
  ;htcondor_cmd += "|;|Rank|=|Kflops"
  htcondor_cmd += "|;|environment|=|PATH="+getenv('PATH')+";LD_LIBRARY_PATH="+ $
                  getenv('LD_LIBRARY_PATH')+";HOME="+getenv('HOME')


  if htcondor eq "transfer_files" then begin
    ; Only when transferring files, we need to specify which those files are:

    ; Locate chips file, since given name could include wildcards (use full path)
    chips_files = file_search(chipsfile,/FULLY_QUALIFY_PATH)
    if n_elements(chips_files) ge 1 then fn_chipsfile = chips_files[0] $
    else begin
      print,"ERROR: No chips file were found!!"
      return
    endelse
    if n_elements(chips_files) gt 1 then print,"WARNING!! More than 1 chip file was found. Using first one!!"
    ; Add request disk (at least 5GB), initial workdir and which 
    ; files to transfer: python script, chips files and ALL files in workdir (./)
    htcondor_cmd += "|;|Request_disk|=|5GB"
    htcondor_cmd += "|;|Initialdir|=|$CHOICE(Process"+workdir+")" 
    htcondor_cmd += "|;|transfer_input_files|=|"+scriptsdir+"/"+"fakered_addstar.py" +","+ fn_chipsfile +",./"
  endif

  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='py',hyperthread=hyperthread,waittime=15,    $
           htcondor=htcondor_cmd, scriptsdir=scriptsdir, pythonbin=pythonbin

endif else $
  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='py',hyperthread=hyperthread,waittime=15,    $
           htcondor=htcondor_cmd, scriptsdir=scriptsdir, pythonbin=pythonbin,/cdtodir



; UPDATE the Lists
;PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
;                    failurelist=failurelist,/silent
; Update the list with all NEW mch files create for each mock (F*M*....mch)
cmd = 'find `pwd` -name "F*M*-*_??.mch" > logs/ALLFRAME.inlist'
print,cmd
SPAWN,cmd,out,errout,EXIT_STATUS=exitcode

if exitcode ne 0 then begin
  ; Check errors (exitcode is NOT 0) 
  printlog,logfile,'' 
  printlog,logfile,"#############"
  printlog,logfile,"# ERROR!!   #"
  printlog,logfile,"#############"
  printlog,logfile,"There was an error when executing script", cmd
  printlog,logfile,"Error was:", errout
  printlog,logfile,'' 
  return
endif



; Make a copy of old apcor.lst and create a new file with all info from partial files
; that were created by Python script
cmd = 'mv -f apcor.lst apcor.lst.bak' 
print,cmd
SPAWN,cmd,out,errout,EXIT_STATUS=exitcode
cmd = 'find . -type f -name "*partial_apcor*.lst" -exec cat {} >> apcor.lst \;'
print,cmd
SPAWN,cmd,out,errout,EXIT_STATUS=exitcode
if exitcode ne 0 then begin
  ; Check errors (exitcode is NOT 0) 
  printlog,logfile,'' 
  printlog,logfile,"#############"
  printlog,logfile,"# WARNING!!   #"
  printlog,logfile,"#############"
  printlog,logfile,"There was an error when executing script", cmd
  printlog,logfile,"Error was:", errout
  printlog,logfile,"Info about Aperture Correction needed in CALIB stage might NOT be valid!!"
  printlog,logfile,'' 
endif

cmd = 'mv -f calib.trans calib.trans.bak' 
print,cmd
SPAWN,cmd,out,errout,EXIT_STATUS=exitcode
cmd = 'find . -type f -name "*partial_calib*.trans" -exec cat {} >> calib.trans \;'
print,cmd
SPAWN,cmd,out,errout,EXIT_STATUS=exitcode
if exitcode ne 0 then begin
  ; Check errors (exitcode is NOT 0) 
  printlog,logfile,'' 
  printlog,logfile,"#############"
  printlog,logfile,"# WARNING!!   #"
  printlog,logfile,"#############"
  printlog,logfile,"There was an error when executing script", cmd
  printlog,logfile,"Error was:", errout
  printlog,logfile,"Info about Transformation Equations needed in CALIB stage might NOT be valid!!"
  printlog,logfile,'' 
endif


printlog,logfile,'FAKERED_ADDSTAR Finished  ',systime(0)

if keyword_set(stp) then stop

end
