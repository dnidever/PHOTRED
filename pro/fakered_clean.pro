;+
;
; FAKERED_CLEAN
;
; Clean up the files.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  No outputs.  Just deleting some of the temporary files.
;
; By D.Nidever  February 2018
;-

pro fakered_clean,redo=redo,stp=stp

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING FAKERED_CLEAN'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'CLEAN'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
inputfile = 'logs/'+thisprog+'.inlist'
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting FAKERED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readpar', 'printlog', 'undefine', 'push', $
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

; PHOTRED parameters
;                       VARIABLE       NAME IN SETUP    SETUP  DEFAULT VALUE   LOGFILE       
redo            = getparam(redo            , 'redo'            , setup, 0             , logfile, /bool)
telescope       = getparam(telescope       , 'telescope'       , setup, '0'           , logfile)
instrument      = getparam(instrument      , 'instrument'      , setup, '0'           , logfile)
hyperthread     = getparam(hyperthread     , 'hyperthread'     , setup, '0'           , logfile, /bool)
nmulti          = getparam(nmulti          , 'nmulti'          , setup, '1'           , logfile)
scriptsdir      = getparam(scriptsdir      , 'scriptsdir'      , setup, ''            , logfile)
sepfielddir     = getparam(sepfielddir     , 'sepfielddir'     , setup, '0'           , logfile, /bool)
sepchipdir      = getparam(sepchipdir      , 'sepchipdir'      , setup, '0'           , logfile, /bool)
starssplit      = getparam(starssplit      , 'starssplit'      , setup, '0'           , logfile, /bool)
;starsfile       = getparam(starsfile       , 'starsfile'       , setup, ''            , logfile)
starscols       = getparam(starscols       , 'starscols'       , setup, ''            , logfile)
datadir         = getparam(datadir         , 'datadir'         , setup, '.'           , logfile)
datatransfer    = getparam(datatransfer    , 'datatransfer'    , setup, 'skip'        , logfile)
chipsfile       = getparam(chipsfile       , 'chipssfile'      , setup, '*chips.fits' , logfile)
maxccdsize      = getparam(maxccdsize      , 'maxccdsize'      , setup, '2048,4096'   , logfile)
magext          = getparam(magext          , 'magext'          , setup, '0'           , logfile)
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
;if starsfile eq '' or file_test(starsfile,/read) ne 1 then begin
;  printlog,logfile,'ERROR: STARSFILE "',starsfile,' NOT FOUND'
;  error='STARSFILE "'+starsfile+' NOT FOUND'
;  return
;endif

; Check if starscols is NOT empty
if starscols eq '' then begin
  printlog,logfile,'ERROR: STARSCOLS NOT DEFINED'
  error='STARSCOLS "'+starsfile+' NOT DEFINED'
  return
endif

; Telescope must be set
if telescope eq '0' then begin
  error = 'TELESCOPE not found in setup file'
  printlog,logfile,error
  return
endif
; Instrument must be set
if instrument eq '0' then begin
  error = 'INSTRUMENT not found in setup file'
  printlog,logfile,error
  return
endif

; Check whether datatransfer is one of the valid methods
if datatransfer ne 'skip' and datatransfer ne 'copy' and datatransfer ne 'move' and datatransfer ne 'link' then datatransfer='skip'

; The main directory
CD,current=maindir

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


; Load "fields" file to get "global" field name
READCOL,maindir+'/fields',field_shnames,field_lnames,format='A,A',/silent

;---------------------------------------------------------------
; Run CHECK_PYTHON.PRO to make sure that you can run PYTHON from IDL
;---------------------------------------------------------------
CHECK_PYTHON,pythontest,pythonbin=pythonbin
if pythontest eq 0 then begin
  print,'PYTHON TEST FAILED.  EXITING'
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

; Get input
;-----------
precursor = 'COMPLETE'
lists = PHOTRED_GETINPUT(thisprog,'COMPLETE.success',redo=redo,ext='phot')
ninputlines = lists.ninputlines

; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines
dirlist = FILE_DIRNAME(inputlines)
baselist = FILE_BASENAME(inputlines,'.phot')
nbaselist = n_elements(baselist)


; Loop over the reference chip files
undefine,outlist,successlist,failurelist
For i=0,nbaselist-1 do begin
  dir1 = dirlist[i]
  base1 = baselist[i]

  mchfile = dir1+'/'+base1+'.mch'
  LOADMCH,mchfile,files
  nfiles = n_elements(files)
  bases = file_basename(files,'.als')
  ; Loop over all of the exposures
  for j=0,nfiles-1 do $
    FILE_DELETE,dir1+'/'+bases[j]+['.add','.fits','.alf'],/allow
  
  ; Reference exposure files
  ;  most of these are ALLFRAME-related files
  FILE_DELETE,dir1+'/'+base1+['.mch','.outlist','.inlist','.shift','.mag','.makemag','_comb.bpm.fits','_comb.fits',$
                              '_comb.mask.fits','_comb.coo','_comb.als.inp','_comb.als','_combs.fits','_comb.log',$
                              '_comb.sex','_comb_sub.cat','_comb_allf.sex','_comb_all.coo','_comb_sub.als','_comb_sub.fits',$
                              '_comb_allf.als','.comb.mch','.nmg','.tfr','.ast','.input','.phot'],/allow

  ; Check that they were successful
  PUSH,successlist,dir1+'/'+base1
  PUSH,outlist,dir1+'/'+base1
  
  BOMB:
Endfor


; UPDATE the Lists
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent

printlog,logfile,'FAKERED_CLEAN Finished  ',systime(0)

if keyword_set(stp) then stop

end
