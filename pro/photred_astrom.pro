;+
;
; PHOTRED_ASTROM
;
; This gets coordinates from the FITS WCS and puts them
; in the photometry files.
;
; INPUTS:
;  =nmulti  Number of parallel jobs to run
;  /redo    Redo files that were already done.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  The final calibrated photometry with accurate astrometry
;
; By D.Nidever  Mar 2008
;-

pro photred_astrom,nmulti=nmulti,redo=redo,stp=stp

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_ASTROM'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'ASTROM'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','head_xyad','hdr2wcstnx','parsetnx','wcstnx_xy2rd',$
         'wcstnxcor','xieta2rd','add_tag','printstr','photred_getinput','photred_updatelists',$
         'photred_loadsetup','push','undefine','printlog','strsplitter','combine_structs','touchzero',$
         'writeline','first_el','mktemp','stress','strep','loadmch']
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


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1
; Minimum number of detections
ndetmin = READPAR(setup,'NDETMIN')
if (ndetmin eq '-1' or ndetmin eq '0' or ndetmin eq '') then undefine,ndetmin else ndetmin=long(ndetmin)

; Hyperthread?
hyperthread = READPAR(setup,'hyperthread')
if hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if strtrim(hyperthread,2) eq '0' then hyperthread=0
if n_elements(hyperthread) eq 0 then hyperthread=0

; Getting NMULTI
if n_elements(nmulti) eq 0 then begin
  nmulti = READPAR(setup,'NMULTI')
  nmulti = long(nmulti)

  ; Use NMULTI_WCS if set
  nmultiwcs = READPAR(setup,'NMULTI_WCS')
  if nmultiwcs ne '0' and nmultiwcs ne '' and nmultiwcs ne '-1' then nmulti=long(nmultiwcs)
endif
nmulti = nmulti > 1  ; must be >=1

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')

; Catalog format to use
catformat = READPAR(setup,'catformat')
if catformat eq '0' or catformat eq '' or catformat eq '-1' then catformat='ASCII'
if catformat ne 'ASCII' and catformat ne 'FITS' then catformat='ASCII'


;###################
; GETTING INPUTLIST
;###################
; INLIST         PHOT files
; OUTLIST        AST files
; SUCCESSLIST    PHOT files

; Get input
;-----------
; Precursors: rename > wcs > split > daophot > match > allframe >
; apcor > *astrom* > calib > combine > deredden > save
precursor = ['ALLFRAME','MATCH']
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo)
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

; Loop through the input PHOT files
FOR i=0,ninputlines-1 do begin
  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)

  cmd1 = "PHOTRED_ASTROM_SINGLE,'"+longfile+"',catformat='"+catformat+"'"
  if keyword_set(redo) then cmd1 += ',/redo'
  if n_elements(ndetmin) gt 0 then cmd1 += ',ndetmin='+strtrim(ndetmin,2)

  PUSH,cmd,cmd1
  PUSH,cmddir,filedir
ENDFOR
ncmd = n_elements(cmd)

if ncmd gt 0 then begin
  cmd = "cd,'"+cmddir+"' & "+cmd  ; go to the directory
  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='astrom',hyperthread=hyperthread,/idle,$
             waittime=1,/cdtodir,scriptsdir=scriptsdir
endif

; Check for success/failures
for i=0,ninputlines-1 do begin
  longfile = inputlines[i]
  file = file_basename(longfile)
  filedir = file_dirname(longfile)
  ending = first_el(strsplit(longfile,'.',/extract),/last)
  if ending eq 'mag' then begin
    base = file_basename(file,'.mag')
  endif else begin
    base = file_basename(file,'.mch')
  endelse 
  astfile = filedir+'/'+base+'.ast'
  ; Check that the file AST file is there
  asttest = FILE_TEST(astfile)
  if (asttest eq 1) then begin
    PUSH,outlist,astfile
    PUSH,successlist,longfile
  endif else begin
    PUSH,failurelist,longfile
    printlog,logfile,astfile,' NOT FOUND'
  endelse
endfor


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


printlog,logfile,'PHOTRED_ASTROM Finished  ',systime(0)

if keyword_set(stp) then stop

end
