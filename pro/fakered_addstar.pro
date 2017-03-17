pro fakered_addstar,redo=redo,stp=stp

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
;-

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



;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS files, it get files from directory
; OUTLIST        FITS files
; SUCCESSLIST    FITS files
undefine,outlist,successlist,failurelist

; Add all fits files to the INLIST
fitsfiles = FILE_SEARCH('*.fits',count=nfitsfiles,/fully)

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

; These are the files to process
procbaselist = FILE_BASENAME(inputlines,'.fits')
procdirlist = FILE_DIRNAME(inputlines)


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
printlog,logfile,'Running addstar python script on '+strtrim(n_elements(procbaselist),2)+' files'
printlog,logfile,''
printlog,logfile,systime(0)

; Make commands for addstar python program
cmd = './addstar.py '+procbaselist
  
; Submit the jobs to the daemon
PBS_DAEMON,cmd,procdirlist,nmulti=nmulti,prefix='dao',hyperthread=hyperthread,waittime=30,/cdtodir



; UPDATE the Lists
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent


printlog,logfile,'FAKERED_ADDSTAR Finished  ',systime(0)

if keyword_set(stp) then stop

end
