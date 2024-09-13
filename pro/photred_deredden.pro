;+
;
; PHOTRED_DEREDDEN
;
; This dereddens photometry for photred
; See deredden.pro
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The final calibrated and dereddened photometry files
;
; By D.Nidever  Mar 2008
;-
pro photred_deredden,nmulti=nmulti,redo=redo,stp=stp


COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_DEREDDEN'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'DEREDDEN'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','first_el','strsplitter','combine_structs',$
         'photred_getinput','photred_updatelists','photred_loadsetup','dust_getval','touchzero',$
         'wcs_getval','wcs_coord2pix','push','undefine','printlog','add_tag','printstr','writeline',$
         'bh_rdfort','djs_int2bin','mktemp','stress','wcs_coord2pix','djs_angpos','djs_ceil','strep']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; Check that DUST_DIR is properly defined in the shell
dustdir = GETENV('DUST_DIR')
if (dustdir eq '') then begin
  print,'DUST_DIR TO SCHLEGEL MAPS NOT DEFINED IN SHELL'
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
; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')
; CATFORMAT
catformat = READPAR(setup,'catformat')
if catformat eq '0' or catformat eq '' or catformat eq '-1' then catformat='ASCII'
if catformat ne 'ASCII' and catformat ne 'FITS' then catformat='ASCII'
; Magnitudes and colors to deredden
todered = READPAR(setup,'TODERED')
if todered eq '-1' or todered eq '0' then undefine,todered
ntodered = n_elements(todered)
; Extinction values to add
toextadd = READPAR(setup,'TOEXTADD')
if toextadd eq '-1' or toextadd eq '0' then undefine,toextadd
ntoextadd = n_elements(toextadd)


;###################
; GETTING INPUTLIST
;###################
; INLIST         CMB files
; OUTLIST        DERED files
; SUCCESSLIST    CMB files

; Get input
;-----------
precursor = 'COMBINE'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='cmb')
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines


; Load the "extinction" file
;--------------------------
printlog,logfile,'Loading the >>extinction<< file'
extfile = curdir+'/extinction'
test = FILE_TEST('extinction')
if (test eq 0) then begin
  scriptsdir = READPAR(setup,'SCRIPTSDIR')
  if scriptsdir eq '-1' or scriptsdir eq '0' then begin
    printlog,logfile,'NO SCRIPTSDIR'
    return
  endif
  test2 = FILE_TEST(scriptsdir+'/extinction')
  if (test2 eq 0) then begin
    printlog,logfile,'NO >>extinction<< file in ',scriptsdir
    return
  endif
  FILE_COPY,scriptsdir+'/extinction','.',/overwrite
  extfile = scriptsdir+'/extinction'
endif

; Load the extinction values
; EXTRATIO is A(filter)/E(B-V), so A(filter) = EXTRATIO * E(B-V)
extstr = IMPORTASCII(extfile,fieldname=['filter','extratio'],comment='#',/noprint)
nextstr = n_elements(extstr)
if (nextstr eq 0) then begin
  printlog,logfile,'>>extinction<< file is EMPTY'
  return
endif

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

; Looping through the input files
FOR i=0,ninputlines-1 do begin
  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  ending = first_el(strsplit(file,'.',/extract),/last)
  base = FILE_BASENAME(file,'.'+ending)  

  cmd1 = "PHOTRED_DEREDDEN_SINGLE,'"+longfile+"','"+extfile+"'"
  if ntodered gt 0 then cmd1 += ",todered='"+todered+"'"
  if ntoextadd gt 0 then cmd1 += ",toextadd='"+toextadd+"'"
  cmd1 += ",catformat='"+catformat+"'"
  
  PUSH,cmd,cmd1
  PUSH,cmddir,filedir
ENDFOR
ncmd = n_elements(cmd)

if ncmd gt 0 then begin
  cmd = "cd,'"+cmddir+"' & "+cmd ; go to the directory
  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='dered',hyperthread=hyperthread,/idle,$
             waittime=1,scriptsdir=scriptsdir
endif

;; Check that the outfiles exist
for i=0,ninputlines-1 do begin
  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  ending = first_el(strsplit(file,'.',/extract),/last)
  base = FILE_BASENAME(file,'.'+ending)  
  deredfile = filedir+'/'+base+'.dered'
  ;; Check that the file dered file is there
  deredtest = FILE_TEST(deredfile)
  if (deredtest eq 1) then begin
    PUSH,outlist,deredfile
    PUSH,successlist,longfile
  endif else begin
    PUSH,failurelist,longfile
    printlog,logfile,deredfile,' NOT FOUND'
  endelse
endfor


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


printlog,logfile,'PHOTRED_DEREDDEN Finished  ',systime(0)

if keyword_set(stp) then stop

end
