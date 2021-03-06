;+
;
; PHOTRED_SAVE
;
; This saves the final files for photred
; see post_processingsmc.pro
;
; INPUTS:
;  /redo      Redo files that were already done.
;  /sumquick  Create the summary file quickly.
;  /stp       Stop at the end of the program.
;
; OUTPUTS:
;  The final calibrated photometry and astrometry files.
;
; By D.Nidever  Mar 2008
;-

pro photred_save,redo=redo,sumquick=sumquick,nmulti=nmulti,stp=stp

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_SAVE'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'SAVE'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','photred_getinput','photred_updatelists',$
         'photred_loadsetup','push','undefine','printlog','importascii','first_el','strsplitter',$
         'touchzero','writeline','mktemp','stress','strep']
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


;; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

;; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')

;; Catalog format to use
catformat = READPAR(setup,'catformat')
if catformat eq '0' or catformat eq '' or catformat eq '-1' then catformat='ASCII'
if catformat ne 'ASCII' and catformat ne 'FITS' then catformat='ASCII'

;; Quick summary
if n_elements(sumquick) eq 0 then begin
  sumquick = READPAR(setup,'sumquick')
  if sumquick eq '0' or sumquick eq '' or sumquick eq '-1' then sumquick=0
endif

;; Getting NMULTI from setup file if not given on command line
if n_elements(nmulti) eq 0 then begin
  nmulti_setup = READPAR(setup,'NMULTI')
  if nmulti_setup ne '0' and nmulti_setup ne '' and nmulti_setup ne '-1' then nmulti=long(nmulti_setup)

  ;; Use NMULTI_SAVE if set
  nmultisave = READPAR(setup,'NMULTI_SAVE')
  if nmultisave ne '0' and nmultisave ne '' and nmultisave ne '-1' then nmulti=long(nmultisave)
  if n_elements(nmulti) eq 0 then nmulti=5  ; by default
  nmulti = nmulti > 1  ; must be >=1
endif

;; Hyperthread?
hyperthread = READPAR(setup,'hyperthread',count=nhyperthread)
if nhyperthread gt 0 and hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if nhyperthread gt 0 and strtrim(hyperthread,2) eq '0' then hyperthread=0
if nhyperthread eq 0 then hyperthread=1  ; use by default

;; Clean intermediate files at the end
clean = READPAR(setup,'CLEAN',count=nclean)
if nclean eq 0 then undefine,clean


;###################
; GETTING INPUTLIST
;###################
; INLIST         DERED files
; OUTLIST        FINAL files
; SUCCESSLIST    DERED files

; Get input
;-----------
precursor = 'DEREDDEN'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='dered')
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

;; Run PHOTRED_SAVE_FIELD.PRO
cmd = "photred_save_field,'"+inputlines+"'"
if keyword_set(sumquick) then cmd+=',/sumquick'
if keyword_set(redo) then cmd+=',/redo'
cmddir = file_dirname(inputlines)
PBS_DAEMON,cmd,cmddir,jobs=jobs,/idle,hyperthread=hyperthread,nmulti=nmulti,prefix='sav',wait=1

;; Check the outputs
undefine,outlist,successlist,failurelist
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  base = FILE_BASENAME(file,'.dered')
  filedir = FILE_DIRNAME(longfile)

  CD,filedir

  ; Test that it exists
  test = FILE_TEST(file)

  ; File exists
  if (test eq 1) then begin

    ; Match it with the proper field name
    sep = '-'  ;'.'
    ishortfield = first_el(strsplit(base,sep,/extract))
    gd = where(shortfields eq ishortfield,ngd)

    ;; A field matched
    if (ngd gt 0) then begin
      ifield = fields[gd[0]]      
      finalfile = ifield+'.final'
      savefile = ifield+'.dat'
      fitsfile = ifield+'.fits'
      gfitsfile = fitsfile+'.gz'

      ;; Copy final files from FIELD subdirectory to "main" directory
      if filedir ne curdir then begin
        ; Rename the final output files
        finalfile = curdir+'/'+finalfile
        savefile = curdir+'/'+savefile
        gfitsfile = curdir+'/'+gfitsfile
      endif

      ;; Check that the FINAL and DAT files exist
      finaltest = FILE_TEST(finalfile)
      savetest = FILE_TEST(savefile)
      fitstest = FILE_TEST(gfitsfile)

      ;; We were successful
      if (finaltest eq 1) and (savetest eq 1) and (fitstest eq 1) then begin
        PUSH,outlist,[finalfile,savefile,gfitsfile]   ; add all files to outputarr
        PUSH,successlist,longfile
      end else begin
        PUSH,failurelist,longfile
        if finaltest eq 0 then printlog,logfile,finalfile,' NOT FOUND'
        if savetest eq 0 then printlog,logfile,savefile,' NOT FOUND'
        if savetest eq 0 then printlog,logfile,gfitsfile,' NOT FOUND'
      endelse

    ;; No field match
    endif else begin
      PUSH,failurelist,longfile
      printlog,logfield,'NO field matched for ',field
    endelse


  ;; File NOT FOUND
  endif else begin
    PUSH,failurelist,longfile
    printlog,logfile,'File ',file,' NOT FOUND'
  endelse

  CD,curdir


  ;#####################
  ; UPDATE the Lists
  ;#####################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,setupdir=curdir,/silent

ENDFOR


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


;;######################
;;  CLEANING UP
;;######################
if keyword_set(clean) then begin
  printlog,logfile,'CLEANING UP.  CLEAN='+strtrim(clean,2)

  ;; Remove FITS files that have resource files
  ;;  only successful ones
  READLIST,curdir+'/logs/DAOPHOT.success',fitsfiles,/unique,/fully,setupdir=curdir,count=nfitsfiles,logfile=logfile,/silent
  for i=0,nfitsfiles-1 do begin
    dir1 = file_dirname(fitsfiles[i])
    base1 = file_basename(fitsfiles[i])
    rfile = dir1+'/.'+base1
    info = file_info(fitsfiles[i])    
    rinfo = file_info(rfile)
    if info.exists eq 1 and info.size gt 1 and rinfo.exists eq 1 then begin
      FILE_DELETE,fitsfiles[i],/allow
      WRITELINE,fitsfiles[i],''   ;; create size=1 file
    endif
  endfor
endif

printlog,logfile,'PHOTRED_SAVE Finished  ',systime(0)

if keyword_set(stp) then stop

end
