;+
;
; PHOTRED_SAVE
;
; This saves the final files for photred
; see post_processingsmc.pro
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The final calibrated photometry and astrometry files.
;
; By D.Nidever  Mar 2008
;-

pro photred_save,redo=redo,stp=stp

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


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

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


; Rename and save the files

undefine,outlist,successlist,failurelist

; Loop through all input files
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

    ; A field matched
    if (ngd gt 0) then begin

      ifield = fields[gd[0]]

      printlog,logfile,''
      printlog,logfile,'=============================='
      printlog,logfile,file,'  ',ifield
      printlog,logfile,'=============================='
      printlog,logfile,''
      printlog,logfile,systime(0)
      
      ; Copy the DERED file to FINAL
      finalfile = ifield+'.final'
      printlog,logfile,'Copying ',file,' -> ',finalfile
      FILE_COPY,file,finalfile,/overwrite
      if file_test(file+'.meta') eq 1 then FILE_COPY,file+'.meta',finalfile+'.meta'

      ; Make the IDL SAVE file
      savefile = ifield+'.dat'
      final = PHOTRED_READFILE(file,meta=meta,count=count)
      printlog,logfile,'Making IDL SAVE file ',savefile
      if n_elements(meta) gt 0 then SAVE,final,meta,file=savefile else SAVE,final,file=savefile

      ; Make FITS binary file
      fitsfile = ifield+'.fits'
      printlog,logfile,'Making FITS binary file ',fitsfile
      MWRFITS,final,fitsfile,/create
      if n_elements(meta) gt 0 then MWRFITS,meta,fitsfile,/silent
      printlog,logfile,'Compressing FITS binary file'
      gfitsfile = fitsfile+'.gz'
      if file_test(gfitsfile) then file_delete,gfitsfile    
      SPAWN,['gzip',fitsfile],/noshell

      ; Copy final files from FIELD subdirectory to "main" directory
      if filedir ne curdir then begin
        printlog,logfile,'Copying final files to main directory ',curdir
        FILE_COPY,[finalfile,savefile,gfitsfile],curdir+'/'+[finalfile,savefile,gfitsfile],/allow_same,/overwrite
        if file_test(finalfile+'.meta') eq 1 then FILE_COPY,finalfile+'.meta',curdir,/allow_same,/overwrite
        ; Rename the final output files
        finalfile = curdir+'/'+finalfile
        savefile = curdir+'/'+savefile
        gfitsfile = curdir+'/'+gfitsfile
      endif

      ; Check that the FINAL and DAT files exist
      finaltest = FILE_TEST(finalfile)
      savetest = FILE_TEST(savefile)
      fitstest = FILE_TEST(gfitsfile)

      ; We were successful
      if (finaltest eq 1) and (savetest eq 1) and (fitstest eq 1) then begin
        PUSH,outlist,filedir+'/'+[finalfile,savefile,gfitsfile]   ; add all files to outputarr
        PUSH,successlist,longfile
      end else begin
        PUSH,failurelist,longfile
        if finaltest eq 0 then printlog,logfile,finalfile,' NOT FOUND'
        if savetest eq 0 then printlog,logfile,savefile,' NOT FOUND'
        if savetest eq 0 then printlog,logfile,gfitsfile,' NOT FOUND'
      endelse

      ; Make Field Summary file
      PHOTRED_FIELDSUMMARY,ishortfield,setupdir=curdir,redo=redo

    ; No field match
    endif else begin
      PUSH,failurelist,longfile
      printlog,logfield,'NO field matched for ',field
    endelse


  ; File NOT FOUND
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


printlog,logfile,'PHOTRED_SAVE Finished  ',systime(0)

if keyword_set(stp) then stop

end
