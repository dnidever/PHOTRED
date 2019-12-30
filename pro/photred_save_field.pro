;+
;
; PHOTRED_SAVE_SAVE
;
; This saves the final files for photred for a single field.
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

pro photred_save_field,input,redo=redo,sumquick=sumquick,stp=stp

COMMON photred,setup

CD,current=curdir

logfile = -1
t0 = systime(1)

;; Load the "fields" file
;;-------------------
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

longfile = input
file = FILE_BASENAME(longfile)
base = FILE_BASENAME(file,'.dered')
filedir = FILE_DIRNAME(longfile)

CD,filedir

;; Test that it exists
test = FILE_TEST(file)

;; File exists
if (test eq 1) then begin

  ;; Match it with the proper field name
  sep = '-'  ;'.'
  ishortfield = first_el(strsplit(base,sep,/extract))
  gd = where(shortfields eq ishortfield,ngd)

  ;; A field matched
  if (ngd gt 0) then begin
    ifield = fields[gd[0]]

    printlog,logfile,''
    printlog,logfile,'=============================='
    printlog,logfile,file,'  ',ifield
    printlog,logfile,'=============================='
    printlog,logfile,''
    printlog,logfile,systime(0)
      
    ;; Copy the DERED file to FINAL
    finalfile = ifield+'.final'
    printlog,logfile,'Copying ',file,' -> ',finalfile
    FILE_COPY,file,finalfile,/overwrite
    if file_test(file+'.meta') eq 1 then FILE_COPY,file+'.meta',finalfile+'.meta',/overwrite

    ;; Make the IDL SAVE file
    savefile = ifield+'.dat'
    final = PHOTRED_READFILE(file,meta,count=count)
    printlog,logfile,'Making IDL SAVE file ',savefile
    if n_elements(meta) gt 0 then SAVE,final,meta,file=savefile else SAVE,final,file=savefile

    ;; Make FITS binary file
    fitsfile = ifield+'.fits'
    printlog,logfile,'Making FITS binary file ',fitsfile
    MWRFITS,final,fitsfile,/create
    if n_elements(meta) gt 0 then MWRFITS,meta,fitsfile,/silent
    printlog,logfile,'Compressing FITS binary file'
    gfitsfile = fitsfile+'.gz'
    if file_test(gfitsfile) then file_delete,gfitsfile    
    SPAWN,['gzip',fitsfile],/noshell

    ;; Copy final files from FIELD subdirectory to "main" directory
    if filedir ne curdir then begin
      printlog,logfile,'Copying final files to main directory ',curdir
      FILE_COPY,[finalfile,savefile,gfitsfile],curdir+'/'+[finalfile,savefile,gfitsfile],/allow_same,/overwrite
      if file_test(finalfile+'.meta') eq 1 then FILE_COPY,finalfile+'.meta',curdir,/allow_same,/overwrite
      ;; Rename the final output files
      finalfile = curdir+'/'+finalfile
      savefile = curdir+'/'+savefile
      gfitsfile = curdir+'/'+gfitsfile
    endif

    ;; Check that the FINAL and DAT files exist
    finaltest = FILE_TEST(finalfile)
    savetest = FILE_TEST(savefile)
    fitstest = FILE_TEST(gfitsfile)

    ;; Make Field Summary file
    PHOTRED_FIELDSUMMARY,ishortfield,setupdir=curdir,redo=redo,quick=sumquick

  ;; No field match
  endif else begin
    printlog,logfield,'NO field matched for ',field
  endelse


;; File NOT FOUND
endif else begin
  printlog,logfile,'File ',file,' NOT FOUND'
endelse

CD,curdir


printlog,logfile,'dt = '+strtrim(systime(1)-t0,2)+' sec.'

if keyword_set(stp) then stop

end
