;+
;
; STDRED_FITDATA
;
; This uses the combined catalogs for each filter and
; fits the transformation equations with FIT_TRANSPHOT.PRO.
; A separate transformation equation is output for each
; filter and for each night.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  A separate transformation equations ".trans" file is
;  output for each filter and for each night.
;
; By D.Nidever  May 2008
;-

pro stdred_fitdata,redo=redo,stp=stp

COMMON photred,setup

print,''
print,'#########################'
print,'RUNNING STDRED_FITDATA'
print,'#########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'FITDATA'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting STDRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','add_tag','printstr','photred_getinput',$
         'photred_updatelists','photred_loadsetup','push','undefine','printlog','stdred_transphot',$
         'writeline','mad','mktemp','mpfit','range','strsplitter','printline','touchzero','first_el',$
         'stress','combine_structs','strep']
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
  PHOTRED_LOADSETUP,setup,count=count,/std
  if count lt 1 then return
endif


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')



;###################
; GETTING INPUTLIST
;###################
; INLIST         PHOT files
; OUTLIST        AST files
; SUCCESSLIST    PHOT files

; Get input
;-----------
precursor = 'COMBINECAT'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='cat')
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines


; ERRLIM
fiterrlim = READPAR(setup,'FITERRLIM')
if fiterrlim eq '0' or fiterrlim eq '' or fiterrlim eq '-1' then undefine,fiterrlim
; SEPCHIP
sepchip = READPAR(setup,'FITSEPCHIP')
if sepchip eq '0' or sepchip eq '' or sepchip eq '-1' then undefine,sepchip else sepchip=1


;##################################################
;#  PROCESSING THE FILES
;##################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'

undefine,outlist,successlist,failurelist,transarr

; Loop through the combined CAT files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  base = FILE_BASENAME(file,'.cat')


  printlog,logfile,''
  printlog,logfile,'================================================='
  printlog,logfile,' Fitting transformation equation for ',file
  printlog,logfile,'================================================='


  ; Test that the CAT file exists
  test = FILE_TEST(file)
  if test eq 1 then nlines = FILE_LINES(file)
  if (test eq 0) or (nlines eq 0) then begin
    PUSH,failurelist,longfile
    if test eq 0 then printlog,logfile,file,' NOT FOUND'
    if test eq 1 and nlines eq 0 then printlog,logfile,file,' HAS 0 LINES'
    goto,BOMB
  endif


  ; Run STDRED_TRANSPHOT.PRO
  printlog,logfile,''
  printlog,logfile,'Running STDRED_TRANSPHOT'
  undefine,tlines,trans,rms
  STDRED_TRANSPHOT,longfile,tlines=tlines,trans=trans,rms=rms,arr=arr,errlim=fiterrlim,$
                   sepchip=sepchip
  ntrans = n_elements(trans)

  ; Success
  if (ntrans gt 0) then begin
    PUSH,successlist,longfile
    PUSH,transarr,trans

  ; Failure
  endif else begin
    printlog,logfile,transfile,' NOT FOUND'
    PUSH,failurelist,longfile
  endelse


  BOMB:


  ;##########################################
  ;#  UPDATING LIST FILES
  ;##########################################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent


END ; filter loop


; Make transformation equations for each night
ui = uniq(transarr.night,sort(transarr.night))
nights = transarr[ui].night
nnights = n_elements(nights)

; Loop through each night
for i=0,nnights-1 do begin

  inight = nights[i]
  printlog,logfile,''
  printlog,logfile,'------------------------------------------------'
  printlog,logfile,'Making transformation equation file for NIGHT=',strtrim(inight,2)
  printlog,logfile,'------------------------------------------------'
  printlog,logfile,''

  gd = where(transarr.night eq inight,ngd)
  itrans = transarr[gd]

  mags = transarr[gd].magname
  printlog,logfile,'Trans. Eqns. for ',strtrim(ngd,2),' magnitude: ',strjoin(mags,', ')

  ; Loop through the elements and add the trans lines
  undefine,transout
  printlog,logfile,''
  for j=0,ngd-1 do begin
    PUSH,transout,itrans[j].magline
    PUSH,transout,itrans[j].errline
    printlog,logfile,itrans[j].magline
    printlog,logfile,itrans[j].errline
  end
  printlog,logfile,''


  ; Write the nightly trans file
  transfile = 'n'+strtrim(inight,2)+'.trans'
  printlog,logfile,'Writing trans file "',transfile,'"'
  WRITELINE,transfile,transout

end



;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'STDRED_FITDATA Finished  ',systime(0)

if keyword_set(stp) then stop


end
