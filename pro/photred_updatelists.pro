;+
;
; PHOTRED_UPDATELISTS
;
; This updates the lists (inputlist, outputlist, successlist and failurelist)
; for a PHOTRED stage.
;
; INPUTS:
;  lists          The lists structure output by photred_getinput.
;  =outlist       The list of new output files.
;  =successlist   The list of new successfully processed files.
;  =failurelist   The list of new failed files.
;  =setupdir      The main directory where photred.setup exists.
;  /silent        Don't print anything
;
; OUTPUTS:
;  The stage's lists will be updated
;
; USAGE:
;  IDL>photred_updatelists,lists,outlist=outlist,successlist=successlist,$
;                          failurelist=failurelist,setupdir=setupdir,silent=silent,stp=stp
;
; By D.Nidever  March 2008
;-

pro photred_updatelists,lists,outlist=outlist,successlist=successlist,$
                        failurelist=failurelist,setupdir=setupdir,silent=silent,stp=stp

; Not enough inputs
nlists = n_elements(lists)
if (nlists eq 0) or (n_elements(setupdir) eq 0) then begin
  print,'Syntax - photred_updatelists,lists,outlist=outlist,successlist=successlist,'
  print,'                             failurelist=failurelist,setupdir=setupdir,silent=silent,stp=stp'
  return
endif


; Check the "lists" structure
type = SIZE(lists,/type)
if (type ne 8) then begin
  print,'LISTS MUST BE A STRUCTURE'
  return
endif
tags = TAG_NAMES(lists)
gdthisprog = where(tags eq 'THISPROG',ngdthisprog)
if (ngdthisprog eq 0) then begin
  print,'LISTS IS MISSING THE "THISPROG" TAG'
  return
endif
gdninputlines = where(tags eq 'NINPUTLINES',ngdninputlines)
if (ngdthisprog eq 0) then begin
  print,'LISTS IS MISSING THE "NINPUTLINES" TAG'
  return
endif
gdnoutputlines = where(tags eq 'NOUTPUTLINES',ngdnoutputlines)
if (ngdthisprog eq 0) then begin
  print,'LISTS IS MISSING THE "NOUTPUTLINES" TAG'
  return
endif
gdnsuccesslines = where(tags eq 'NSUCCESSLINES',ngdnsuccesslines)
if (ngdthisprog eq 0) then begin
  print,'LISTS IS MISSING THE "NSUCCESSLINES" TAG'
  return
endif


; Is "thisprog" a valid stage?
;-------------------------------
thisprog = strupcase(lists.thisprog)
stages = ['RENAME','WCS','SPLIT','DAOPHOT','MATCH','ALLFRAME','APCOR',$
          'CALIB','ASTROM','COMBINE','DEREDDEN','SAVE','HTML',        $
          'APERPHOT','DAOGROW','MATCHCAT','COMBINECAT','FITDATA',     $ ; STDRED  stages
          'COMPLETE','ADDSTAR','CLEAN',$                                ; FAKERED stages
          'ZEROPOINT']                                                  ; delvered stages
; NOT a valid stage
stageind = where(stages eq thisprog,nstageind)
if (nstageind eq 0) then begin
  print,thisprog,' IS NOT A VALID STAGE'
  return
endif


; List filenames
logfile = 'logs/'+lists.thisprog+'.log'
inputfile = 'logs/'+lists.thisprog+'.inlist'
outputfile = 'logs/'+lists.thisprog+'.outlist'
successfile = 'logs/'+lists.thisprog+'.success'
failurefile = 'logs/'+lists.thisprog+'.failure'


;##########################################
;#  UPDATING LIST FILES
;##########################################
if not keyword_set(silent) then begin
  printlog,logfile,''
  printlog,logfile,'-----------------------'
  printlog,logfile,'UPDATING THE LISTS'
  printlog,logfile,'-----------------------'
  printlog,logfile,''
endif


; Inputlines
undefine,inputlines
ninputlines = lists.ninputlines
if (ninputlines) gt 0 then inputlines=lists.inputlines

;----------------------
; Updating SUCCESS list
;----------------------
nsuccesslist = n_elements(successlist)
if not keyword_set(silent) then $
  printlog,logfile,strtrim(nsuccesslist,2),' files successfully processed'
if (nsuccesslist gt 0) then begin

  ; Adding to the OLD successlist
  undefine,oldsuccesslines,newsuccesslines
  if (lists.nsuccesslines gt 0) then oldsuccesslines = lists.successlines
  PUSH,newsuccesslines,oldsuccesslines
  PUSH,newsuccesslines,successlist

  ; Get unique ones
  ui = UNIQ(newsuccesslines,sort(newsuccesslines))
  ui = ui[sort(ui)]
  newsuccesslines = newsuccesslines[ui]
  nnewsuccesslines = n_elements(newsuccesslines)

  ;; Make the paths relative to SETUPDIR
  flen = strlen(setupdir)
  if strmid(setupdir,0,1,/reverse) ne '/' then flen+=1
  gd = where(stregex(newsuccesslines,'^'+setupdir,/boolean) eq 1,ngd)
  if ngd gt 0 then newsuccesslines[gd] = strmid(newsuccesslines[gd],flen)

  ; Write the success file
  WRITELINE,successfile,newsuccesslines

  ;out = successlist+'  SUCCESS'
  ;printlog,logfile,out,/logonly
endif


;------------------
; Updating OUTLIST
;------------------
noutlist = n_elements(outlist)
if (noutlist gt 0) then begin

  ; Adding to the OLD outlist
  undefine,oldoutputlines,newoutputlines
  if (lists.noutputlines gt 0) then oldoutputlines = lists.outputlines
  PUSH,newoutputlines,oldoutputlines
  PUSH,newoutputlines,outlist

  ; Get unique ones
  ui = UNIQ(newoutputlines,sort(newoutputlines))
  ui = ui[sort(ui)]
  newoutputlines = newoutputlines[ui]
  nnewoutputlines = n_elements(newoutputlines)

  nnew = nnewoutputlines - lists.noutputlines
  if not keyword_set(silent) then begin
    printlog,logfile,strtrim(nnewoutputlines,2),' files in '+thisprog+'.outlist'
    printlog,logfile,strtrim(nnew,2),' files were ADDED to '+thisprog+'.outlist'
  endif

  ;; Make the paths relative to SETUPDIR
  flen = strlen(setupdir)
  if strmid(setupdir,0,1,/reverse) ne '/' then flen+=1
  gd = where(stregex(newoutputlines,'^'+setupdir,/boolean) eq 1,ngd)
  if ngd gt 0 then newoutputlines[gd] = strmid(newoutputlines[gd],flen)

  ; Write the output file
  WRITELINE,outputfile,newoutputlines

  ;out = 'Added '+outlist+' to OUTLIST'
  ;printlog,logfile,out,/logonly
endif


;----------------
; Updating INLIST
;----------------
; Some successful ones
if (nsuccesslist gt 0) then begin
  if not keyword_set(silent) then $
    printlog,logfile,strtrim(nsuccesslist,2),' files removed from '+thisprog+'.inlist'

  ; Some files left over
  if (nsuccesslist lt ninputlines) then begin
    ; Match them
    MATCH,inputlines,successlist,ind1,ind2,count=nmatch

    if (nmatch gt 0) then begin
      newinputlines = inputlines
      REMOVE,ind1,newinputlines

      ;; Make the paths relative to SETUPDIR
      flen = strlen(setupdir)
      if strmid(setupdir,0,1,/reverse) ne '/' then flen+=1
      gd = where(stregex(newinputlines,'^'+setupdir,/boolean) eq 1,ngd)
      if ngd gt 0 then newinputlines[gd] = strmid(newinputlines[gd],flen)

      WRITELINE,inputfile,newinputlines
    endif else begin
      printlog,logfile,'NO MATCH BETWEEN INPUTLINES AND NEW SUCCESSLIST'
    endelse

  ; All succeeded, "empty" inlist
  endif else begin
    ; Make a zero length file
    FILE_DELETE,inputfile,/allow
    TOUCHZERO,inputfile
  endelse
endif


;----------------------
; Updating FAILURE list
;----------------------
nfailurelist = n_elements(failurelist)
if not keyword_set(silent) then $
  printlog,logfile,strtrim(nfailurelist,2),' files FAILED'

; Checking to see if any past failures were successful this time
remnewsuccess = 0
if (lists.nfailurelines gt 0 and nsuccesslist gt 0) then begin
  MATCH,lists.failurelines,successlist,ind1,ind2,count=nind1
  if nind1 gt 0 then remnewsuccess=1
endif

if (nfailurelist gt 0) or (remnewsuccess eq 1) then begin

  ; Adding to the OLD failure list
  undefine,oldfailurelines,newfailurelines
  if (lists.nfailurelines gt 0) then oldfailurelines = lists.failurelines
  PUSH,newfailurelines,oldfailurelines
  PUSH,newfailurelines,failurelist

  ; Get unique ones
  ui = UNIQ(newfailurelines,sort(newfailurelines))
  ui = ui[sort(ui)]
  newfailurelines = newfailurelines[ui]
  nnewfailurelines = n_elements(newfailurelines)

  ; Removing PAST failures that were successful this time
  if (nsuccesslist gt 0) then begin
    MATCH,newfailurelines,successlist,ind1,ind2,count=nind1

    if (nind1 gt 0) then begin
      if nind1 lt nnewfailurelines then REMOVE,ind1,newfailurelines
      if nind1 eq nnewfailurelines then undefine,newfailurelines
      nnewfailurelines = n_elements(newfailurelines)
    endif
  endif ; some successful files

  ;; Make the paths relative to SETUPDIR
  if n_elements(newfailurelines) gt 0 then begin
    flen = strlen(setupdir)
    if strmid(setupdir,0,1,/reverse) ne '/' then flen+=1
    gd = where(stregex(newfailurelines,'^'+setupdir,/boolean) eq 1,ngd)
    if ngd gt 0 then newfailurelines[gd] = strmid(newfailurelines[gd],flen)
  endif

  ; Write the failure file, should be okay even if no lines
  WRITELINE,failurefile,newfailurelines

endif

;stop

if keyword_set(stp) then stop

end
