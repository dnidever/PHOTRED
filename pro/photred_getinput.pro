function photred_getinput,thisprog,precursor,redo=redo,stp=stp,error=error,$
                          extension=extension,noempty=noempty

;+
;
; PHOTRED_GETINPUT
;
; This gets the input for a given PHOTRED stage.  The outputs from the precursor stage
; are copied/moved over and the output list and successlist are checked to make sure
; that a file has not already been processed (unless /redo is set).
;
; INPUTS:
;  thisprog     The stage to get the inputs for (e.g. 'DAOPHOT','MATCH', etc.)
;  precursor    The precursor stage to get outputs from (e.g. 'SPLIT.output','DAOPHOT', etc.).
;                 This can be an array.  If there are no files in the first one, then the next
;                   one is checked, and so on.
;                 If no ending is specified then ".outlist" is used.
;                 The filename should not include "logs/" at the beginning, this will be
;                   automatically prepended.
;                 This is NOT a required input (i.e. RENAME has no precursor).
;  /redo        Redo files already processed
;  /stp         Stop at the end of the program
;  =extension   Only accept inputs with this extension (i.e. 'als').  Do
;                 not include the dot.
;
; OUTPUTS:
;  lists        An IDL structure that includes all of the list information: precursor, prestage,
;                 prefile, inputlines, ninputlines, outputlines, noutputlines, successlines,
;                 and nsuccesslines.
;  =error       The error, if one occured, otherwise undefined.
;
; USAGE:
;  IDL>lists = photred_getinput('DAOPHOT','SPLIT.output')
;
; By D.Nidever  March 2008
;-

COMMON photred,setup

undefine,error

; Not enough inputs
;-------------------
nthisprog = n_elements(thisprog)
if (nthisprog eq 0) then begin
  print,'Syntax - lists = photred_getinput(thisprog,precursor,redo=redo)'
  return,{ninputlines:-1}
endif

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'PHOTRED_GETINPUT ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return,{ninputlines:-1}
endif


; Is "thisprog" a valid stage?
;-------------------------------
thisprog = strupcase(thisprog)
stages = ['RENAME','WCS','SPLIT','DAOPHOT','MATCH','ALLFRAME','APCOR',$
          'CALIB','ASTROM','COMBINE','DEREDDEN','SAVE','HTML',$
          'APERPHOT','DAOGROW','MATCHCAT','COMBINECAT','FITDATA']   ; STDRED stages
; NOT a valid stage
stageind = where(stages eq thisprog,nstageind)
if (nstageind eq 0) then begin
  print,thisprog,' IS NOT A VALID STAGE'
  return,{ninputlines:-1}
endif


; Does the logs/ directory exist?
testlogs = FILE_TEST('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
logfile = 'logs/'+thisprog+'.log'
inputfile = 'logs/'+thisprog+'.inlist'
outputfile = 'logs/'+thisprog+'.outlist'
successfile = 'logs/'+thisprog+'.success'
failurefile = 'logs/'+thisprog+'.failure'
; If the files don't exist create them
if file_test(logfile) eq 0 then TOUCHZERO,logfile
if file_test(inputfile) eq 0 then TOUCHZERO,inputfile
if file_test(outputfile) eq 0 then TOUCHZERO,outputfile
if file_test(successfile) eq 0 then TOUCHZERO,successfile
if file_test(failurefile) eq 0 then TOUCHZERO,failurefile


; Check that all of the required programs are available
progs = ['undefine','readline','readlist','readpar','strsplitter','writeline','push',$
         'printlog','touchzero']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return,{ninputlines:-1}
endif



;############################################
;#  DEALING WITH LIST FILES
;############################################
printlog,logfile,''
printlog,logfile,'--------------------'
printlog,logfile,'CHECKING THE LISTS'
printlog,logfile,'--------------------'
printlog,logfile,''


; CHECK LISTS
;-----------------
READLIST,inputfile,inputlines,/exist,/unique,/fully,count=ninputlines,logfile=logfile,/silent
printlog,logfile,strtrim(ninputlines,2),' files in '+thisprog+'.inlist'

; Load the output list
READLIST,outputfile,outputlines,/unique,/fully,count=noutputlines,logfile=logfile,/silent
printlog,logfile,strtrim(noutputlines,2),' files in '+thisprog+'.outlist'

; Load the success list
READLIST,successfile,successlines,/unique,/fully,count=nsuccesslines,logfile=logfile,/silent
printlog,logfile,strtrim(nsuccesslines,2),' files in '+thisprog+'.success'

; Load the failure list
READLIST,failurefile,failurelines,/unique,/fully,count=nfailurelines,logfile=logfile,/silent
printlog,logfile,strtrim(nfailurelines,2),' files in '+thisprog+'.failure'


; CREATING INLIST
;----------------
nprecursor = n_elements(precursor)
flag = 0
if nprecursor eq 0 then flag=1
count = 0

WHILE (flag eq 0) do begin

  ; Getting the file name
  pre = FILE_BASENAME(precursor[count])
  parts = strsplit(pre,'.',/extract)
  prestage = strupcase(parts[0])
  nparts = n_elements(parts)
  if nparts eq 1 then ending='.outlist' else ending='.'+parts[1]
  prefile = 'logs/'+prestage+ending
  printlog,logfile,'PRESTAGE = ',prestage

  ; Check that this is a valid stage
  prestageind = where(stages eq prestage,nprestageind)
  if (nprestageind eq 0) then begin
    print,prestage,' IS NOT A VALID STAGE'
    goto,BOMB
  endif

  ; Test that the file exists
  pretest = FILE_TEST(prefile)
  if (pretest eq 1) then begin
    ; Read list, make sure the files exist!!
    READLIST,prefile,poutputlines,/exist,/unique,/fully,count=npoutputlines,logfile=logfile,/silent
    printlog,logfile,strtrim(npoutputlines,2),' files in '+prestage+ending+' file that exist'
  endif else npoutputlines=0

  ; Some files in PRECURSOR.outlist
  ; Move/Copy files from PRECURSOR outlist to CURRENT inlist
  if (npoutputlines gt 0) then begin

    ; Write to input file right away
    WRITELINE,inputfile,poutputlines,/append

    ; "Empty" PRECURSOR output
    ; ONLY if it is an "outlist" and not RENAME
    if (prestage ne 'RENAME') and (ending eq '.outlist') and not keyword_set(noempty) then begin
      printlog,logfile,'EMPTYING '+prestage+ending
      FILE_DELETE,prefile
      SPAWN,'touch '+prefile,out
    endif

    ; Add these to the inputlist
    PUSH,inputlines,poutputlines
    ninputlines = n_elements(inputlines)

    ; Remove redundant names
    ui = UNIQ(inputlines,sort(inputlines))
    ui = ui[sort(ui)]
    inputlines = inputlines[ui]
    ninputlines = n_elements(inputlines)

    ; End now
    flag = 1

  ; The PRECURSOR outlist is empty
  endif else begin
    if pretest eq 0 then $
      printlog,logfile,'NO FILES in ',prestage+ending
  endelse


  BOMB:

  ; Have we exhausted the precursor list
  if (count eq (nprecursor-1)) then flag=1

  count++

ENDWHILE



;----------------------------------------------------------------
; DO NOT OVERWRITE/REDO Files already done (unless /REDO is set)
;----------------------------------------------------------------

; Remove all files in the outlist from the inlist, unless REDO is set
;-------------------------------------------------------------------------
if (ninputlines gt 0) and (noutputlines gt 0) then begin
  MATCH,inputlines,outputlines,ind1,ind2,count=nind1

  if (nind1 gt 0) then begin
    printlog,logfile,strtrim(nind1,2),' files in '+thisprog+'.outlist are also in '+thisprog+'.inlist'

    ; REDOING these
    if keyword_set(redo) then begin
      printlog,logfile,'REDO set.  Files in '+thisprog+'.outlist *NOT* removed from '+thisprog+'.inlist'

    ; Not redoing these
    endif else begin
      printlog,logfile,strtrim(nind1,2),' files in '+thisprog+'.outlist removed from '+thisprog+'.inlist'
      if nind1 lt ninputlines then REMOVE,ind1,inputlines
      if nind1 eq ninputlines then undefine,inputlines
      ninputlines = n_elements(inputlines)
    endelse
  endif
endif


; Remove all files in the success list from the inlist, unless REDO is set
;-------------------------------------------------------------------------
if (ninputlines gt 0) and (nsuccesslines gt 0) then begin
  MATCH,inputlines,successlines,ind1b,ind2b,count=nind1b

  if (nind1b gt 0) then begin
    printlog,logfile,strtrim(nind1b,2),' files in '+thisprog+'.success are also in '+thisprog+'.inlist'

    ; REDOING these
    if keyword_set(redo) then begin
      printlog,logfile,'REDO set.  Files in '+thisprog+'.success *NOT* removed from '+thisprog+'.inlist'

    ; Not redoing these
    endif else begin
      printlog,logfile,strtrim(nind1b,2),' files in '+thisprog+'.success removed from '+thisprog+'.inlist'
      if nind1b lt ninputlines then REMOVE,ind1b,inputlines
      if nind1b eq ninputlines then undefine,inputlines
      ninputlines = n_elements(inputlines)
    endelse
  endif
endif

; Checking the EXTENSION
;-------------------------
if (n_elements(extension) gt 0 and ninputlines gt 0) then begin

  extarr = strarr(ninputlines)
  for i=0,ninputlines-1 do extarr[i]=first_el(strsplit(inputlines[i],'.',/extract),/last)
  gdinp = where(extarr eq extension,ngdinp)

  ; Some endings matched
  if (ngdinp gt 0) then begin
    inputlines = inputlines[gdinp]
    ninputlines = ngdinp

  ; None matched
  endif else begin
    undefine,inputlines
    ninputlines = 0
  endelse

  ndiff = ninputlines-ngdinp
  if ndiff gt 0 then printlog,logfile,strtrim(ndiff,2),' DO NOT HAVE THE REQUIRED >>',extension,'<< EXTENSION'

end


; Writing INLIST
;---------------
; Some input files
if (ninputlines gt 0) then begin

  WRITELINE,inputfile,inputlines
  printlog,logfile,strtrim(ninputlines,2),' input files'

; No input files, empty it
endif else begin

  printlog,logfile,'NO input files'
  FILE_DELETE,inputfile
  SPAWN,'touch '+inputfile,out
endelse


; Making LIST structure
;------------------------
;cmd = 'lists = {thisprog:thisprog'
;if nprecursor gt 0 then cmd=cmd+', precursor:precursor'
;if n_elements(prestage) gt 0 then cmd=cmd+', prestage:prestage'
;if n_elements(prefile) gt 0 then cmd=cmd+', prefile:prefile'
;if ninputlines gt 0 then cmd=cmd+', inputlines:inputlines'
;cmd = cmd+', ninputlines:ninputlines'
;if noutputlines gt 0 then cmd=cmd+', outputlines:outputlines'
;cmd = cmd+', noutputlines:noutputlines'
;if nsuccesslines gt 0 then cmd=cmd+', successlines:successlines'
;cmd = cmd+', nsuccesslines:nsuccesslines'
;if nfailurelines gt 0 then cmd=cmd+', failurelines:failurelines'
;cmd = cmd+', nfailurelines:nfailurelines}'
;dum = EXECUTE(cmd)

lists = {thisprog:thisprog}
if nprecursor gt 0 then lists = CREATE_STRUCT(lists,'precursor',precursor)
if n_elements(prestage) gt 0 then lists = CREATE_STRUCT(lists,'prestage',prestage)
if n_elements(prefile) gt 0 then lists = CREATE_STRUCT(lists,'prefile',prefile)
if ninputlines gt 0 then lists = CREATE_STRUCT(lists,'inputlines',inputlines)
lists = CREATE_STRUCT(lists,'ninputlines',ninputlines)
if noutputlines gt 0 then lists = CREATE_STRUCT(lists,'outputlines',outputlines)
lists = CREATE_STRUCT(lists,'noutputlines',noutputlines)
if nsuccesslines gt 0 then lists = CREATE_STRUCT(lists,'successlines',successlines)
lists = CREATE_STRUCT(lists,'nsuccesslines',nsuccesslines)
if nfailurelines gt 0 then lists = CREATE_STRUCT(lists,'failurelines',failurelines)
lists = CREATE_STRUCT(lists,'nfailurelines',nfailurelines)


if keyword_set(stp) then stop

return,lists

end
