;+
;
; FAKERED_SUMMARY
;
; This gives a summary of the progress of the FAKERED stages
;
; INPUTS:
;  none
;
; OUTPUTS:
;  It prints information to the screen
;
; USAGE:
;  IDL>fakered_summary
;
; By D.Nidever   Feb 2017, mostly copy of stdred_summary.pro
;-

pro fakered_summary,stp=stp

;
stages = ['COPY','ADDSTAR','DAOPHOT','MATCH','ALLFRAME','ASTROM',$
          'CALIB','COMBINE','COMPLETE']  
nstages = n_elements(stages)

; Does the logs directory exist
test = FILE_TEST('logs',/directory)
if test eq 0 then begin
  print,'NO "logs" directory'
  return
endif

CD,current=curdir
print,''
print,'FAKERED SUMMARY for ',curdir
print,''
print,'--------------------------------------------'
print,'STAGE    INLIST   OUTLIST  SUCCESS  FAILURE '
print,'--------------------------------------------'

; Loop through the stages
for i=0,nstages-1 do begin

  stage = stages[i]

  undefine,inputlines,outputlines,successlines,failurelines

  ninputlines=0
  noutputlines=0
  nsuccesslines=0
  nfailurelines=0

  inputfile = 'logs/'+stage+'.inlist'
  outputfile = 'logs/'+stage+'.outlist'
  successfile = 'logs/'+stage+'.success'
  failurefile = 'logs/'+stage+'.failure'

  ; INLIST
  intest = FILE_TEST(inputfile)
  if intest eq 1 then $
  READLIST,inputfile,inputlines,/unique,/fully,count=ninputlines,/silent
  if intest eq 0 then intext='--' else intext = strtrim(ninputlines,2)

  ; OUTLIST
  outtest = FILE_TEST(outputfile)
  if outtest eq 1 then $
  READLIST,outputfile,outputlines,/unique,count=noutputlines,/silent
  if outtest eq 0 then outtext='--' else outtext = strtrim(noutputlines,2)

  ; SUCCESS LIST
  successtest = FILE_TEST(successfile)
  if successtest eq 1 then $
  READLIST,successfile,successlines,/unique,count=nsuccesslines,/silent
  if successtest eq 0 then successtext='--' else successtext = strtrim(nsuccesslines,2)

  ; FAILURE LIST
  failuretest = FILE_TEST(failurefile)
  if failuretest eq 1 then $
  READLIST,failurefile,failurelines,/unique,count=nfailurelines,/silent
  if failuretest eq 0 then failuretext='--' else failuretext = strtrim(nfailurelines,2)


  ; Printing
  format = '(A-11,A-9,A-9,A-9,A-9,A-9)'
  print,format=format,stage,intext,outtext,successtext,failuretext

endfor

print,'--------------------------------------------'

if keyword_set(stp) then stop

end
