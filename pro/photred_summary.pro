pro photred_summary,stp=stp

;+
;
; PHOTRED_SUMMARY
;
; This gives a summary of the progress of the PHOTRED stages
;
; INPUTS:
;  none
;
; OUTPUTS:
;  It prints information to the screen
;
; USAGE:
;  IDL>photred_summary
;
; By D.Nidever   April 2008
;-

;
stages = ['RENAME','SPLIT','WCS','DAOPHOT','MATCH','ALLFRAME',$
          'APCOR','ASTROM','CALIB','COMBINE','DEREDDEN','SAVE','HTML']
nstages = n_elements(stages)

; Does the logs directory exist
test = FILE_TEST('logs',/directory)
if test eq 0 then begin
  print,'NO "logs" directory'
  return
endif

CD,current=curdir
print,''
print,'PHOTRED SUMMARY for ',curdir
print,''
print,'-----------------------------------------------------'
print,'STAGE   INLIST   OUTLIST  SUCCESS  FAILURE  COMPLETED'
print,'-----------------------------------------------------'

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

  ; How many completed, for DAOPHOT and ALLFRAME
  completetext=''
  ; DAOPHOT, Check how many have ".als" and "a.als" files
  if stage eq 'DAOPHOT' and (ninputlines gt 0 or nsuccesslines gt 0) then begin
    undefine,files
    PUSH,files,inputlines
    PUSH,files,successlines
    nfiles = n_elements(files)

    ui = uniq(files,sort(files))
    files = files[ui]
    nfiles = n_elements(files)

    if (nfiles gt 0) then begin
      dirs = FILE_DIRNAME(files)
      base = FILE_BASENAME(files,'.fits')

      completarr = intarr(nfiles)+1

      for j=0,nfiles-1 do begin
        alsfile = dirs[j]+'/'+base[j]+'.als'
        aalsfile = dirs[j]+'/'+base[j]+'a.als'
        ; Check that this file has an ALS file
        alstest = FILE_TEST(alsfile)
        if alstest eq 1 then alslines=FILE_LINES(alsfile) else alslines=0
        ; Check the A.ALS file
        aalstest = FILE_TEST(aalsfile)
        if aalstest eq 1 then aalslines=FILE_LINES(aalsfile) else aalslines=0

        if (alstest eq 0 or alslines lt 3) then completarr[j]=0
        if (aalstest eq 0 or aalslines lt 3) then completarr[j]=0
      end

      ncomplete = long(total(completarr))
      completetext = strtrim(ncomplete,2)
    endif

  endif  ; DAOPHOT complete

  ; ALLFRAME, Check how many have a ".mag" file
  if stage eq 'ALLFRAME' and (ninputlines gt 0 or nsuccesslines gt 0) then begin
    undefine,files
    PUSH,files,inputlines
    PUSH,files,successlines
    nfiles = n_elements(files)

    if (nfiles gt 0) then begin
      dirs = FILE_DIRNAME(files)
      base = FILE_BASENAME(files,'.mch')

      completarr = intarr(nfiles)+1

      for j=0,nfiles-1 do begin
        magfile = dirs[j]+'/'+base[j]+'.mag'
        ; Check that this file has an MAG file
        magtest = FILE_TEST(magfile)
        if magtest eq 1 then maglines=FILE_LINES(magfile) else maglines=0
        if (magtest eq 0 or maglines lt 3) then completarr[j]=0
      end

      ncomplete = long(total(completarr))
      completetext = strtrim(ncomplete,2)
    endif

  endif  ; ALLFRAME complete


  ; Printing
  format = '(A-10,A-9,A-9,A-9,A-9,A-9)'
  print,format=format,stage,intext,outtext,successtext,failuretext,completetext


end

print,'-----------------------------------------------------'

if keyword_set(stp) then stop

end
