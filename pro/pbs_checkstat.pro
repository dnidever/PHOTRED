;+
;
; PBS_CHECKSTAT
;
; This checks the status of PBS jobs
; If no jobs are found in queue then an empty
; statstr structure is returned.
;
; INPUTS:
;  =jobid   Specific JOBID to check.
;  /stp     Stop at the end of the program.
;  /hyperthread  Not on a PBS machine but one with multipe hyperthreaded
;                  processors running simultaneously.
;
; OUTPUTS:
;  statstr  Stat structure.
;
; USAGE:
;  IDL>pbs_checkstat,statstr,jobid=jobid,stp=stp
;
; By D.Nidever   February 2008
;-

pro pbs_checkstat,statstr,jobid=jobid,stp=stp,hyperthread=hyperthread

njobid = n_elements(jobid)

; PBS
;--------
If not keyword_set(hyperthread) then begin

  addon = ''
  if njobid gt 0 then addon=' '+jobid[0]

  SPAWN,'qstat'+addon,out,errout
  dum = where(out ne '',nout)
  ;if nout gt 2 and strmid(out[0],0,3) eq 'Job' then statlines=out[2:*]
  gd = where(stregex(out,'^'+jobid,/boolean) eq 1,ngd)
  if ngd gt 0 then begin
    statlines = reform(out[gd[0]])
  endif else begin
    undefine,statlines
  endelse

  nstat = n_elements(statlines)
  ; Some jobs in queue
  if nstat gt 0 then begin
    arr = strsplitter(statlines,' ',/extract)
    dumdum = {jobid:'',name:'',user:'',timeuse:'',status:'',queue:''}
    statstr = replicate(dumdum,nstat)
    statstr.jobid = reform(arr[0])
    statstr.name = reform(arr[1])
    statstr.user = reform(arr[2])
    statstr.timeuse = reform(arr[3])
    statstr.status = reform(arr[4])
    statstr.queue = reform(arr[5])
  ; No jobs in queue
  endif else begin
    statstr = {jobid:'',name:'',user:'',timeuse:'',status:'',queue:''}
  endelse

; Hyperthreaded.  Need a jobid
;-----------------------------
Endif else begin

  ; No JOBID input
  if njobid eq 0 then begin
    print,'Need JOBID with /hyperthread'
    statstr = {jobid:'',name:'',user:'',timeuse:'',status:'',queue:''}
    return
  endif

  ;SPAWN,['ps','-p',strtrim(jobid[0],2)],out,errout,/noshell
  SPAWN,['ps','-o','pid,user,etime,command','-p',strtrim(jobid[0],2)],out,errout,/noshell
  if n_elements(out) eq 0 then $
    SPAWN,['ps','-o','pid,user,etime,command','-p',strtrim(jobid[0],2)],out,errout,/noshell

  ; can put in the column that you want
  ; ps -o etime -p jobid

  ; Nothing returned
  if n_elements(out) eq 0 then begin
    statstr = {jobid:strtrim(jobid[0],2),name:'?',user:'?',timeuse:'?',status:'?',queue:'hyperthread'}
    return
  endif

  out = strtrim(out,2)
  gd = where(stregex(out,'^'+strtrim(jobid,2),/boolean) eq 1,ngd)
  if ngd gt 0 then begin
    statlines = reform(out[gd[0]])
  endif else begin
    undefine,statlines
  endelse

  nstat = n_elements(statlines)
  ; Some jobs in queue
  if nstat gt 0 then begin
    arr = strsplitter(statlines,' ',/extract)
    dumdum = {jobid:'',name:'',user:'',timeuse:'',status:'',queue:''}
    statstr = replicate(dumdum,nstat)
    statstr.jobid = reform(arr[0])
    statstr.user = reform(arr[1])
    ; CAN'T get the name.
    ;statstr.name = reform(arr[1])
    statstr.timeuse = reform(arr[2])
    statstr.status = 'R'
    statstr.queue = 'hyperthread'
  ; No jobs in queue
  endif else begin
    statstr = {jobid:'',name:'',user:'',timeuse:'',status:'',queue:''}
  endelse


Endelse

;stop

if keyword_set(stp) then stop

end
