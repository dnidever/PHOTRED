;+
;
; PBS_DAEMON
;
; This program is the daemon for photred.  It controls pleione jobs.
; See daophot_setup.pro for how to make the pleione scripts
;
; If we are running on Pleione, multiple jobs is turned on, and multiple
; jobs are input then run the daemon, otherwise just run single jobs
;
; NOTE:  If you want to "kill" all of the jobs create a file "killpbs"
; in the same directory that PBS_DAEMON is being run in and all of the
; PBS jobs will be killed.
;
; INPUTS:
;  input     A string array with the IDL commands (i.e. jobs) to be run.
;  dirs      The directories in which the commands are to be run.
;  /idle     This is an IDL command, otherwise a SHELL command.
;  =prefix   The prefix for the PBS script names
;  =nmulti   How many nodes to run these jobs on.  Default is 8.
;  /hyperthread  Not on a PBS server but one that has multiple processors
;                 hyperthreaded.  Run multiple jobs at the same time on
;                 the same server.
;  =waittime  Time to wait between checking the running jobs.  Default
;                is 60 sec.
;  /htcondor   Submit the jobs using HT Condor.
;  /htc_idlvm  Use the IDL Virtual Machine with HT Condor.
;  =pythonbin  The path to the PYTHON binary.
;  =scriptsdir The directory that contains the IDL .sav files (when
;                using the IDL Virtual Machine) and the HT Condor
;                scripts.
;
; OUTPUTS:
;  Jobs are run.
;  =jobs       The jobs structure with information about each job.
;
; USAGE:
;  IDL>pbs_daemon,input,dirs,jobs=jobs,idle=idle,prefix=prefix,nmulti=nmulti
;
; By D.Nidever   February 2008
;    Antonio Dorta  April-July 2017,  added HT Condor support
;-

pro pbs_daemon,input,dirs,jobs=jobs,idle=idle,prefix=prefix,nmulti=nmulti,  $
               hyperthread=hyperthread,waittime=waittime,cdtodir=cdtodir,   $
               htcondor=htcondor, htc_idlvm=htcondor_idlvm,                 $
               pythonbin=pythonbin, scriptsdir=scriptsdir

; How many input lines
ninput = n_elements(input)
if ninput eq 0 then begin
  print,'Syntax - pbs_daemon,input,dirs,jobs=jobs,idle=idle,prefix=prefix,nmulti=nmulti,'
  print,'                    hyperthread=hyperthread'
  return
endif

if not keyword_set(htcondor) then htcondor='0' $
else if htcondor eq "" then htcondor='0'

; Current directory
CD,current=curdir

ndirs = n_elements(dirs)
if ndirs eq 0 then dirs = replicate(curdir,ninput)
if ndirs eq 1 then dirs = replicate(dirs,ninput)

; Defaults
if n_elements(nmulti) eq 0 then nmulti=8           ; number of jobs to submit at a time
if n_elements(waittime) eq 0 then waittime=60          ; wait time
waittime = waittime > 1

; What host
host = getenv('HOST')
pleione = stregex(host,'pleione',/boolean,/fold_case)
hyades = stregex(host,'hyades',/boolean,/fold_case)

; On Pleione replace /export/local/ with /local/
; the nodes use /local/
if (pleione eq 1) then begin
  bd = where(stregex(dirs,'^/export/local',/boolean) eq 1,nbd)
  if nbd gt 0 then dirs[bd] = strmid(dirs[bd],7,200)
endif

; Which IDL are we using?
if keyword_set(idle) then begin
  SPAWN,'which idl',out,errout
  if STRPOS(out[0],'aliased to') ne -1 then begin
    aliasto = first_el(strsplit(out[0],' ',/extract),/last)
    idlprog = FILE_SEARCH(aliasto,count=nidlprog)
    ; aliased to another program, e.g. gdl
    if nidlprog eq 0 then begin
      spawn,'which '+aliasto,out2,errout2
      idlprog = FILE_SEARCH(out2[0],count=nidlprog)
    endif
  ; no alias
  endif else begin
    idlprog = FILE_SEARCH(out[0],count=nidlprog)
  endelse
  if (nidlprog eq 0) then begin
    print,'IDL PROGRAM NOT AVAILABLE'
    return
  endif
endif

; Create RUNBATCH and IDLBATCH if using /hyperthread
if keyword_set(hyperthread) then begin
  if not keyword_set(idle) and FILE_TEST('runbatch') eq 0 then begin
    undefine,lines
    push,lines,"if test $# -eq 0"
    push,lines,"then"
    push,lines,"  echo 'Syntax - runbatch program'"
    push,lines,"else"
    push,lines,"  echo 'Log file: '$1'.log'"
    push,lines,"  ( nohup  $1 > $1.log 2>&1 ) &"
    push,lines,"  echo $!"
    push,lines,"fi"
    WRITELINE,'runbatch',lines
    FILE_CHMOD,'runbatch','755'o
  endif
  if keyword_set(idle) and FILE_TEST('idlbatch') eq 0 then begin
    undefine,lines
    push,lines,"if test $# -eq 0"
    push,lines,"then"
    push,lines,"  echo 'Syntax - idlbatch idl.batch'"
    push,lines,"else"
    push,lines,"  echo 'Log file: '$1'.log'"
    push,lines,"  ( nohup "+idlprog+" < $1 > $1.log 2>&1 ) &"
    push,lines,"  echo $!"
    push,lines,"fi"
    WRITELINE,'idlbatch',lines
    FILE_CHMOD,'idlbatch','755'o
  endif
endif

print,'---------------------------------'
print,' RUNNING PBS_DAEMON for ',strtrim(ninput,2),' JOB(S)'
print,'---------------------------------'
print,'Host=',host
print,'Nmulti=',strtrim(nmulti,2)

;--------
; DAEMON
;--------
IF (ninput gt 1) and (nmulti gt 1) and ((pleione eq 1) or (hyades eq 1) or (hyperthread eq 1) or (htcondor ne '0')) then begin

  ;-------------
  ; HTCONDOR
  ;-------------
  if htcondor ne '0' then begin

    ;---------------------------------------------------------------
    ; Run CHECK_PYTHON.PRO to make sure that you can run PYTHON from IDL
    ;---------------------------------------------------------------
    CHECK_PYTHON,pythontest,pythonbin=pythonbin
    if pythontest eq 0 then begin
      print,'PYTHON TEST FAILED.  EXITING'
      return
    endif

    print,"Checking HTCondor..."
    CHECK_PYHTCONDOR,pyhtcondortest,pythonbin=pythonbin
    if pyhtcondortest eq 0 then begin
      print,'HTCondor TEST FAILED.  EXITING'
      return
    endif

    if strlowcase(strmid(htcondor,0,7)) eq '#shared' then print,"Using shared directories (file transfer disabled!)"
    ; Generate

    basename = maketemp('htcondor-')+"-"
    if not keyword_set(htcondor_idlvm) then htcondor_idlvm = '0'

    for i=0, ninput-1 do begin

      cmd = input[i]
      cmdfile = basename+strtrim(i,2)+'.sh'

      undefine,lines
      push,lines,'#!/bin/bash'

      commas = strpos(cmd,';')
      if commas[0] ne -1 then begin
        temp = strsplit(cmd,';',/extract)
        cmd = temp
      endif
      if strlowcase(strmid(htcondor,0,7)) eq '#shared' then begin
        push,lines,'cd '+dirs[i]
      endif

      for j=0,n_elements(cmd)-2 do push,lines,cmd[j]
      last_cmd = cmd[n_elements(cmd)-1]
      ; Last value of cmd array should have the command to execute
      ; run it with IDL when idle is set
      if keyword_set(idle) then begin
        last_cmd = STRJOIN(STRSPLIT(last_cmd, "'", /EXTRACT), '"')
        if htcondor_idlvm eq '0' then begin  
          last_cmd = idlprog + " <<< '" + last_cmd + "'"
          htcondor += "|;|concurrency_limits|=|idl:22" ; XXX XXX XXX FIXME!!
        endif else begin
          push,lines, 'ln -s ' + scriptsdir + '/fakered_allframe_wrapper.sav fakered_allframe_wrapper.sav'
          last_cmd = htcondor_idlvm + " fakered_allframe_wrapper.sav '" + last_cmd + "' 0 1"
        endelse
      endif
     
      push,lines,last_cmd
      WRITELINE, cmdfile, lines
      FILE_CHMOD, cmdfile, /U_EXECUTE

    endfor

    executable = basename + "$(Process).sh"

    htcondor_cmd = pythonbin + ' ' + scriptsdir + '/htcondor_run.py '
    ; Running python to manage HTCondor jobs 

    print,systime(0)
    print,"SENDING JOBS TO HTCondor... Please wait!!!"
    undefine,out,errout

    
    SPAWN, htcondor_cmd + 'condor_submit ' + strtrim(ninput,2) + "  '" + executable + "' '" + htcondor  + "'", out, errout

    if n_elements(errout) gt 1 then begin
      print,systime(0)
      print,"There was an error when submitting jobs: " 
      for xx=0,n_elements(errout)-1 do print,"OUTERR: ", errout[xx]
      return
    endif else begin
      clusterId = out
      print,ninput," jobs submitted to HTCondor clusterId ", clusterID
    endelse
   
    count = 0
    flag = 0
    JOBSTATUS = ['Idle', 'Running', 'Removed', 'Completed', 'Held', 'Transferring Output', 'Suspended']
    while (flag eq 0) do begin
      print,systime(0)

      if FILE_TEST('killhtcondor') then begin
        print,"File killhtcondor FOUND... Killing ALL jobs process of clusterId " + clusterId
        print,""
        undefine,out,errout
        SPAWN, htcondor_cmd + 'condor_rm ' + clusterId, out, errout
      endif
 
      undefine,out,errout
      SPAWN, htcondor_cmd + 'condor_q ' + clusterId, out, errout
      if n_elements(errout) gt 1 then begin
      print,"There was an error when checking HTCondor queue. Retrying... " 
      for xx=0,n_elements(errout)-1 do print,"OUTERR: ", errout[xx]
      endif

      if out ne 0 then begin
        condor_queue = strsplit(out,';',/extract)
        condor_info =strtrim(ninput,2) + " total HTCondor jobs. " + condor_queue[0] + " remaining -> "
        for i=1,n_elements(condor_queue)-1 do begin
          condor_status = strsplit(condor_queue[i],":",/extract)
          if condor_status[0] le n_elements(JOBSTATUS) then                                 $
            condor_info += JOBSTATUS[condor_status[0]-1] + ": " + condor_status[1] + "  "    $
          else                                                                               $
            condor_info += "OTHER: " + condor_status[1] + "  "
        endfor
        print,condor_info

        ; Wait 
        ;--------------
        print,'Waiting '+strtrim(waittime,2)+'s'
        print,""
        wait,waittime

      endif else begin
        print,'All jobs are done. Continuing...'
        flag = 1
      endelse


      ; Increment the counter
      count++

    endwhile

  ;-------------
  ; NON-HTCONDOR
  ;-------------
  Endif else begin

    ; Keep submitting jobs until nmulti is reached
    ;
    ; Check every minute or so to see how many jobs are still
    ; running.  If it falls below nmulti and more jobs are left then
    ; submit more jobs
    ;
    ; Don't return until all jobs are done.


    ; Start the "jobs" structure
    ; id will be the ID from Pleione
    dum = {jobid:'',input:'',name:'',scriptname:'',submitted:0,done:0}
    jobs = replicate(dum,ninput)
    jobs.input = input
    njobs = ninput

    ; Loop until all jobs are done
    ; On each loop check the pleione queue and figure out what to do
    count = 0.
    flag = 0
    WHILE (flag eq 0) DO BEGIN

      print,''
      print,systime(0)

      ; Check disk space
      ;-----------------
      ;SPAWN,'df -k '+dirs[0],dfout,dfouterr
      ; Linux
      ;SPAWN,'df -B 1M '+dirs[0],dfout,dfouterr
      ; Mac OS X
      SPAWN,'df -m '+dirs[0],dfout,dfouterr
      dfarr = strsplitter(dfout,' ',/extract)
      nline = n_elements(dfout)
      ;;available = float(reform(dfarr[3,1]))
      ;available = float(reform(dfarr[2,nline-1]))
      available = float(reform(dfarr[3,nline-1]))
      ;available = float(reform(dfarr[2,nline-1]))
      ;available = available[0]/1000.0     ; convert to MB
      print,''
      print,strtrim(available[0],2),' MB of disk space available'
      ; Not enough disk space available
      if available[0] lt 100. then begin
        print,'NOT enough disk space available'
        return
      endif


      ; Check for kill file
      ;--------------------
      kill = FILE_TEST('killpbs')
      if (kill eq 1) then begin
        sub = where(jobs.submitted eq 1 and jobs.done eq 0,nsub)

        print,'Kill file found.  Killing all ',strtrim(nsub,2),' PBS job(s)'

        for i=0,nsub-1 do begin
          ; Killing the job
          print,'Killing ',jobs[sub[i]].name,'  JobID=',jobs[sub[i]].jobid

          if not keyword_set(hyperthread) then begin
            SPAWN,'qdel '+jobs[sub[i]].jobid,out,errout
          endif else begin
            SPAWN,'kill -9 '+jobs[sub[i]].jobid,out,errout
          endelse
        endif

        ; Remove the kill file
        print,'Deleting kill file "killpbs"'
        FILE_DELETE,'killpbs',/allow,/quiet

        ; Now exit
        goto,BOMB

      endif


      ; Check the status of jobs that are running/submitted
      ;----------------------------------------------------
      print,'--Checking queue--'
      sub = where(jobs.submitted eq 1 and jobs.done eq 0,nsub)
      for i=0,nsub-1 do begin

        ; Checking status
        jobid = jobs[sub[i]].jobid
        PBS_CHECKSTAT,statstr,jobid=jobid,hyperthread=hyperthread

        ; Done
        ; Should probably check the output files too
        if statstr.jobid eq '' then begin
          print,'Input ',strtrim(sub[i]+1,2),' ',jobs[sub[i]].name,' JobID=',jobs[sub[i]].jobid,' FINISHED'
          jobs[sub[i]].done=1
        endif

        ; Check for errors as well!! and put in jobs structure
      endfor
      

      ; How many jobs are still in the queue
      ;-------------------------------------
      dum = where(jobs.submitted eq 1 and jobs.done eq 0,Ninqueue)

      ; How many have not been done yet
      ;--------------------------------
      dum = where(jobs.submitted eq 0,Nnosubmit)

      ; What is the current summary
      ;----------------------------------
      dum = where(jobs.done eq 1,nfinished)
      print,'--Jobs Summary--'
      print,strtrim(ninput,2),' total, ',strtrim(nfinished,2),' finished, ',strtrim(Ninqueue,2),$
            ' running, ',strtrim(Nnosubmit,2),' left'


      ; Need to Submit more jobs
      ;-------------------------
      Nnew = (nmulti-ninqueue) > 0
      Nnew = Nnew < Nnosubmit
      If (Nnew gt 0) then begin

        ; Get the indices of new jobs to be submitted
        nosubmit = where(jobs.submitted eq 0)
        newind = nosubmit[0:nnew-1]

        print,''
        print,'--Updating Queue--'
        print,strtrim(ninqueue,2),' JOB(S) running, out of ',strtrim(nmulti,2),' Maximum'
        print,'Submitting ',strtrim(nnew,2),' more job(s)'

        ; Loop through the new submits
        For i=0,nnew-1 do begin

          print,''
          cmd = jobs[newind[i]].input
          if keyword_set(idle) then cmd='IDL>'+cmd
          print,'Input ',strtrim(newind[i]+1,2),'  Command: >>',cmd,'<<'

          ; Make PBS script
          undefine,name,scriptname
          PBS_MAKESCRIPT,jobs[newind[i]].input,dir=dirs[newind[i]],name=name,scriptname=scriptname,$
                         prefix=prefix,idle=idle,hyperthread=hyperthread

          ; Check that the script exists
          test = FILE_TEST(scriptname)

          ; Submitting the job
          if not keyword_set(hyperthread) then begin
            SPAWN,'qsub '+scriptname[0],out,errout
          endif else begin
            if keyword_set(idle) then batchprog='./idlbatch' else batchprog='./runbatch'
            if keyword_set(cdtodir) then cd,dirs[newind[i]]
            SPAWN,batchprog+' '+scriptname[0],out,errout
            if keyword_set(cdtodir) then cd,curdir
          endelse

          ; Check that there weren't any errors
          dum = where(errout ne '',nerror)

          if nerror gt 0 then begin
            print,"ERRORS WHEN EXECUTING!!"
            for xx=0,nerror-1 do print,errout[xx]
            return
          endif

          ; Getting JOBID
          jobid = reform(out)
          if keyword_set(hyperthread) then jobid = reform(out[1])

          ; Printing info
          print,'Submitted ',scriptname[0],'  JobID=',jobid[0]

          ; Updating the jobs structure
          jobs[newind[i]].submitted = 1
          jobs[newind[i]].jobid = jobid[0]
          jobs[newind[i]].name = name[0]
          jobs[newind[i]].scriptname = scriptname[0]

        Endfor  ; submitting new jobs loop
      Endif  ; new jobs to submit


      ; Are we done?
      ;-------------
      dum = where(jobs.done eq 1,ndone)
      if ndone eq njobs then flag=1


      ; Wait a minute
      ;--------------
      print,''
      print,'Waiting '+strtrim(waittime,2)+'s'
      if flag eq 0 then wait,waittime

      ; Increment the counter
      count++

    ENDWHILE  ; large loop until all jobs are done

  ENDELSE  ; NON-HTCONDOR

;------------
; NON-DAEMON
;------------
ENDIF ELSE BEGIN


  ; Loop through the jobs
  FOR i=0,ninput-1 do begin
    
    ; CD to the directory
    ;--------------------
    cd,dirs[i]


    ; Run the IDL command
    if keyword_set(idle) then begin

      ; Execute the command
      ;--------------------
      ; Do we want to do this with idlbatch so that
      ; there is a nice log file
      print,''
      print,'RUNNING IDL>',input[i]
      dum = EXECUTE(input[i])


    ; Run SHELL command
    endif else begin

      ; Execute the command
      print,''
      print,'RUNNING ',input[i]
      SPAWN,input[i]
      ;SPAWN,input[i],out,errout
      ;printline,out
      ;printline,errout

    endelse

    cd,curdir  ; back to original directory

  ENDFOR

ENDELSE


BOMB:

;stop

if keyword_set(stp) then stop

end
