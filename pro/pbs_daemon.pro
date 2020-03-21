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
;  input        A string array with the IDL commands (i.e. jobs) to be run.
;  dirs         The directories in which the commands are to be run.
;  =lockfiles   List of "lock" files to use indicate that a job is being
;                 worked on either by this server or another one.
;  =donefiles   List of "done" files to indicate a file is finished and
;                 should not be redone.
;  /idle        This is an IDL command, otherwise a SHELL command.
;  =prefix      The prefix for the PBS script names
;  =nmulti      How many nodes to run these jobs on.  Default is 8.
;  /hyperthread Not on a PBS server but one that has multiple processors
;                 hyperthreaded.  Run multiple jobs at the same time on
;                 the same server.
;  =scriptsdir  The directory that contains the PHOTRED scripts.
;  =waittime    Time to wait between checking the running jobs.  Default
;                 is 5 sec.
;  =verbose     Currently this only controls if the command is printed.
;
; OUTPUTS:
;  Jobs are run.
;
; USAGE:
;  IDL>pbs_daemon,input,dirs,jobs=jobs,idle=idle,prefix=prefix,nmulti=nmulti
;
; By D.Nidever   February 2008
;-

pro pbs_daemon,input,dirs,jobs=jobs,idle=idle,prefix=prefix,nmulti=nmulti,  $
               hyperthread=hyperthread,waittime=waittime,cdtodir=cdtodir,   $
               htcondor=htcondor, htc_idlvm=htcondor_idlvm,                 $
               pythonbin=pythonbin, scriptsdir=scriptsdir,verbose=verbose,  $
               lockfiles=lockfiles,donefiles=donefiles

t0 = systime(1)
  
; How many input lines
ninput = n_elements(input)
if ninput eq 0 then begin
  print,'Syntax - pbs_daemon,input,dirs,jobs=jobs,idle=idle,prefix=prefix,nmulti=nmulti,'
  print,'                    hyperthread=hyperthread'
  return
endif

; Lockfiles input but not donefiles
if n_elements(lockfiles) gt 0 and n_elements(donefiles) eq 0 then begin
  print,'Must use DONEFILES with LOCKFILES'
  return
endif

if not keyword_set(htcondor) then htcondor='0' $
else if htcondor eq "" then htcondor='0'
if n_elements(verbose) eq 0 then verbose=1

; Current directory
CD,current=curdir

if not keyword_set(scriptsdir) then scriptsdir=curdir  ; use absolute path

ndirs = n_elements(dirs)
if ndirs eq 0 then dirs = replicate(curdir,ninput)
if ndirs eq 1 then dirs = replicate(dirs,ninput)

; Defaults
if n_elements(nmulti) eq 0 then nmulti=8              ; number of jobs to submit at a time
if n_elements(waittime) eq 0 then waittime=5          ; wait time
waittime = waittime > 0.01

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
  ; Try GDL
  if (nidlprog eq 0) then begin
    SPAWN,'which gdl',out,errout
    idlprog = FILE_SEARCH(out[0],count=nidlprog)
  endif
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
  if keyword_set(idle) then begin
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


IF (nmulti gt 1) and ((pleione eq 1) or (hyades eq 1) or (hyperthread eq 1) or (htcondor ne '0')) then begin


  if htcondor ne '0' then begin
    ;-------------
    ; HTCONDOR
    ;-------------

    ;---------------------------------------------------------------
    ; Run CHECK_PYTHON.PRO and CHECK_PYHTCONDOR to make sure that you can run PYTHON/HTCondor from IDL
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

    ; HTCondor files will begin with prefix htcondor-...
    basename = maketemp('htcondor-')+"-"

    ; Check if the htcondor command begins with "shared" to use shared directories 
    ;---------------------------------------------------------------
    if strcmp(htcondor, '#shared', 7, /FOLD_CASE) then print,"Using shared directories (file transfer disabled!)"
    ; Generate

    ; Check if we are going to use the IDL Virtual Machine
    if not keyword_set(htcondor_idlvm) then htcondor_idlvm = '0'

    ;---------------------------------------------------------------
    ; Create each job from command list (cmd)
    ;---------------------------------------------------------------
    for i=0, ninput-1 do begin

      cmd = input[i]
      cmdfile = basename+strtrim(i,2)+'.sh'

      ; Write bash script file with all commands
      undefine,lines
      push,lines,'#!/bin/bash'
      ; If there are several commands per line (separated by ";"), split them in different lines
      commas = strpos(cmd,';')
      if commas[0] ne -1 then begin
        temp = strsplit(cmd,';',/extract)
        cmd = temp
      endif
 
      ; Check if we are going to transfer files or use shared directory
      if strcmp(htcondor, '#shared', 7, /FOLD_CASE) then $
        push,lines,'cd '+dirs[i]

      ; Print each command
      for j=0,n_elements(cmd)-2 do push,lines,cmd[j]
      last_cmd = cmd[n_elements(cmd)-1]
      ; Last value of cmd array should have the command to execute
      ; run it with IDL when idle is set
      if keyword_set(idle) then begin
        last_cmd = STRJOIN(STRSPLIT(last_cmd, "'", /EXTRACT), '"')

        ; Are we going to use IDL VM??
        if htcondor_idlvm eq '0' then begin  
          ; No Virutal Machine, use IDL
          last_cmd = idlprog + " <<< '" + last_cmd + "'"
        endif else begin
          ; Use the IDL Virtual Machine
          push,lines, 'ln -s ' + scriptsdir + '/fakered_allframe_wrapper.sav fakered_allframe_wrapper.sav'
          last_cmd = htcondor_idlvm + " fakered_allframe_wrapper.sav '" + last_cmd + "' 0 1"
        endelse
      endif
    
      ; Write command file and make it executable 
      push,lines,last_cmd
      WRITELINE, cmdfile, lines
      FILE_CHMOD, cmdfile, /U_EXECUTE

    endfor

    ; Get executable (script filename)
    executable = basename + "$(Process).sh"

    ; Get HTCondor script that manages the submissions and queue
    ; Running python to manage HTCondor jobs 
    htcondor_cmd = pythonbin + ' ' + scriptsdir + '/htcondor_run.py '

    ;---------------------------------------------------------------
    ; SUBMIT JOBS!! 
    ;---------------------------------------------------------------
    print,systime(0)
    print,"SENDING JOBS TO HTCondor... Please wait!!!"
    undefine,out,errout

    SPAWN, htcondor_cmd + 'condor_submit ' + strtrim(ninput,2) + "  '" + executable + "' '" + htcondor  + "'", out, errout

    ; -----------------------------------
    ; Check submission
    ; -----------------------------------
    if out eq '0' then begin   ; ClusterID is 0, there was an error
      print,''
      print,systime(0)
      print,"THERE WAS AN ERROR WHEN SUBMITTING JOBS."
      print,"If you have used htcondor_cmd in fakered.setup to specify"
      print,"HTCondor submit commands, please, check they are correct." 
      print,"Also check the syntax is valid:"
      print," - use |=| to assign values to commands"
      print," - use |;| to separate commands"
      print," - and do NOT use extra blanks or other separators."
      print,"EXAMPLE (fakered.setup):"
      print,"htcondor         1"
      print,"htcondor_cmd     request_memory|=|5GB|;|request_disk|=|200GB"
      print,''
      print,'You can still write .c to continue the execution'
      stop
    endif
    if n_elements(errout) gt 1 then begin   ; There were errors when submitting
      print,''
      print,systime(0)
      print,"THERE WAS AN ERROR WHEN SUBMITTING JOBS."
      for xx=0,n_elements(errout)-1 do print,"OUTERR: ", errout[xx]
      print,''
      print,'You can still write .c to continue the execution'
      stop
    endif

    clusterId = out
    print,ninput," jobs submitted to HTCondor clusterId ", clusterID
  

    ;---------------------------------------------------------------
    ; Check queue, wait until completion 
    ;---------------------------------------------------------------
    count = 0
    flag = 0
    JOBSTATUS = ['Idle', 'Running', 'Removed', 'Completed', 'Held', 'Transferring Output', 'Suspended']
    while (flag eq 0) do begin
      print,systime(0)

      ; Check if file 'killhtcondor' exists. If so, REMIOVE all jobs of this cluster
      if FILE_TEST('killhtcondor') then begin
        ; If file killhtcondor is found, remove all jobs
        print,"File killhtcondor FOUND... Killing ALL jobs process of clusterId " + clusterId
        print,""
        undefine,out,errout
        SPAWN, htcondor_cmd + 'condor_rm ' + clusterId, out, errout
      endif

      ; Get info from the queue 
      undefine,out,errout
      SPAWN, htcondor_cmd + 'condor_q ' + clusterId, out, errout
      if n_elements(errout) gt 1 then begin
      print,"There was an error when checking HTCondor queue. Retrying... " 
      for xx=0,n_elements(errout)-1 do print,"OUTERR: ", errout[xx]
      endif

      ; There are still some pending jobs. Print info from queue
      if out ne 0 then begin
        condor_queue = strsplit(out,';',/extract)
        condor_info =strtrim(ninput,2) + " total HTCondor jobs. " + condor_queue[0] + " remaining -> "
        for i=1,n_elements(condor_queue)-1 do begin
          condor_status = strsplit(condor_queue[i],":",/extract)
          if condor_status[0] le n_elements(JOBSTATUS) then                                  $
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


  endif else begin
  ;-------------
  ; NON-HTCONDOR
  ;-------------


  ; Keep submitting jobs until nmulti is reached
  ;
  ; Check every minute or so to see how many jobs are still
  ; running.  If it falls below nmulti and more jobs are left then
  ; submit more jobs
  ;
  ; Don't return until all jobs are done.


  ; Start the "jobs" structure
  ; id will be the ID from Pleione
  schema = {jobid:'',cmd:'',dir:'',name:'',scriptname:'',submitted:0,done:0}
  if n_elements(lockfiles) gt 0 then schema=create_struct(schema,'lockfile','','donefile','','locked',0)
  jobs = replicate(schema,ninput)
  jobs.cmd = input
  jobs.dir = dirs
  if n_elements(lockfiles) gt 0 then begin
    jobs.lockfile = lockfiles
    jobs.donefile = donefiles
  endif
  njobs = ninput

  ; Loop until all jobs are done
  ; On each loop check the pleione queue and figure out what to do
  count = 0LL
  endflag = 0
  WHILE (endflag eq 0) DO BEGIN

    print,''
    print,systime(0)

    ; Check disk space
    ;-----------------
    ;SPAWN,'df -k '+dirs[0],dfout,dfouterr
    ; Linux
    ;SPAWN,'df -B 1M '+dirs[0],dfout,dfouterr
    ; Mac OS X
    available = -1.0
    SPAWN,'df -m '+dirs[0],dfout,dfouterr
    if n_elements(dfout) gt 1 and dfouterr[0] eq '' then begin
      ;; sometimes the information is split across two lines, e.g.
      ;; Filesystem           1M-blocks      Used Available Use% Mounted on
      ;;  dlnfs.datalab.noao.edu:/net/dl1
      ;;                      366257664 333705629  32552035  92% /net/dl1
      if n_elements(dfout) eq 3 then line=dfout[1]+'  '+dfout[2] else line = first_el(dfout,/last) 
      dfarr = strsplit(line,' ',/extract)
      if n_elements(dfarr) ge 4 then available = float(dfarr[3])
    endif
    if available gt 0 then begin
      print,''
      print,strtrim(available[0],2),' MB of disk space available'
      ; Not enough disk space available
      if available[0] lt 100. then begin
        print,'NOT enough disk space available'
        return
      endif
    endif else print,'Cannot determine available disk space'


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
      end

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
        ; Create DONE file
        if n_elements(donefiles) gt 0 then begin
          if file_test(jobs[sub[i]].donefile) eq 0 then begin
             doneline = 'Finished '+systime(0)+'  PID='+jobs[sub[i]].jobid+'  HOST='+host+'  '+jobs[sub[i]].scriptname
             WRITELINE,jobs[sub[i]].donefile,doneline            
          endif
        endif
        ; Remove lockfile
        if n_elements(lockfiles) gt 0 then $
          if jobs[sub[i]].lockfile ne '' then FILE_DELETE,jobs[sub[i]].lockfile,/allow,/quiet
      endif


      ; Check for errors as well!! and put in jobs structure

    endfor ; check jobs loop


    ; How many jobs are still in the queue
    ;-------------------------------------
    dum = where(jobs.submitted eq 1 and jobs.done eq 0,Ninqueue)

    ; How many have not been done yet
    ;--------------------------------
    dum = where(jobs.submitted eq 0 and jobs.done eq 0,Nnosubmit)

    ; What is the current summary
    ;----------------------------------
    dum = where(jobs.done eq 1,ndone)       ; we finished these
    dum = where(jobs.done gt 0,nfinished)   ; finished by all processes
    print,'--Jobs Summary--'
    print,strtrim(ninput,2),' total, ',strtrim(ndone,2),'/',strtrim(nfinished,2),' done/finished, ',strtrim(Ninqueue,2),$
          ' running, ',strtrim(Nnosubmit,2),' left'


    ; Need to Submit more jobs
    ;-------------------------
    Nnew = 0 > (nmulti-ninqueue) < Nnosubmit
    Nnew = Nnew < Nnosubmit
    If (Nnew gt 0) then begin

      ; Get the indices of new jobs to be submitted
      nosubmit = where(jobs.submitted eq 0 and jobs.done eq 0,n_nosubmit)

      print,''
      print,'--Updating Queue--'
      print,strtrim(ninqueue,2),' JOB(S) running, out of ',strtrim(nmulti,2),' Maximum'
      print,'Submitting ',strtrim(nnew,2),' more job(s)'

      ; Loop through the new submits
      submitflag = 0
      i = 0L
      newsubmitted = 0L
      WHILE (submitflag eq 0) do begin

        ; Checking DONE and LOCK files
        if n_elements(lockfiles) gt 0 then begin
           ; Check if the DONEFILE exists
           test = file_test(jobs[nosubmit[i]].donefile)
           if test eq 1 then begin
              jobs[nosubmit[i]].done = 2
              goto,BOMB_SUBMIT
           endif
           ; Check the lockfile
           lock = LOCKFILE(jobs[nosubmit[i]].lockfile)
           if lock eq 0 then begin
             jobs[nosubmit[i]].locked = 1
             goto,BOMB_SUBMIT
           endif
        endif
         
        print,''
        cmd = jobs[nosubmit[i]].cmd
        if keyword_set(idle) then cmd='IDL>'+cmd
        if verbose gt 0 then print,'Input ',strtrim(nosubmit[i]+1,2),'  Command: >>',cmd,'<<'

        ; Make PBS script
        undefine,name,scriptname
        PBS_MAKESCRIPT,jobs[nosubmit[i]].cmd,dir=jobs[nosubmit[i]].dir,name=name,scriptname=scriptname,$
                       prefix=prefix,idle=idle,hyperthread=hyperthread

        ; Check that the script exists
        test = FILE_TEST(scriptname)

        ; Submitting the job
        if not keyword_set(hyperthread) then begin
          SPAWN,'qsub '+scriptname[0],out,errout
        endif else begin
          if keyword_set(idle) then batchprog=scriptsdir+'/idlbatch' else batchprog=scriptsdir+'/runbatch'
          if keyword_set(cdtodir) then cd,jobs[nosubmit[i]].dir
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
        jobs[nosubmit[i]].submitted = 1
        jobs[nosubmit[i]].jobid = jobid[0]
        jobs[nosubmit[i]].name = name[0]
        jobs[nosubmit[i]].scriptname = scriptname[0]

        newsubmitted++  ; increment
        
        BOMB_SUBMIT:

        ; Are we done?
        if (newsubmitted eq nnew) or (i eq n_nosubmit-1) then submitflag=1
        i++                     ; increment 
      Endwhile  ; submitting new jobs loop
      if newsubmitted eq 0 then print,'No jobs to submit'
    Endif   ; new jobs to submit


    ; Are we done?
    ;-------------
    dum = where(jobs.done gt 0,ndone)
    if n_elements(lockfiles) gt 0 then dum = where(jobs.done gt 0 or jobs.locked eq 1,ndone)
    if ndone eq njobs then endflag=1


    ; Wait a minute
    ;--------------
    if endflag eq 0 then begin
      print,''
      print,'Waiting '+strtrim(waittime,2)+'s'
      wait,waittime
    endif
      
    ; Increment the counter
    count++

    ;stop

  ENDWHILE

;-------------
; NON-HTCONDOR
;-------------
 ENDELSE

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

print,''
print,'dt = ',strtrim(systime(1)-t0,2),' sec'

BOMB:

;stop

if keyword_set(stp) then stop

end

