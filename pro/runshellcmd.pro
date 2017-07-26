;+
;
; RUNSHELLCMD
;
; Function to run a shell command using ILD SPAWN
; A default value can be set.
;
; INPUTS:
;  cmd        Shell command to be executed
;  msg        If set, message to be printed when there are errors
;  logfile    Log file to write message to (stdout will be used if no logfile)
;  /printcmd  Print command     
;  /printout  Print output (stdout)
;  /quiet     Silent run, NO messages will be printed (even if errors)
;
; OUTPUTS:
;  exitcode     exit status after running command (returned)
;  output       messages printed on stdout when running command
;  errors       messages printed on stderr when running command
;
; USAGE:
;  IDL> exitcode = runshellcmd('ls -l','Error when getting list',logfile=logfile,output=out,errors=err)
;
; By Antonio Dorta <adorta@iac.es> July 2017
;-

function runshellcmd, cmd, msg=msg, logfile=logfile, output=output, errors=errors, quiet=quiet, printcmd=printcmd,printout=printout

  ; Run command
  if not keyword_set(quiet) and keyword_set(printcmd) then                       $
    if keyword_set(logfile) then printlog,logfile,"Running shell command: ",cmd else print,"Running shell command",cmd
  SPAWN, cmd, output, errors, EXIT_STATUS=exitcode

  ; Check if execution was successful 
  if exitcode ne 0 then begin

    ; There were errors (exitcode != 0). Print them if not quiet!
    if not keyword_set(quiet) then begin
      if keyword_set(logfile) then begin
        printlog,logfile,'' 
        printlog,logfile,"#############"
        printlog,logfile,"# ERROR!!   #"
        printlog,logfile,"#############"
        printlog,logfile,"There was an error when executing command: ", cmd 
        printlog,logfile,"Error was:", errors
        if keyword_set(msg) then printlog,logfile,msg
        printlog,logfile,'' 
      endif else begin
        print,'' 
        print,"#############"
        print,"# ERROR!!   #"
        print,"#############"
        print,"There was an error when executing command: ", cmd 
        print,"Error was:", errors
        if keyword_set(msg) then print,msg
        print,'' 
      endelse
    endif
  endif
  
  if not keyword_set(quiet) and keyword_set(printout) then $
    if keyword_set(logfile) then printlog,logfile,output else print,output

  return,exitcode
end
