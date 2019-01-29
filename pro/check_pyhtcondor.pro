;+
;
; CHECK_PYTHON
;
; This program checks that you can run PYTHON from IDL.
;
; INPUTS:
;  =pythonbin   The path to the PYTHON binary
;  /silent    Don't print anything to the screen.
;  /stp       Stop at the end of the program.
;
; OUTPUTS:
;  test       1 if the test was successful, 0 if it was not, and -1 if there was an error.
;  out        The PYTHON output from the test
;
; USAGE:
;  IDL>check_python,pythontest,pythonbin=pythonbin
;
; By Antonio Dorta <adorta@iac.es> based on check_iraf by D.Nidever  Apr.2017
;-

pro check_pyhtcondor,test,out,pythonbin=pythonbin,silent=silent,stp=stp

test = 0         ; bad to start with

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   test = -1            ; There was an error
   PHOTRED_ERRORMSG
   CATCH, /CANCEL 
   return
endif



; Check that we have all the required programs
progs = ['maketemp','writeline','undefine']
progtest = PROG_TEST(progs)
if min(progtest) eq 0 then return

; Input PYTHON directory
if n_elements(pythonbin) gt 0 and pythonbin ne 'python' then begin

  filetest = file_test(pythonbin,/executable)
  if filetest eq 0 then begin
    if not keyword_set(silent) then $
      print,'EXECUTABLE ',pythonbin,' DOES NOT EXIST.  TRYING python INSTEAD.'
      pythonbin='python'
  endif

endif

CD,current=curdir


; Write a test HTCondor PYTHON script'
push,cmd,'# -*- coding: utf-8 -*-'
push,cmd,'########################################################'
push,cmd,'#'
push,cmd,'# PURPOSE: Simple HTCondor Python test.'
push,cmd,'#          It checks Python bindings are installed'
push,cmd,'#          And whether it is possible to submit one job'
push,cmd,'#'
push,cmd,'#  AUTHOR: Antonio Dorta <adorta@iac.es>'
push,cmd,'#    DATE: 2017-06-20'
push,cmd,'#'
push,cmd,'########################################################'
push,cmd,''
push,cmd,'# Import HTCondor Python'
push,cmd,'import htcondor'
push,cmd,'import classad'
push,cmd,''
push,cmd,'N = 1 '
push,cmd,'schedd = htcondor.Schedd()'
push,cmd,'sub = htcondor.Submit({"executable": "/bin/hostname"})'
push,cmd,'clusterId = 0'
push,cmd,'# Submit job'
push,cmd,'with schedd.transaction() as txn:'
push,cmd,'  clusterId = sub.queue(txn, N)'
push,cmd,''
push,cmd,'if clusterId:'
push,cmd,'  # Jobs submitted!! Remove it!'
push,cmd,'  schedd.act(htcondor.JobAction.Remove, "ClusterId == %d" % clusterId)'
push,cmd,'  print ("hello world")'
push,cmd,'else:'
push,cmd,'  print ("ERROR: No jobs were submitted")'
push,cmd,''
cmdfile = maketemp('temp','.py')
WRITELINE,cmdfile,cmd

; Running PYTHON
undefine,out
SPAWN,pythonbin+' '+curdir+'/'+cmdfile,out,errout

; Return to original directory
CD,curdir

; Erasing the temporary files
FILE_DELETE,cmdfile,/allow,/quiet


; See if it printed the message
gd = where(stregex(out,'hello world',/boolean) eq 1,ngd)
if ngd gt 0 then test = 1          ; YES, it worked

; Explain how to fix the login.cl file
if test eq 0 and not keyword_set(silent) then begin
  print,''
  print,'ERROR: Running a simple HTCondor job failed!'
  print,'Python bindings for HTCondor are used to manage HTCondor jobs.'
  print,'Please, check pythonbin in setup file and make sure it is properly set.'
  print,'Also check needed python modules and packages: import htcondor, classad'
  print,'This script must be executed on a machine with HTCondor installed and'
  print,'with a running Schedd where jobs can be submitted.'
  print,'Python returned:'
  print,out
  print,errout
endif


if keyword_set(stp) then stop

end
