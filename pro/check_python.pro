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

pro check_python,test,out,pythonbin=pythonbin,silent=silent,stp=stp

test = 0         ; bad to start with

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'CHECK_PYTHON ERROR: ', !ERROR_STATE.MSG  
   test = -1            ; There was an error
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


; Write a test PYTHON script'
push,cmd,'# -*- coding: utf-8 -*-'
push,cmd,'from __future__ import division'
push,cmd,'import os'
push,cmd,'import glob '
push,cmd,'import numpy as np'
push,cmd,'import sys'
push,cmd,'import time'
push,cmd,'from subprocess import Popen, PIPE '
push,cmd,'from astropy.io.fits import getdata'
push,cmd,'from math import ceil, log10, sqrt '
push,cmd,'from numpy.linalg import inv'
push,cmd,'import random'
push,cmd,''
push,cmd,'print("hello world")'

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
  print,'Running PYTHON from IDL failed!'
  print,'Please check pythonbin in setup file and make sure it is properly set'
  print,'Also check needed modules and packages:'
	print,'os, glob, numpy, sys, time, subprocess, astropy.io.fits, math, random'
  print,'Python returned:'
  print,out
  print,errout
endif


if keyword_set(stp) then stop

end
