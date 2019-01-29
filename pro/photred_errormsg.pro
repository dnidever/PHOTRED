;+
;
; PHOTRED_ERRORMSG
;
; Print out the error message and traceback
;
; INPUTS:
;  =logfile  The name of a logfile to print to using printlog.pro.
;
; OUTPUTS:
;  Error messages and traceback information printed to the screen.
;
; USAGE:
; IDL>photred_errormsg
;
; By D.Nidever  Jan 2019, based on suggestion from David Fanning here
;   http://www.idlcoyote.com/widget_tips/widgettrace.html and his cgErrorMsg.pro program.
;-

pro photred_errormsg,logfile=logfile

;; This prints out an error message and tracekback information

if n_elements(logfile) eq 0 then logfile=-1

help, Calls=callStack
callingRoutine = (Str_Sep(StrCompress(callStack[1])," "))[0]
help, /Last_Message, Output=traceback
printlog,logfile,''
printlog,logfile, 'Traceback Report from ' + StrUpCase(callingRoutine) + ':'
printlog,logfile, ''
for j=0,n_elements(traceback)-1 do printlog,logfile, "     " + traceback[j]

end
