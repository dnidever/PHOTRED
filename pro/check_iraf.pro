pro check_iraf,test,out,irafdir=irafdir,silent=silent,stp=stp

;+
;
; CHECK_IRAF
;
; This program checks that you can run IRAF from IDL.
;
; INPUTS:
;  =irafdir   The path to the IRAF directory
;              where login.cl should be located.  If
;              this is not entered then it assumes that
;              ~/iraf/ is the IRAF directory.
;  /silent    Don't print anything to the screen.
;  /stp       Stop at the end of the program.
;
; OUTPUTS:
;  test       1 if the test was successful, 0 if it was not, and -1 if there was an error.
;  out        The IRAF output from the test
;
;
; By D.Nidever  Nov.2007
;-

test = 0         ; bad to start with

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'CHECK_IRAF ERROR: ', !ERROR_STATE.MSG  
   test = -1            ; There was an error
   CATCH, /CANCEL 
   return
endif



; Check that we have all the required programs
progs = ['maketemp','writeline','undefine']
progtest = PROG_TEST(progs)
if min(progtest) eq 0 then return

; Input IRAF directory
if n_elements(irafdir) gt 0 then begin

  diriraf = file_search(irafdir,/fully_qualify_path)

  dirtest = file_test(diriraf,/directory)
  if dirtest eq 0 then begin
    if not keyword_set(silent) then $
      print,'DIRECTORY ',irafdir,' DOES NOT EXIST.  TRYING ~/iraf/ INSTEAD.'
    undefine,diriraf
  endif

endif

; No IRAF directory yet
if n_elements(diriraf) eq 0 then begin

  diriraf = file_search('~/iraf/',/fully_qualify_path)
  dirtest = file_test(diriraf,/directory)
  if dirtest eq 0 then begin
    if not keyword_set(silent) then $
      print,'NO IRAF DIERECTORY. RETURNING'
    return
  endif

end


CD,current=curdir


; Write a test IRAF script
push,cmd,'print("hello world")'
push,cmd,'logout'
cmdfile = maketemp('temp','.cl')
WRITELINE,cmdfile,cmd

; Goto the IRAF directory
CD,diriraf

; Running IRAF
undefine,out
SPAWN,'cl < '+curdir+'/'+cmdfile,out,errout

; Return to original directory
CD,curdir

; Erasing the temporary files
FILE_DELETE,cmdfile,/allow,/quiet


; See if it printed the message
gd = where(stregex(out,'hello world',/boolean) eq 1,ngd)
if ngd gt 0 then test = 1          ; YES, it worked

; Explain how to fix the login.cl file
if test eq 0 and not keyword_set(silent) then begin
  print,'Running IRAF from IDL failed!'
  print,'Please edit your login.cl file so that it does not'
  print,'print anything to the screen on IRAF startup.'
  print,'The most likely cause are the 9 lines after'
  print,'"Set the terminal type".'
  print,'On PLEIONE the 4 lines after'
  print,'"# Delete any old MTIO lock (magtape position) files."'
  print,'can be problematic.'
end


if keyword_set(stp) then stop

end
