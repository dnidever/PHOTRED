;+
; 
; IRAF_RUN
;
; This runs an IRAF script.
;
; INPUTS:
;  scriptname   The absolute filename of the IRAF script
;               It must end with "logout".
;  irafdir      The IRAF home directory.
;  /silent      Don't print anything to the screen
;
; OUTPUTS:
;  The IRAF script will be executed and the output printed
;  the screen (unless /silent is set).
;  =out         The IRAF output
;  =error       The error message if there was one.
;
; USAGE:
;  IDL>iraf_run,scriptname,irafdir,out=out
;
; By D.Nidever  August 2007
;  Similar to a perl script written by Armin Rest
;-

pro iraf_run,scriptname,irafdir,silent=silent,out=out,error=error

undefine,error

nscriptname = n_elements(scriptname)
nirafdir = n_elements(irafdir)

; Not enough inputs
if nscriptname eq 0 or nirafdir eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - iraf_run,scriptname,irafdir'
  return
endif

; Important directories
if n_elements(irafdir) eq 0 then irafdir='~/iraf/'
irafdir = FILE_SEARCH(irafdir,/fully_qualify,count=nirafdir)
if nirafdir eq 0 then begin
  error = 'NO IRAF DIRECTORY'
  print,'NO IRAF DIRECTORY'
  return
endif
CD,current=curdir

; Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL
;---------------------------------------------------------------
CHECK_IRAF,iraftest,irafdir=irafdir
if iraftest eq 0 then begin
  error = 'IRAF TEST FAILED'
  print,'IRAF TEST FAILED.  EXITING'
  return
endif


; Go to IRAF home directory
cd,irafdir

; Message to the screen
if not keyword_set(silent) then begin
  print,'Running IRAF...'
  print,'with IRAF home directory: ',irafdir
  print,'Executing script: ',scriptname
endif

; Execute the script
spawn,'cl < '+scriptname,out

; Print the output to the screen
if not keyword_set(silent) then $
  printline,out

; Go back to original directory
cd,curdir

if keyword_set(stp) then stop

end
