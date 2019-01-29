;+
;
; PRINTLOG
;
; This prings messages to the screen and to a logfile
;
; INPUTS:
;  logfile   The name of the logfile.  If logfile is not a string then
;             the messages are only printed to the screen.
;  d1-d50    Strings input to print.
;  =format   A format string.
;  /logonly  Print only to the logfile and not to the screen 
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The messages are printed to the screen and the logfile
;  =error    The error, if one occured, otherwise undefined.
;
; USAGE:
;  IDL>printlog,'log.txt','This needs to be logged'
;
; By D.Nidever  Feb. 2008   parts of it copied from Markwardt's PRINTLOG.PRO
;-

pro printlog,  logfile, d1,  d2,  d3,  d4,  d5,  d6,  d7,  d8,  d9, d10, $
               d11, d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, $
               d23, d24, d25, d26, d27, d28, d29, d30, d31, d32, d33, d34, $
               d35, d36, d37, d38, d39, d40, d41, d42, d43, d44, d45, d46, $
               d47, d48, d49, d50, d51, d52, d53, d54, d55, d56, d57, d58, $
               d59, d60, format=format,stp=stp,$
               logonly=logonly0,error=error

undefine,error

; Not enough inputs
nlogfile = n_elements(logfile)
if nlogfile eq 0 then begin
  print,'Syntax - printlog,logfile,messages,format=format,logonly=logonly,stp=stp'
  return
endif

type = size(logfile,/type)
if keyword_set(logonly0) then logonly=1 else logonly=0
if type ne 7 then logonly=0

; Too many parameters
np = n_params()-1   ; minus the logfile
if np GT 50 then $
  message, 'PRINTLOG ERROR: number of parameters to PRINTLOG cannot exceed 50'


; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   error = !ERROR_STATE.MSG
   PHOTRED_ERRORMSG
   CATCH, /CANCEL 
   return
endif



; Making command to create the string
if np gt 0 then begin
  ;cmd = string(lindgen(np)+1, $
  ;             format='("str = string(",50("D",I0,:,","))')
  ;if n_elements(format) GT 0 then cmd = cmd + ",format=format(0))" $
  ;else cmd = cmd + ")"
  ;
  ;str = ''
  ;result = execute(cmd)
  ;
  ;; Error
  ;if result NE 1 then begin
  ;  print,'ERROR:  Returning'
  ;  return
  ;endif

  ; Make a new format array
  if n_elements(format) gt 0 then begin
    ;fmtarr = strarr(np)
    format2 = REPSTR(format,'(','')
    format2 = REPSTR(format2,')','')
    format2arr = strsplit(format2,',',/extract)
    nformat2arr = n_elements(format2arr)
    for i=0,nformat2arr-1 do begin
      char = format2arr[i]
      chararr = strarr(strlen(char))
      for j=0,n_elements(chararr)-1 do chararr[j]=string( (byte(char))[j] )
      charnum = valid_num(chararr)
      firstletterind = first_el(where(charnum eq 0))
      ; There is a number out front, i.e. 2F10.2
      if firstletterind gt 0 then begin
        num = long(strmid(char,0,firstletterind))
        left = strmid(format2arr[i],firstletterind)
        push,fmtarr,replicate(left,num)
      endif else begin
        push,fmtarr,format2arr[i]
      endelse
    end
    fmtarr = '('+fmtarr+')'
    nfmtarr = n_elements(fmtarr)

    ; Compare to the number of elements to print
    nprint = 0
    for i=0,np-1 do $
      nprint = nprint + n_elements((SCOPE_VARFETCH('D'+strtrim(i+1,2))))

    ; Format string does not match input
    if nfmtarr ne nprint then begin
      print,'Format string does not match input'
      error = 'Format string does not match input'
      return
    endif

  ; No format string input
  endif else begin
    fmtarr = replicate('(1(A," "))',np)  ; This makes a space after each one
    ;if np eq 1 then fmtarr = ''
    ;fmtarr = replicate('',np)
  endelse

  ; Loop through each string element to print
  ; Some inputs may be arrays
  str = ''
  pcount = 0
  nfmtarr = n_elements(fmtarr)
  for i=0,np-1 do begin
    temp = (SCOPE_VARFETCH('D'+strtrim(i+1,2)))
    ntemp = n_elements(temp)
    for j=0,ntemp-1 do begin
      str=str+STRING( temp[j], format=fmtarr[pcount<(nfmtarr-1)])
      pcount++
    end
  end


  ; Only one input, print as a column
  if (np eq 1 and n_elements(format) eq 0) then begin
    temp = (SCOPE_VARFETCH('D1'))
    str = STRING(temp,format='(A)')
  endif


; Nothing input, empty string
endif else begin
  str = ''
endelse


; Print to logfile
if (type eq 7) then begin

  if not keyword_set(logonly) then PRINT,str,format='(A)'
  WRITELINE,logfile,str,/append

; Only print to screen
endif else begin

  PRINT,str,format='(A)'

endelse

if keyword_set(stp) then stop

end
