;+
; This reads in all the lines from a file
;
; INPUTS:
;  file        File name.
;  =comment    Comment string.  Lines starting with this
;              character will be skipped
;  =nlineread  Only read this many lines.
;  /noblank    Remove blank lines and lines with only whitespace.
;  /stp        Stop at the end of the program
;
; OUTPUTS:
;  out         String array of the lines in the file
;  =count      The number of lines read in.  This is set to -1
;                if there was an error.
;
; USAGE:
;  IDL>readline,'file.txt',lines
;
; By D.Nidever  Feb.2007
;-

pro readline,file,out,comment=comment,nlineread=nlineread,count=count,noblank=noblank,stp=stp

undefine,out
count = 0

if n_elements(file) eq 0 then begin
  print,'Syntax - readline,file,out,comment=comment,nlineread=nlineread,count=count,'
  print,'                  noblank=noblank,stp=stp'
  return
endif

test = file_test(file)
if test eq 0 then begin
  print,'FILE ',file,' DOES NOT EXIST'
  return
endif


; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   undefine,out
   count = -1                ; There was an error
   PHOTRED_ERRORMSG
   CATCH, /CANCEL 
   return
endif

nlines = file_lines(file)
if nlines eq 0 then begin
  undefine,out
  count=0
  return
endif

OPENR,unit,/get_lun,file

line = ''
;out = strarr(10000L)
out = strarr(nlines)
count = 0.
endflag = 0
;while(~EOF(unit)) do begin
while (endflag eq 0) do begin

  ; Add more elements if necessary
  ;if count gt n_elements(out)-1 then out = [out,strarr(10000L)]

  ; Read the line
  readf,unit,line

  ; Checking for comment string
  if n_elements(comment) gt 0 then begin
    len = strlen(comment)
    dum = strmid(strtrim(line,2),0,len)
    
    ; Not commented
    if (dum ne comment) then begin

      ; Put the line in the array
      out[count] = line

      ; Increment
      count++

    end ; not commented

  endif else begin

    ; Put the line in the array
    out[count] = line

    ; Increment
    count++

  endelse ; no comment string

  ; Are we stopping
  if EOF(unit) eq 1 then endflag=1
  if keyword_set(nlineread) then if count ge nlineread then endflag=1

end

CLOSE,unit
FREE_LUN,unit


; Only want the elements that were read in
out = out[0:count-1]


; Remove blank lines
if keyword_set(noblank) then begin

  gd = where(strtrim(out,2) ne '',ngd)
  if ngd gt 0 then out=out[gd] else undefine,out

endif

count = n_elements(out)

if keyword_set(stp) then stop

end
