pro writeline,file,lines,stp=stp,prepend=prepend,append=append

;+
;
; WRITELINE
;
; This writes a string array to a file. Reverse of READLINE.PRO
;
; INPUTS:
;  file      File name.
;  lines     String array
;  =comment  Comment string.  Lines starting with this
;            character will be skipped
;  /prepend  Prepend to file (i.e. at the beginning)
;  /append   Append to file (i.e at the end)
;  /stp      Stop at the end of the program
;
; OUTPUTS:
;  The output written to the file
;
; USAGE:
;  IDL>writeline,'test.txt',lines
;
; By D.Nidever  Feb.2007
;-

if n_elements(file) eq 0 then begin
  print,'Syntax - writeline,file,lines,stp=stp'
  return
endif

info = file_info(file)

; File does *not* exist
if (info.exists eq 0) then begin

  OPENW,unit,/get_lun,file

  nlines = n_elements(lines)
  for i=0.,nlines-1 do printf,unit,lines[i]

  CLOSE,unit
  FREE_LUN,unit

; File exists
endif else begin

  ; Prepend
  if keyword_set(prepend) and not keyword_set(append) then begin

    ; Write new lines to a temporary file
    ;temp1file = maketemp('temptemp','temp')
    temp1file = MKTEMP('temp')
    OPENW,unit,/get_lun,temp1file

    nlines = n_elements(lines)
    for i=0.,nlines-1 do printf,unit,lines[i]

    CLOSE,unit
    FREE_LUN,unit

    ; Move original file to second temporary file
    ;temp2file = maketemp('temptemp2','temp')
    temp2file = MKTEMP('temp')
    FILE_MOVE,file,temp2file,/over

    ; Concatenate them, new lines first
    SPAWN,'cat '+temp1file+' '+temp2file+' > '+file,out

    ; Erase temporary file
    FILE_DELETE,temp1file
    FILE_DELETE,temp2file

  endif

  ; Normal operation or append
  if not keyword_set(prepend) then begin

    OPENW,unit,/get_lun,file,append=append

    nlines = n_elements(lines)
    for i=0.,nlines-1 do printf,unit,lines[i]

    CLOSE,unit
    FREE_LUN,unit

  endif


endelse

if keyword_set(stp) then stop

end
