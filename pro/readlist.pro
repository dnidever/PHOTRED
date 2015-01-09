pro readlist,filename,lines,exist=exist,unique=unique,count=count,$
             fully=fully,logfile=logfile,silent=silent,stp=stp

;+
;
; READLIST
;
; This reads in a list of filenames for photred
;
; INPUTS:
;  filename   The name of the list file to load
;  /exist     Test that all of the files exist
;  /unique    Only return unique filenames
;  /fully     Return fully qualified paths (i.e. full absolute paths).
;  =logfile   A logfile to print messages to as well as the screen
;  /silent    Don't print anything.
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  lines      The filenames read in
;  =count     The number of lines read in
;
; USAGE:
;  IDL>readlist,'inputlist.txt',lines,/test,/unique,/fully,logfile=logfile,count=count
;
; By D.Nidever  Feb 2008
;-

; Not enough inputs
nfilename = n_elements(filename)
if nfilename eq 0 then begin
  print,'Syntax - readlist,filename,lines,test=test,unique=unique,fully=fully,'
  print,'                  logfile=logfile,count=count,stp=stp'
  return
endif

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1

undefine,lines,count

; Does the file exist
test = file_test(filename)
if test eq 0 then begin
  if not keyword_set(silent) then $
    printlog,logf,'FILE ',filename,' NOT FOUND'
  count=0
  return
endif

; Read in the file
if numlines(filename) gt 0 then $
;READLINE,filename,lines,comment='#'
;lines = strtim(lines,2)
READCOL,filename,lines,format='A',comment='#',/silent
nlines = n_elements(lines)

; No lines
if nlines eq 0 then begin
  if not keyword_set(silent) then $
    printlog,logf,'NO LINES IN FILE ',filename
  count=0
  return
endif


; Fully qualify path
if keyword_set(fully) and (nlines gt 0) then begin
  old_lines = lines

  lines = FILE_SEARCH(lines,/fully_qualify)
  gd = where(strtrim(lines,2) ne '',ngd)

  ; Some good files
  if ngd gt 0 then begin
    lines = lines[gd]
    nlines = ngd
    count = ngd
  endif else begin
    if not keyword_set(silent) then $
      printlog,logf,'None of the files exist'
    undefine,lines
    nlines = 0
    count=0
  endelse

endif

; Unique files
if keyword_set(unique) and (nlines gt 0) then begin
  ui = uniq(lines,sort(lines))
  si = sort(ui)   ; Keep them in the original order
  ui = ui[si]
  lines = lines[ui]
endif

; Test the files exist
if keyword_set(exist) and (nlines gt 0) then begin
  test = file_test(lines)
  gd = where(test eq 1,ngd)

  ; Some good files
  if ngd gt 0 then begin
    lines = lines[gd]
    nlines = 0
    count = ngd
  endif else begin
    if not keyword_set(silent) then $
      printlog,logf,'None of the files exist'
    undefine,lines
    nlines = 0
    count=0
  endelse

endif


count = n_elements(lines)

if keyword_set(stp) then stop

end
