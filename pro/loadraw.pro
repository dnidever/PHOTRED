pro loadraw,filename,phot,head,mastable=mastable,stp=stp

;+
;
; LOADRAW
;
; This loads the DAOPHOT/DAOMASTER RAW file
;
; INPUTS:
;  filename   The name of the ALS file
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  phot       A structure with the ALS data
;  head       A string array with the ALS header
;
; USAGE:
;  IDL>loadraw,'obj2153_1.raw',phot,head
;
; By D. Nidever   Feb. 2008 (basically a copy of loadals.pro
;-

; Not enough inputs
if n_elements(filename) eq 0 then begin
  print,'Syntax - loadraw,filename,phot,head,stp=stp'
  return
endif

test = file_test(filename)
if test eq 0 then begin
  print,'FILE ',filename,' DOES NOT EXIST'
  phot=-1
  return
endif

; Is this an ALS or PHOT file
line1='' & line2='' & line3='' & line4='' & line5=''
openr,unit,/get_lun,filename
readf,unit,line1
readf,unit,line2
readf,unit,line3
;readf,unit,line4
;readf,unit,line5
close,unit
free_lun,unit


; This is an ALS file 
arr1 = strsplit(line1,' ',/extract)
if arr1[0] eq 'NL' and strtrim(line3,2) eq '' then begin

  ; Figure out the number of columns
  line1='' & line2='' & line3='' & line4='' & line5=''
  openr,unit,/get_lun,filename
  readf,unit,line1
  readf,unit,line2
  readf,unit,line3

  ; First line for the first star
  readf,unit,line4
  ;arr4 = strsplit(line4,' ',/extract)
  ;narr4 = n_elements(arr4)
  ;ncol = narr4
  instr = line4

  ; Check for continuation lines
  endflag = 0
  nstarline = 1
  continuation = 0
  WHILE (endflag ne 1) and ~eof(unit) do begin

    line4 = ''
    readf,unit,line4

    ; If there are too many frames/columns then these
    ; go on separate lines and lead with 27 blank spaces

    ; This is a continuation line
    if strtrim(strmid(line4,0,15),2) eq '' then begin
      trial = strmid(line4,33,1)
      if trial eq ' ' then nspaces=24 else nspaces=25
      instr1 = strmid(line4,nspaces)
      instr += instr1
      ;arr4 = strsplit(line4,' ',/extract)
      ;narr4 = n_elements(arr4)
      ;ncol += narr4
      nstarline++
      continuation=1
    endif else endflag=1
  ENDWHILE
  close,unit
  free_lun,unit

  ;; Now parts the long line for a single star with a formatted read
  ;fmt = '(I7,2A9,'+strtrim(2*nmag+2,2)+'F9.4)'
  nmag = (strlen(instr)-(7+2*9+2*9)) / 9 / 2
  ncol = nmag*2+5

  ; ID  X  Y  MAG1  ERR1  MAG2  ERR2 ...  CHI SHARP
  ;nmag = (ncol-5)/2

  ; Stars in this file
  numstar = (FILE_LINES(filename)-3L )/nstarline

  ; LOADING THE DATA
  ;------------------

  ;; Cannot use importascii.pro to load the .raw files because
  ;; sometimes the columns run together (e.g. large negative numbers)
  ;; and then importascii.pro has problems.  We need a formatted read.

  ;; mastable is where everything is stored, id, x, y, unsolved magnitudes, chi, sharp
  mastable = fltarr(numstar,2*nmag+5)

  ;; Reading in the magnitude file
  openr, unit,/get_lun,filename

  line1='' & lines2='' & lines3=''
  readf,unit,line1
  readf,unit,line2
  readf,unit,line3
  head = [line1,line2]

  ;; Loop through the stars
  for j=0l,numstar-1 do begin
    ;instr=' '
    instr=''
    instr1=' '
    inline = fltarr(2*nmag+5)
      
    ;; Loop through the lines per star
    for k=0l,nstarline-1 do begin
      readf, unit, instr1
      if k gt 0 then begin
        ;; There are leading spaces, 24 or 25
        ;; Use the first character AFTER the first column to figure out
        ;;   how many spaces we need to strip off
        ;; The MAG/ERR columns are F9.4, so if there are 25 spaces then
        ;;   the 34th character (index=33) should be the last digit of MAG1 (25+9=34)
        ;;   and be a number.  If there are 24 leading space then the
        ;;   34th character will right after MAG1 and be a space.
        trial = strmid(instr1,33,1)
        if trial eq ' ' then nspaces=24 else nspaces=25
        instr1 = strmid(instr1,nspaces)
      endif
      ;if k gt 0 then instr1=strmid(instr1,25) ; 2nd and later lines have 25 leading spaces
      instr += instr1
    endfor

    ; WRITE (1,111) IDMAST(IMASTR), POSIT1, POSIT2,
    ; .            ((DATUM(J,I), J=1,2), I=1,NFRM), SUMCHI, SUMSHP
    ;        111       FORMAT (I7, 2A9, 12F9.4:/ (25X, 12F9.4))

    ;; formatted read
    fmt = '(I7,2A9,'+strtrim(2*nmag+2,2)+'F9.4)'
    id = 0L & x='' & y='' & all=fltarr(2*nmag+2)
    reads,instr,id,x,y,all,format=fmt
    inline[0] = id
    inline[1] = x
    inline[2] = y
    inline[3:2*nmag+5-1] = all

    ; old method, causes problems when there are no spaces between values
    ;reads,instr,inline

    mastable[j,0:2*nmag+5-1] = inline[0:2*nmag+5-1]
  endfor

  close, unit
  free_lun,unit

  ; Now transfer to structure
  fieldtypes = [3,lonarr(4+2*nmag)+4]
  fieldnames = ['ID','X','Y']
  for i=1,nmag do fieldnames = [fieldnames,'MAG'+strtrim(i,2),'ERR'+strtrim(i,2)]
  fieldnames = [fieldnames,'CHI','SHARP']

  mastable2 = transpose(mastable)
  if numstar gt 1 then begin
    phot = ARR2STR(mastable2,fieldnames=fieldnames,fieldtypes=fieldtypes,/noprint)
  endif else begin
    ; arr2str needs a 2D array, make it [Ncol, 1]
    phot = ARR2STR(reform(mastable2,n_elements(mastable2),1),fieldnames=fieldnames,fieldtypes=fieldtypes,/noprint)
  endelse


; Not a RAW file
endif else begin

  print,'This is NOT a RAW file'
  return

endelse

if keyword_set(stp) then stop

end

