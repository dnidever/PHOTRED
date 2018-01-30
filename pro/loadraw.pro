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
readf,unit,line4
readf,unit,line5
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
  arr4 = strsplit(line4,' ',/extract)
  narr4 = n_elements(arr4)
  ncol = narr4

  ; Check for continuation lines
  endflag = 0
  nstarline = 1
  continuation = 0
  WHILE (endflag ne 1) do begin

    line4 = ''
    readf,unit,line4

    ; If there are too many frames/columns then these
    ; go on separate lines and lead with 27 blank spaces

    ; This is a continuation line
    if strtrim(strmid(line4,0,15),2) eq '' then begin
      arr4 = strsplit(line4,' ',/extract)
      narr4 = n_elements(arr4)
      ncol += narr4
      nstarline++
      continuation=1
    endif else endflag=1
  ENDWHILE
  close,unit
  free_lun,unit

  ; ID  X  Y  MAG1  ERR1  MAG2  ERR2 ...  CHI SHARP
  nmag = (ncol-5)/2


  ; Stars in this file
  numstar = (FILE_LINES(filename)-3L )/nstarline

  ; LOADING THE DATA
  ;------------------

  ; nfiles < 12
  ; Each star has data on 1 line
  ;if (nmag lt 12) then begin
  if (continuation eq 0) then begin

    fieldtypes = [3,lonarr(ncol-1)+4]
    fieldnames = ['ID','X','Y']
    for i=1,nmag do fieldnames = [fieldnames,'MAG'+strtrim(i,2),'ERR'+strtrim(i,2)]
    fieldnames = [fieldnames,'CHI','SHARP']

    phot = importascii(filename,fieldtype=fieldtypes,fieldnames=fieldnames,skip=3,/noprint)
    head = [line1,line2]

  ; nfiles >= 12
  ; Each star has data on 2 lines
  ; Most of this code was copied from photcalib.pro
  endif else begin

    ; mastable is where everything is stored, id, x, y, unsolved magnitudes, chi, sharp
    mastable = fltarr(numstar,2*nmag+5)

    ; Reading in the magnitude file
    openr, unit,/get_lun,filename

    line1='' & lines2='' & lines3=''
    readf,unit,line1
    readf,unit,line2
    readf,unit,line3
    head = [line1,line2]

    ; Loop through the stars
    for j=0l,numstar-1 do begin
      ;instr=' '
      instr=''
      instr1=' '
      inline = fltarr(2*nmag+5)
      
      ; Loop through the lines per star
      for k=0l,nstarline-1 do begin
        readf, unit, instr1
        if k gt 0 then instr1=strmid(instr1,25) ; 2nd and later lines have 25 leading spaces
        instr += instr1
      endfor

      ; WRITE (1,111) IDMAST(IMASTR), POSIT1, POSIT2,
      ; .            ((DATUM(J,I), J=1,2), I=1,NFRM), SUMCHI, SUMSHP
      ;        111       FORMAT (I7, 2A9, 12F9.4:/ (25X, 12F9.4))

      ; formatted read
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
    fieldtypes = [3,lonarr(nmag-1)+4]
    fieldnames = ['ID','X','Y']
    for i=1,nmag do fieldnames = [fieldnames,'MAG'+strtrim(i,2),'ERR'+strtrim(i,2)]
    fieldnames = [fieldnames,'CHI','SHARP']

    mastable2 = transpose(mastable)
    phot = ARR2STR(mastable2,fieldnames=fieldnames,fieldtypes=fieldtypes,/noprint)

  endelse ; nfiles >=12


; Not a RAW file
endif else begin

  print,'This is NOT a RAW file'
  return

endelse

if keyword_set(stp) then stop

end

