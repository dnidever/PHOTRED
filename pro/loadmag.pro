;+
;
; LOADMAG
;
; This loads the PHOTRED .mag files.
;
; INPUTS:
;  file    The .mag filename.
;
; OUTPUTS:
;  phot    The photometry structure
;  =error  The error message if one occurred.
;
; USAGE:
;  IDL>loadmag,'file.mag',phot
;
; By D.Nidever  made a separate program  July 2017
;-

pro loadmag,file,phot,error=error

undefine,error
undefine,phot
  
; Not enough inputs
if n_elements(file) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - loadmag,file,phot,error=error'
  return
endif

; Load the photometry file
;-------------------------
; This is copied from LOADRAW.PRO

; Figure out the number of columns
line1='' & line2='' & line3='' & line4='' & line5=''
openr,unit,/get_lun,file
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
WHILE (endflag ne 1) do begin

  line4 = ''
  readf,unit,line4

  ; If there are too many frames/columns then DAOMASTER puts
  ; these on separate lines and lead with ~27 blank spaces
  ;  Not sure if this is needed for MAG files as well

  ; This is a continuation line
  if strtrim(strmid(line4,0,15),2) eq '' then begin
    arr4 = strsplit(line4,' ',/extract)
    narr4 = n_elements(arr4)
    ncol += narr4
    nstarline++
  endif else endflag=1
ENDWHILE
close,unit
free_lun,unit

; Stars in this file
numstar = (FILE_LINES(file)-3L )/nstarline

; Number of observations
nextra = 7
numobs = (ncol-nextra)/2
;nextra = ncol - 2*numobs    ; nextra includes, ID, X, Y, CHI, SHARP, etc.

; mastable is where everything is stored, id, x, y, unsolved magnitudes, chi, sharp
mastable = fltarr(numstar,2*numobs+nextra)

; Reading in the magnitude file
get_lun,unit
openr, unit, file

line=''
readf,unit,line
readf,unit,line
readf,unit,line

; Loop through the stars
for j=0.,numstar-1 do begin
  instr=''
  instr1=''
  inline = fltarr(2*numobs+nextra)

  ; Loop through the lines per star
  for k=0l,nstarline-1 do begin
    readf, unit, instr1
    ;instr += instr1
    ; remove extra 25 spaces at the beginning of extra/wrap lines
    if k eq 0 then instr+=instr1 else instr+=strmid(instr1,25)
  endfor
  ; We need to use the formatted read because sometimes there are
  ; NO spaces between the numbers in the columns.
  ; ***BUT if the daomaster format changes then this will be MESSED UP!!!!***
  ; MAKEMAG.PRO uses a slightly different format than daomaster
  ;  more space for larger integers
  ;  fmt='(A1,I8,2F9.3,'+strtrim(nfiles*2+2,2)+'F9.4)'
  fmt='(I9,2F9.3,'+strtrim(ncol-3,2)+'F9.4)'   ; makemag output
  reads,instr,inline,format=fmt
  mastable[j,0:2*numobs+nextra-1] = inline[0:2*numobs+nextra-1]
endfor

close, unit
free_lun,unit

; Now create the structure
schema = {id:0L,x:0.0,y:0.0}
for m=0,numobs-1 do $
  schema = CREATE_STRUCT(schema,'MAG'+strtrim(m+1,2),0.0,'MAG'+strtrim(m+1,2)+'ERR',0.0)
schema = CREATE_STRUCT(schema,'CHI',0.0,'SHARP',0.0,'FLAG',0L,'PROB',0.0)
phot = REPLICATE(schema,numstar)
ncopy = 2*numobs+nextra
for j=0,ncopy-1 do phot.(j)=reform(mastable[*,j])

end
