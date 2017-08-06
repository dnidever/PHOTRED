;+
;
; READ_TRANS
;
; Read in a photometric transformation file.
;
; INPUTS:
;  transfile  This gives the transformation information needed
;             to calibrate the raw photometry.  Normally there
;             is a different one of these for every night.
;
;             There need to be two lines per band.
;             FIRST line:  band name,  color name, transformation
;             equation coefficients (zero point, airmass, color
;             airmass*color, color^2)
;             SECOND line: errors in the transformation equation
;             coefficients
;
;     This is an example transfile:
;     M    M-T  -0.9990    0.1402     -0.1345    0.0000   0.0000
;               1.094E-02  5.037E-03  2.010E-03  0.0000   0.0000
;     T    M-T  -0.0061    0.0489     0.0266     0.0000   0.0000
;               6.782E-03  3.387E-03  1.374E-03  0.0000   0.0000
;     D    M-D  1.3251     0.1403     -0.0147    0.0000   0.0000
;               1.001E-02  5.472E-03  2.653E-02  0.0000   0.0000
;
;      If the transfile has chip information then it should
;      look like this:
;  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
;              0.0040   -0.0000    0.0001   0.0000   0.0000
; 
;  2  G  G-R  -0.3617    0.1713   -0.1193   0.0000   0.0000
;              0.0039   -0.0000    0.0001   0.0000   0.0000
; 
;  3  G  G-R  -0.3457    0.1713   -0.1193   0.0000   0.0000
;              0.0039   -0.0000    0.0001   0.0000   0.0000
;
;      If the transfile has night and chip information then it should
;      look like this:
;  55975  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
;                     0.0040   -0.0000    0.0001   0.0000   0.0000
; 
;  55975  2  G  G-R  -0.3617    0.1713   -0.1193   0.0000   0.0000
;                     0.0039   -0.0000    0.0001   0.0000   0.0000
; 
;  55975  3  G  G-R  -0.3457    0.1713   -0.1193   0.0000   0.0000
;                     0.0039   -0.0000    0.0001   0.0000   0.0000
;
;      The transformation information can also be input for
;      individual file:
;  F5-00517150_43  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
;                           0.0040   -0.0000    0.0001   0.0000   0.0000
;
;  F5-00517150_44  G  G-R  -0.4284    0.1767   -0.1261   0.0000   0.0000
;                           0.0030    0.0020    0.0001   0.0000   0.0000
;
;  F5-00517150_45  G  G-R  -0.5082    0.1801   -0.1215   0.0000   0.0000
;                           0.0025    0.0023    0.0001   0.0000   0.0000
;
;  Each pair of lines can use a different format, i.e. the four format
;  types can be mixed in a file.  But this is not recommended since
;  it might be confusing what transformation to use for a given file
;  since there could be multiple "matches".
;
;  /silent    Don't print anything to the screen.
;  =logfile   The name of a logfile to write messages to.
;  /stp       Stop at the end of the program.
;
; OUTPUTS:
;  trans      The transformation structure.  NIGHT and CHIP will
;               always be included even if they are "blank".
;
; USAGE:
;  IDL>read_trans,'n1.trans',trans
;
; By D.Nidever  Feb.2013
;    Added Night+chip format March 2015
;    Added filename format  July 2017
;-

pro read_trans,transfile,trans,silent=silent,logfile=logfile,error=error,stp=stp

undefine,trans

if n_elements(transfile) eq 0 then begin
  print,'Syntax - read_trans,transfile,trans,silent=silent,logfile=logfile,error=error,stp=stp'
  return
endif

if file_test(transfile) eq 0 then begin
  print,transfile,' NOT FOUND'
  return
endif

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1


;# #####################################################
;# READ THE TRANSFORMATION FILE
;# Two lines per band.
;# First line:  band name,  color name, transformation equation coefficients
;# Second line: errors in the transformation equation coefficients

; If this has chip-specific transformations then the lines will be
; First line:  chip,  band name, color name, trans eqns.
; second line:  errors
;  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
;              0.0040   -0.0000    0.0001   0.0000   0.0000
; 
;  2  G  G-R  -0.3617    0.1713   -0.1193   0.0000   0.0000
;              0.0039   -0.0000    0.0001   0.0000   0.0000
; 
;  3  G  G-R  -0.3457    0.1713   -0.1193   0.0000   0.0000
;              0.0039   -0.0000    0.0001   0.0000   0.0000
;
; If the file has night and chip information then the lines will be
; First line: night, chip,  band name, color name, trans eqns.
; second line:  errors
;  55975  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
;                     0.0040   -0.0000    0.0001   0.0000   0.0000
; 
;  55975  2  G  G-R  -0.3617    0.1713   -0.1193   0.0000   0.0000
;                     0.0039   -0.0000    0.0001   0.0000   0.0000
; 
;  55975  3  G  G-R  -0.3457    0.1713   -0.1193   0.0000   0.0000
;                     0.0039   -0.0000    0.0001   0.0000   0.0000
;
; If the file has information on individual files then the lines will be
; First line: filename,  band name, color name, trans eqns.
; second line:  errors
;  F5-00517150_43  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
;                           0.0040   -0.0000    0.0001   0.0000   0.0000
;
;  F5-00517150_44  G  G-R  -0.4284    0.1767   -0.1261   0.0000   0.0000
;                           0.0030    0.0020    0.0001   0.0000   0.0000
;
;  F5-00517150_45  G  G-R  -0.5082    0.1801   -0.1215   0.0000   0.0000
;                           0.0025    0.0023    0.0001   0.0000   0.0000
;

; What options?
; 1 - single night, single chip (no NIGHT or CHIP information)
; 2 - single night, multiple chips (no NIGHT but CHIP information)
; 3 - multiple nights, multiple chips (NIGHT and CHIP information)
; 4 - individual files
optcase = -1
; I DON'T THINK "optcase" IS ACTUALLY USED FOR ANYTHING ANYMORE

openr,unit,/get_lun,transfile

while (~EOF(unit)) do begin

  trans1 = {night:-1L,chip:-1,file:'',band:'',color:'',colband:'',colsign:0,zpterm:0.0d,amterm:0.0d,colterm:0.0d,$
            amcolterm:0.0d,colsqterm:0.0d,zptermsig:0.0d,amtermsig:0.0d,coltermsig:0.0d,$
            amcoltermsig:0.0d,colsqtermsig:0.0d}

  ; Reading in the transformation coefficients line
  line = ''
  readf,unit,line

  ; Not a blank line
  if strtrim(line,2) ne '' and strmid(line,0,1) ne '#' then begin
    arr = strsplit(line,' ',/extract)
    narr = n_elements(arr)

    ; This is the format.  NIGHT (MJD), CHIP, BAND, COLOR, ZPTERM, AMTERM,
    ;                       COLORTERM, AMCOLTERM, AMSQTERM
    ;  55975  1  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
    ;                     0.0040   -0.0000    0.0001   0.0000   0.0000
    ; NIGHT and CHIP are optional.
    ; For specific files it looks like this:
    ;  F5-00517150_43  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
    ;                           0.0040   -0.0000    0.0001   0.0000   0.0000

    
    ; Are NIGHT and CHIP supplied?
    case narr of
       ; ---NIGHT and CHIP---
       9: begin

         ; Make sure NIGHT and CHIP are numbers
         if valid_num(arr[0],night) and valid_num(arr[1],chip) then begin
           trans1.night = long(night)
           trans1.chip = long(chip)
           arr = arr[2:*]
         ; Parsing error
         endif else begin
           error = 'NIGHT and CHIP must be numbers.'
           if not keyword_set(silent) then printlog,logf,error
           return
         endelse
         optcase = 3 > optcase

       end
       ; ---CHIP or FILENAME---
       8: begin

         ; CHIP, first value is a number
         if valid_num(arr[0],chip) then begin
           trans1.chip = long(chip)
           arr = arr[1:*]
           optcase = 2 > optcase
         ; FILENAME, first value is NOT a number
         endif else begin
           trans1.file = strtrim(arr[0],2)
           arr = arr[1:*]
           optcase = 4 > optcase
         endelse

       end
       ; ---Neither---
       7: begin
          ; not much to do
          optcase = 1 > optcase
       end
       ; ---Not enough values---
       else: begin
         error = 'Need at least 7 values in the TRANS line.'
         if not keyword_set(silent) then printlog,logf,error
         return
       end
    endcase
       
    ; Parse the rest of the line    
    trans1.band = arr[0]
    trans1.color = arr[1]
    trans1.zpterm = arr[2]
    trans1.amterm = arr[3]
    trans1.colterm = arr[4]
    trans1.amcolterm = arr[5]
    trans1.colsqterm = arr[6]

    ; Reading in the error line
    ;---------------------------
    line2 = ''
    readf,unit,line2
    arr2 = strsplit(line2,' ',/extract)
    narr2 = n_elements(arr2)
    
    ; Need at least 5 terms
    if narr2 lt 5 then begin
      error = 'Need at least 5 values in the ERROR line.'
      if not keyword_set(silent) then printlog,logf,error
      return
    endif

    ; Parse the error line
    trans1.zptermsig = arr2[0]
    trans1.amtermsig = arr2[1]  
    trans1.coltermsig = arr2[2] 
    trans1.amcoltermsig = arr2[3]
    trans1.colsqtermsig = arr2[4]

    ; Add to final transformation structure
    push,trans,trans1

  endif

endwhile

close,unit
free_lun,unit

; Leave in NIGHT and CHIP columns even if they are "blank".  For
; consistency sake

;; No chip information, strip CHIP
;gdnight = where(trans.chip ge 0,ngdchip)
;gdchip = where(trans.chip ge 0,ngdchip)
;if ngdchip eq 0 then begin
;  oldtrans = trans
;  trans = replicate({band:'',color:'',colband:'',colsign:0,zpterm:0.0d,amterm:0.0d,colterm:0.0d,$
;            amcolterm:0.0d,colsqterm:0.0d,zptermsig:0.0d,amtermsig:0.0d,coltermsig:0.0d,$
;            amcoltermsig:0.0d,colsqtermsig:0.0d},n_elements(trans))
;  STRUCT_ASSIGN,oldtrans,trans
;  undefine,oldtrans
;endif

ntrans = n_elements(trans)


; Figure out the colband and colsign for each band/chip
for i=0,ntrans-1 do begin

  band = strtrim(trans[i].band,2)

  col = strcompress(trans[i].color,/remove_all)

  ; Splitting up the two bands
  arr = strsplit(col,'-',/extract)

  ind = where(arr eq band,nind)

  ; colsign = 1 means band - colband
  if (ind[0] eq 0) then begin
    trans[i].colband = arr[1]
    trans[i].colsign = 1
  endif

  ; colsign = 2 means colband - band
  if (ind[0] eq 1) then begin
    trans[i].colband = arr[0]
    trans[i].colsign = 2
  endif

  if (ind[0] eq -1) then begin
    trans[i].colband = ''
    trans[i].colsign = -1
  endif

endfor


; Print the transformation equations
if not keyword_set(silent) then begin
  printlog,logf,' TRANSFORMATION EQUATIONS'
  printlog,logf,'--------------------------------------------------------------------------------'
  printlog,logf,'  NIGHT/CHIP/FILE  BAND COLOR ZERO-POINT  AIRMASS   COLOR     AIR*COL   COLOR^2 '
  printlog,logf,'--------------------------------------------------------------------------------'
  for i=0,ntrans-1 do begin
    form1 = '(I10,I6,A6,A7,F10.4,F10.4,F10.4,F10.4,F10.4)'
    form1f = '(A-16,A6,A7,F10.4,F10.4,F10.4,F10.4,F10.4)'
    ; FILENAME
    if trans[i].file ne '' then $
      printlog,logf,format=form1f,trans[i].file,'  '+trans[i].band,trans[i].color,trans[i].zpterm,$
                        trans[i].amterm,trans[i].colterm,trans[i].amcolterm,trans[i].colsqterm
    ; NO Filename
    if trans[i].file eq '' then $
      printlog,logf,format=form1,trans[i].night,trans[i].chip,'  '+trans[i].band,trans[i].color,trans[i].zpterm,$
                        trans[i].amterm,trans[i].colterm,trans[i].amcolterm,trans[i].colsqterm
    form2 = '(A29,F10.4,F10.4,F10.4,F10.4,F10.4)'
    printlog,logf,format=form2,'',trans[i].zptermsig,trans[i].amtermsig,trans[i].coltermsig,$
                      trans[i].amcoltermsig,trans[i].colsqtermsig
  endfor
  printlog,logf,'--------------------------------------------------------------------------------'
  printlog,logf,''
endif

if keyword_set(stp) then stop

end
