pro loadaper,filename,phot,head,stp=stp

;+
;
; LOADAPER
;
; This loads the DAOPHOT aperture photometry file
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
;  IDL>loadaper,'obj2153_1.als',phot,head
;
; By D. Nidever   January 2008
;-

; Not enough inputs
if n_elements(filename) eq 0 then begin
  print,'Syntax - loadaper,filename,phot,head,stp=stp'
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


; This is an APER file 
arr1 = strsplit(line1,' ',/extract)
arr2 = strsplit(line2,' ',/extract)
if arr1[0] eq 'NL' and strtrim(arr2[0],2) eq '2' and strtrim(line3,2) eq '' then begin

  ; Making the header
  head = [line1,line2]

  ; How many columns
  arr5 = strsplit(line5,' ',/extract)
  ncol = n_elements(arr5)
  naper = ncol-3

  ; How many stars
  nlines = FILE_LINES(filename)
  nstars = (nlines-3.0)/3.0

  ; Starting the photometry structure
  cmd = 'dum={id:0L, x:0.0, y:0.0, sky:0.0, skysig:0.0, skyskew:0.0, mag:fltarr('+strtrim(naper,2)+')'
  cmd = cmd+', err:fltarr('+strtrim(naper,2)+')}'
  out = EXECUTE(cmd)
  phot = replicate(dum,nstars)

  ; This is what a Daophot aperture photometry file looks like
  ; NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
  ;  2  2040  2047  1837.1 26000.0  114.56    4.00    2.20    2.86    1.98
  ; 
  ;
  ;      1    3.000    1.000   97.999
  ;      2055.039 41.29  0.00  9.9999
  ;
  ;      2  397.000    1.000   97.999
  ;      2066.848 37.25  0.12  9.9999
  ;
  ;  The columns are: ID, X, Y, Mag1, Mag2, etc..
  ;                   Sky, St.Dev. of sky, skew of sky, Mag1err, Mag2err, etc.
  ; 

  OPENR,unit,/get_lun,filename

  ; Read first three lines
  dum=''
  readf,unit,dum
  readf,unit,dum
  readf,unit,dum

  fmt1 = '(I7,F9.3,F9.3,'+strtrim(naper,2)+'F9.3)'
  fmt2 = '(F14.3,F6.2,F6.2,F8.4'
  if naper gt 1 then fmt2=fmt2+','+strtrim(naper-1,2)+'F9.4)' else fmt2=fmt2+')'

  count = 0LL
  WHILE (~EOF(unit)) do begin

    blank=''
    id=0L & x=0.0 & y=0.0
    sky=0.0 & skysig=0.0 & skyskew=0.0
    mag = fltarr(naper)
    err = fltarr(naper)
    readf,unit,blank   ; First line is blank
    readf,unit,id,x,y,mag,format=fmt1
    readf,unit,sky,skysig,skyskew,err,format=fmt2

    phot[count].id = id
    phot[count].x = x
    phot[count].y = y
    phot[count].sky = sky
    phot[count].skysig = skysig
    phot[count].skyskew = skyskew
    
    phot[count].mag = mag
    phot[count].err = err

    count++

  END

  CLOSE,unit
  FREE_LUN,unit


; This is a PHOT file
endif else begin

  print,'This is NOT a Daophot aperture photometry output file'
  return

  ;; Read the photometry file
  ;phot = importascii(filename,/header,/noprint)

endelse

if keyword_set(stp) then stop

end

