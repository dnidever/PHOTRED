;+
;
; LOADALS
;
; This loads the DAOPHOT photometry file
;
; INPUTS:
;  filename   The name of the ALS file
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  phot       A structure with the ALS data
;  head       A string array with the ALS header
;  =count     The number of entries.  count=0 if there
;               was a problem.
;
; USAGE:
;  IDL>loadals,'obj2153_1.als',phot,head
;
; By D. Nidever   January 2007
;-

pro loadals,filename,phot,head,count=count,stp=stp

count=0

; Not enough inputs
if n_elements(filename) eq 0 then begin
  print,'Syntax - loadals,filename,phot,head,count=count,stp=stp'
  return
endif

test = file_test(filename)
if test eq 0 then begin
  print,'FILE ',filename,' DOES NOT EXIST'
  phot=-1
  return
endif

nlines = file_lines(filename)
if nlines lt 4 then begin
  print,'No sources in '+filename
  phot = -1
  return
endif

; Is this an ALS or PHOT file
line1='' & line2='' & line3='' & line4=''
openr,unit,/get_lun,filename
readf,unit,line1
readf,unit,line2
readf,unit,line3
readf,unit,line4
close,unit
free_lun,unit

; This is an ALS file 
arr1 = strsplit(line1,' ',/extract)
if arr1[0] eq 'NL' and strtrim(line3,2) eq '' then begin

  ;arr4 = strsplit(line4,' ',/extract)
  ;narr4 = n_elements(arr4)
  ;nmag = (narr4-5)/2

  ;fieldtypes = [3,lonarr(narr4-1)+4]
  ;fieldnames = ['ID','X','Y']
  ;for i=1,nmag do fieldnames = [fieldnames,'MAG'+strtrim(i,2),'MERR'+strtrim(i,2)]
  ;fieldnames = [fieldnames,'CHI','SHARP']

  fields = ['ID','X','Y','MAG','ERR','SKY','ITER','CHI','SHARP']
  types = [3,4,4,4,4,4,4,4,4]
  phot = importascii(filename,fieldtype=types,fieldnames=fields,skip=3,/noprint)

  head = [line1,line2]

; This is a PHOT file
endif else begin

  print,'This is NOT an ALLSTAR output file'
  return

  ;; Read the photometry file
  ;phot = importascii(filename,/header,/noprint)

endelse

count = n_elements(phot)

if keyword_set(stp) then stop

end

