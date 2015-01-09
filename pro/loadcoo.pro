pro loadcoo,filename,str,head,stp=stp

;+
;
; LOADCOO
;
; This loads a DAOPHOT FIND coo file
;
; INPUTS:
;  filename   The name of the COO file
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  str        A structure with the COO data
;  head       A string array with the COO header
;
; USAGE:
;  IDL>loadcoo,'obj2153_1.coo',str,head
;
; By D. Nidever   Feb. 2008 (copied from loadals.pro)
;-

; Not enough inputs
if n_elements(filename) eq 0 then begin
  print,'Syntax - loadcoo,filename,str,head,stp=stp'
  return
endif

test = file_test(filename)
if test eq 0 then begin
  print,'FILE ',filename,' DOES NOT EXIST'
  str=-1
  return
endif

; Is this a COO file
line1='' & line2='' & line3='' & line4=''
openr,unit,/get_lun,filename
readf,unit,line1
readf,unit,line2
readf,unit,line3
readf,unit,line4
close,unit
free_lun,unit

arr4 = strsplit(line4,' ',/extract)
ncol = n_elements(arr4)

; This is a COO file 
arr1 = strsplit(line1,' ',/extract)
if arr1[0] eq 'NL' and strtrim(line3,2) eq '' then begin

  fields = ['ID','X','Y','MAG','SHARP','ROUND','ROUND2']
  ; ROUND2 is the Roundness index based on the marginal distribution
  types = [3,5,5,4,4,4,4]
  str = importascii(filename,fieldtype=types,fieldnames=fields,skip=3,/noprint)

  head = [line1,line2]

; This is NOT a COO file
endif else begin

  print,'This is NOT a DAOPHOT/FIND output file'
  return

endelse

if keyword_set(stp) then stop

end

