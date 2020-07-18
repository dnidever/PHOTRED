pro loadmch,mchfile,files,trans,magoff,count=count,stp=stp

;+
;
; LOADMCH
;
; This loads a DAOMATCH/DAOMASTER mch file
;
; INPUTS:
;  mchfile  The MCH filename
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  files    The list of files in the MCH file
;  trans    The transformation equations in a Nfilesx6 array.
;  magoff   The magnitude offsets and errors array.
;  =count   The number of files.
;
; USAGE:
;  IDL>loadmch,'ccd1001.mch',files,trans,count=count
;
; By D.Nidever   February 2008
;-

count=0

; Enough inputs
nmchfile = n_elements(mchfile)
if nmchfile eq 0 then begin
  print,'Syntax - loadmch,mchfile,files,trans,magoff,stp=stp'
  return
endif

; Test the file
test = file_test(mchfile)
if test eq 0 then begin
  print,mchfile,' NOT FOUND'
  return
endif

openr,unit,/get_lun,mchfile

while(~EOF(unit)) do begin
  line=''
  readf,unit,line
  push,lines,line
endwhile

close,unit
free_lun,unit

; Creating the trans array
lines2 = repstr(lines,"'",'')
nlines = n_elements(lines)
arr = strsplit(lines2[0],' ',/extract)
ntrans = n_elements(arr)-3    ; first line is the samee, last two are mag offsets

; Getting the file names
arr2 = strsplitter(lines2,' ',/extract)
files = reform(arr2[0,*])

; Initializing the array
trans = dblarr(nlines,ntrans)
; Filling the aray
for i=0,nlines-1 do begin
  arr = strsplit(lines2[i],' ',/extract)
  trans[i,*] = arr[1:ntrans]
endfor

; The magnitude offsets
magoff = float(arr2[ntrans+1:*,*])

count = nlines

if keyword_set(stp) then stop

end
