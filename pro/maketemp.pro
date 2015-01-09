function maketemp,prefix,suffix,ndig=ndig0,directory=directory

;+
;
; MAKETEMP.PRO
;
; PURPOSE:
;  This function creates a temporary filename.  It also checks
;  that this file does not exist already.  Seven random digits
;  are added to the prefix.
;
;  WARNING!!! Currently if two instances of MAKETEMP.PRO are run
;  at exactly the same time then they can output the *SAME*
;  "random" filename.  Not sure how to get around this.
;  Use MKTEMP.PRO instead!
;
; INPUTS:
;  prefix    The prefix for the temporary filename (e.g. "temp")
;  suffix    The suffix for the temporary filename (e.g. ".txt")
;  =ndig     How many digits to use. Default is 8.
;  /directory  This is for a directory.
;
; OUTPUTS:
;  temp      The filename of a temporary file
;
; USAGE:
;  IDL>tmp = maketemp('temp','.txt')
;
; By D.Nidever   August 2005
;-

dum = 'bad'

if n_elements(ndig0) eq 0 then ndig=8 else ndig=ndig0
ndig = ndig > 1

; Making sure it doesn't exist already
while (dum(0) ne '') do begin

  if n_elements(prefix) eq 0 then prefix='temp'
  if n_elements(suffix) eq 0 then suffix=''

  file = prefix
  for i=0,ndig-1 do file=file+strtrim(rndint(),2)
  file = file+suffix

  ; Checking if it exists already
  ;dum = FILE_SEARCH(file+'*')
  ;if keyword_set(directory) then dum = FILE_SEARCH(file+'*',/test_directory)
  dum = FILE_TEST(file)
  if keyword_set(directory) then dum = FILE_TEST(file,/directory)

end

; Create an empty file
openw,unit,file,/get_lun
close,unit
free_lun,unit

;stop

return,file

end
