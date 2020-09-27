;+
;
; LOADOPT
;
; Load a DAOPHOT option file.
;
; INPUTS:
;  file    Option file name.
;
; OUTPUTS:
;  optstr  A structure with the value in the option file.
;  =error  The error, if one occurred.
;
; USAGE:
;  IDL>loadopt,optfile,optstr
;
; By D. Nidever  July 2020
;-

pro loadopt,file,optstr,error=error

undefine,optstr
nfile = n_elements(file)
;; Not enough inputs
if nfile eq 0 then begin
  print,'Syntax - loadopt,optfile,optstr,error=error'
  error = 'Not enough inputs'
  return
endif

;; Check that the file exists
if file_test(file) eq 0 then begin
  error = file+' NOT FOUND'
  print,error
  return
endif

;; Read in the lines
READLINE,file,optlines,count=noptlines
if noptlines eq 0 then begin
  error = file+' IS EMPTY'
  print,error
  return
endif
;; Remove any blank lines
bd = where(strtrim(optlines,2) eq '',nbd)
if nbd eq noptlines then begin
  error = file+' has only blank lines'
  print,error
  return
endif
if nbd gt 0 then REMOVE,bd,optlines
noptlines = n_elements(optlines)

;; Only keep lines with = in them
gd = where(strpos(optlines,'=') ne -1,ngd)
if ngd eq 0 then begin
  error = 'No key/value pairs found'
  print,error
  return
endif
optlines = optlines[gd]
noptlines = n_elements(optlines)

;; Split the lines into key/value pairs
;; RE =     0.01
;; GA =    28.67
;; LO =     7.00
dum = strsplitter(optlines,'=',/extract)
keys = reform(strtrim(dum[0,*],2))
values = reform(strtrim(dum[1,*],2))
nkeys = n_elements(keys)

;; Make the structure, assume they are all floats
optstr = create_struct(keys[0],float(values[0]))
for i=1,nkeys-1 do $
   optstr = create_struct(optstr,idl_validname(keys[i],/convert_all),float(values[i]))

end
