pro photred_loadsetup,setup,count=count,std=std,setupdir=setupdir,fake=fake,stp=stp

;+
;
; PHOTRED_LOADSETUP
;
; INPUTS:
;  /std       Load the stdred.setup file.
;  /fake      Load the fakered.setup file.
;  =setupdir  The directory in which to look for the setup file.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  setup      The setup file.  It is a string array with
;               dimensions of 2xN_parameters.  READPAR can be
;               used to read the parameters.
;  =count     The number of parameters.  count is -1 if there
;               was a problem.
;
; USAGE:
;  IDL>photred_loadsetup,setup,count=count,std=std,fake=fake,setupdir=setupdir,stp=stp
;
; By D.Nidever  March 2008
;-

count = 0
setup = strarr(2,2)+'-1'         ; no setup

if n_elements(setupdir) eq 0 then setupdir = '.'  ; default setup directory

; Type of setup file
setupfile = 'photred'
if keyword_set(std) then setupfile = 'stdred'
if keyword_set(fake) then setupfile = 'fakered'

; LOAD THE SETUP FILE
;--------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
setupfiles = FILE_SEARCH(setupdir+'/'+setupfile+'.*setup',count=nsetupfiles)
if (nsetupfiles lt 1) then begin
  print,'NO '+strupcase(setupfile)+' SETUP FILE'
  count = -1             ; there was a problem
  return
endif
if (nsetupfiles gt 1) then begin
  print,'MORE THAN ONE '+strupcase(setupfile)+' SETUP FILE FOUND'
  print,setupfiles
  count = -1              ; there was a problem
  return
endif

; Read the setup file
READLINE,setupfiles[0],lines,comment='#',count=nlines

; Process the lines
if nlines gt 0 then begin
  lines2 = strsplitter(lines,' ',/extract)
  lines2 = strtrim(lines2,2)
  sz = size(lines2)
  ncol = sz[1]
  npar = sz[2]

  ; Only one, add a blank second column
  if ncol lt 2 then begin 
    setup = strarr(2,ncol)
    setup[0,*] = lines2
  endif
  ; Two or more columns, only take the first two
  if ncol ge 2 then setup = lines2[0:1,*]

  count = npar

; No lines to process
endif else begin
  print,setupfiles[0],' HAS NOT LINES TO PROCESS'
  count = 0
  return
endelse

if keyword_set(stp) then stop

end
