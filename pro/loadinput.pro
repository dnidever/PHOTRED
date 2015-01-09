pro loadinput,input0,list,comment=comment,count=count,exist=exist,stp=stp

;+
;
; NAME:
;  LOADINPUT
;
; PURPOSE:
;  This program can be used to load command-line inputs.
;  The input can be:
;   (1) an array list of files, i.e. ['one.txt','two.txt']
;   (2) a globbed list, i.e. "*.txt"
;   (3) a comma separated list, i.e.  'one.txt,two.txt'
;   (4) the name of a file that contains a list, i.e. "list.txt"
;         This can NOT be used in combination with any of the other three options.
;
; INPUTS:
;  input    The input list.  There are three possibilities
;            (1) an array list of files
;            (2) a globbed list, i.e. "*.txt"
;            (3) a comma separated list, i.e.  'one.txt,two.txt'
;            (4) the name of a file that contains a list, i.e. "@list.txt"
;                  This can NOT be used in combination with any of the other three options.
;  =comment Comment string to use when loading a list file. By default
;             comment='#'
;  /exist   Files must exist
;  /stp     Stop at the end of the program
;
; OUTPUTS:
;  list     The list of files
;  =count   The number of elements in list
;
; USAGE:
;  IDL>loadinput,'*.fits',list,count=count
;
; By D.Nidever   April 2007
;-

undefine,list
count = 0

; Not enough inputs
if n_elements(input0) eq 0 then begin
  print,'Syntax - loadinput,input,list,comment=comment,count=count,exist=exist,stp=stp'
  return
endif
input = input0

; Empty string input
if strtrim(input[0],2) eq '' then begin
  undefine,list
  count = 0
  return
endif


; A list file was input
if strmid(input[0],0,1) eq '@' then begin

  inp = strmid(input[0],1)

  ; Loading the files
  ;readcol,inp,list,format='A',/silent,comment='#'
  if n_elements(comment) eq 0 then comment='#'
  readline,inp,list,comment=comment,/noblank
  nlist = n_elements(list)

endif else begin

  ; Break up comma lists
  ninput = n_elements(input)
  undefine,temp
  for i=0.,ninput-1 do begin
    dum = strsplit(input[i],',',/extract)
    push,temp,dum
  end
  input = temp
  undefine,temp
  if n_elements(input) eq 1 then input=input[0]

  ; Probably an array of filenames
  if n_elements(input) gt 1 then begin

    ; Try to expand as wildcards
    wildind = where(strpos(input,'*') ne -1 or strpos(input,'?') ne -1,nwild)
    normind = where(strpos(input,'*') eq -1 and strpos(input,'?') eq -1,nnorm)
    
    ; Some wildcards
    if (nwild gt 0) then begin

      undefine,list
      if nnorm gt 0 then push,list,input[normind]
      wild = input[wildind]

      ; Loop through the wildcards
      for i=0.,nwild-1 do begin
        wfiles = file_search(wild[i],count=nwfiles)
        if nwfiles gt 0 then push,list,wfiles
      end
      
      nlist = n_elements(list)

    ; No wildcards
    endif else begin

      list = input
      nlist = n_elements(list)

    endelse


  ; A globbed list or only one file
  endif else begin
    list = file_search(input,count=nlist)
    if nlist eq 0 then list=input    ; Maybe just one filename
    ;nlist = n_elements(list)
  endelse

endelse

; The files must exist
if keyword_set(exist) then begin
  test = file_test(list)
  gd = where(test eq 1,ngd)
  if ngd gt 0 then begin
    list = list[gd]
    count = ngd
  endif else begin
    undefine,list
    count = 0
  endelse
endif

count = n_elements(list)
if count eq 1 then list=list[0]

if keyword_set(stp) then stop

end
