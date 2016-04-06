;+
;
; PRINTSTR
;
;  This program prints an IDL structure to an ASCII file
;
; INPUTS:
;  str         The structure to output
;  file        The name of the output ASCII file to write to
;  /noheader   Don't add a header line
;  =dbformat   Create a dbformat file
;  /silent     Don't print anything to the screen.
;  /stp        Stop at the end
;
; OUTPUTS:
;  The structure is printed to an ASCII file.
;
; USAGE:
;  IDL>printstr,str,'filename.txt',noheader=noheader,dbformat=dbformat
;               silent=silent,stp=stp
;
;  Written by D.Nidever  May 2006
;-

pro printstr_dummy
FORWARD_FUNCTION formatize
end

;------------------

function formatize,type

; This function gets the correct format code for a given data type

if not keyword_set(type) then typ=-1 else typ=type

case typ of
0:  form = ''         ; Undefined
1:  form = 'A4'       ; Byte
2:  form = 'I9'       ; Int
3:  form = 'I12'      ; Long
;3:  form = 'I13'      ; Long
;4:  form = 'G15.7'    ; Float
;5:  form = 'G20.13'   ; Double
4:  form = 'F13.4'    ; Float
5:  form = 'F15.6'    ; Double
6:  form = 'G15.7'    ; Complex
7:  form = 'A22'      ; String
8:  form = ''         ; Struct
9:  form = 'G20.13'   ; DComplex
10: form = ''         ; Pointer
11: form = ''         ; Objref
12: form = 'I9'       ; Uint
13: form = 'I13'      ; Ulong
14: form = 'I22'      ; Long64
15: form = 'I22'      ; Ulong64
else: form = ''
endcase

return,form

end

;---------------------------------------------------------

pro printstr,str,file,noheader=noheader,dbformat=dbformat,stp=stp,silent=silent

; This program prints a structure to an ASCII file

if not keyword_set(file) then file='idl.txt'

if not keyword_set(str) then begin
  print,'Syntax - printstr,str,file,noheader=noheader,dbformat=dbformat,'
  print,'                  silent=silent,stp=stp'
  return
endif


tags = tag_names(str)
ntags = n_elements(tags)

for i=0,ntags-1 do begin

  sz = size(str(0).(i))
  nsz = n_elements(sz)
  type = sz(nsz-2)
  dim = sz(0)
 
  ;if i ne ntags-1 then add = ',' else add=''
  add = ','

  ; Figuring out what the format is and the names
  case dim of
    0: begin
         push,names,tags(i)
         ;format = format+formatize(type)+add
         push,format,'('+formatize(type)+')'
       end
    1: begin
         n = sz(1)
         for j=0,n-1 do begin
           push,names,tags(i)+strtrim(j,2)
           ;format = format+formatize(type)+add
           push,format,'('+formatize(type)+')'
         end
       end

    2: begin
         n1 = sz(1)
         n2 = sz(2)
         for j=0,n1-1 do begin
           for k=0,n2-1 do begin
             push,names,tags(i)+strtrim(j,2)+strtrim(k,2)
             ;format = format+formatize(type)+add
             push,format,'('+formatize(type)+')'
           end  ; for k
         end  ; for j
       end ; 2
  endcase

end ; for i


openw,unit,/get_lun,file

; Printing the header
if not keyword_set(noheader) then begin
  head = '# '

  if keyword_set(dbformat) then head = head+'     KEY     '

  nnames = n_elements(names)
  for i=0,nnames-1 do begin
    len = strlen(names[i])                     ; length of name
    ind = strpos(format[i],'.')                ; position of . in format string

    ; Floating point
    if ind[0] ne -1 then begin
      spaces = long(strmid(format[i],2,ind-2))   ; 
    endif

    ; Integer or Character
    if ind[0] eq -1 then begin
      len2 = strlen(format[i])
      spaces = long(strmid(format[i],2,len2-3))      
    endif

    ; Spaces before and after the name
    nfirst = floor(float(spaces-len)*0.5) > 1
    nlast = ( spaces - (len+nfirst) ) > 1
    first = string(replicate(32B,nfirst))
    last = string(replicate(32B,nlast))
    head = head+first+names[i]+last

  end

  printf,unit,head

endif

nrec = n_elements(str)


; Are there any tags that have are arrays?
normal =1
for i=0,ntags-1 do begin
  sz = size(str(0).(i))
  nsz = n_elements(sz)
  type = sz(nsz-2)
  dim = sz(0)
  if dim gt 0 then normal=0
end


; Normal array
IF keyword_set(normal) then begin

  ; Start 2D output structure
  nstr = n_elements(str)
  outarr = strarr(ntags,nrec)
  format = strarr(ntags)

  ; Make a string array and put in tag by tag
  ; formatting each.
  for i=0,ntags-1 do begin
    type = size(str[0].(i),/type)
    fmt = formatize(type)

    ; String
    if type eq 7 then begin
      lenarr = strlen(str.(i))
      maxlen = max(lenarr)
      fmt = 'A'+strtrim(maxlen+2,2)
    endif

    ; Float
    if type eq 4 then begin
      ;maximum = long64(max(str.(i))+1.0)
      maximum = long64(max(abs(str.(i)))+1.0)
      maxlen = strlen(strtrim(maximum,2))
      totlen = (maxlen+5+2) > 11
      fmt = 'F'+strtrim(totlen,2)+'.4'
    endif

    ; Double
    if type eq 5 then begin
      ;maximum = long64(max(str.(i))+1.0)
      maximum = long64(max(abs(str.(i)))+1.0)
      maxlen = strlen(strtrim(maximum,2))
      totlen = (maxlen+7+2) > 11
      fmt = 'F'+strtrim(totlen,2)+'.6'
    endif

    format[i] = fmt
    outarr[i,*] = string(str.(i),format='('+fmt+')')
  endfor 

  ; Make 1D output structure
  outarr1 = strarr(nrec)
  for i=0.,nrec-1 do outarr1[i]=strjoin(reform(outarr[*,i]),'') 

  ; Write it to the file
  for i=0.,nrec-1 do printf,unit,reform(outarr1[i])

  close,unit
  free_lun,unit


; Some tags are arrays
ENDIF else begin


  ; Looping through records
  for i=0LL,nrec-1 do begin

    cnt = 0
    text = ''

    if keyword_set(dbformat) then begin
      text = text + string( i+1, FORM='('+formatize(3)+')')
    endif

    ; Looping through tags
    for j=0,ntags-1 do begin

      sz = size(str(0).(j))
      nsz = n_elements(sz)
      type = sz(nsz-2)
      dim = sz(0)
 
      case dim of
        0: begin
             text = text + string( str[i].(j), FORM=format[cnt++])
           end
        1: begin
             n = sz(1)
             for k=0,n-1 do begin
               text = text + string( (str[i].(j))(k), FORM=format[cnt++])        
             end ; for k
           end

        2: begin
             n1 = sz(1)
             n2 = sz(2)
             for k=0,n1-1 do begin
                for l=0,n2-1 do begin
                  text = text + string( (str[i].(j))(k,l), FORM=format[cnt++])
               end  ; for l
             end  ; for k
           end  ; 2
      endcase

    end ; for j, tag loop

    printf,unit,text

  end  ; for i

ENDELSE

close,unit
free_lun,unit

if not keyword_set(silent) then $
  print,'STRUCTURE written to ASCII file: ',file

; DBFORMAT file
if keyword_set(dbformat) then begin
  if (string(dbformat) eq '1' or size(dbformat,/type) ne 7) then $
    dbfile = 'idl.dbformat' else dbfile=dbformat

  openw,unit,/get_lun,dbfile

  printf,unit,'KEY            int     %13d   key'

  nnames = n_elements(names)
  for i=0,nnames-1 do begin
    form = format(i)
    name = names(i)
    len = strlen(name)
    space1 = string(replicate(32B,15-len))
    space2 = string(replicate(32B,4))
    lo = strpos(form,'.')
    length = long(strmid(form,2,lo-2))

    ; Float or Double
    ;if strmid(form,1,1) eq 'G' or strmid(form,1,1) eq 'F' then begin
    if strmid(form,0,1) eq 'G' or strmid(form,0,1) eq 'F' then begin
      ;hi = strpos(form,')')
      typ1 = 'real'
      ;typ2 = '%'+strmid(form,2,hi-2)+'f'
      typ2 = '%'+strmid(form,1)+'f'
    endif
    ; Integer
    ;if strmid(form,1,1) eq 'I' then begin
    if strmid(form,0,1) eq 'I' then begin
      len = strlen(form)
      ;length = long(strmid(form,2,len-3))      
      length = long(strmid(form,1,len-1))      
      typ1 = 'int'
      typ2 = '%'+strtrim(length,2)+'d'
    endif
    ; Character
    ;if strmid(form,1,1) eq 'A' then begin
    if strmid(form,0,1) eq 'A' then begin
      len = strlen(form)
      ;length = long(strmid(form,2,len-3))
      length = long(strmid(form,1,len-1))
      typ1 = 'char'
      typ2 = '%'+strtrim(length,2)+'s  '+strtrim(length,2)
    endif
    printf,unit,name,space1,typ1,space2,typ2

  end

  close,unit
  free_lun,unit

  if not keyword_set(silent) then $
    print,'DBFORMAT file written to: ',dbfile

  ; Creating keys file
  openw,unit,/get_lun,'idl.keys'
  printf,unit,'KEY '+strtrim(names(0),2)
  printf,unit,strtrim(names(0),2)+' '+strtrim(names(1),2)
  printf,unit,strtrim(names(1),2)+' '+strtrim(names(2),2)
  close,unit
  free_lun,unit

  print,'KEYS file written to: idl.keys'  

end

if keyword_set(stp) then stop

end

