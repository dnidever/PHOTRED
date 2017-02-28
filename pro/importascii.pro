;+
;
; IMPORTASCII.PRO
;
; This imports an ascii file by using the IDL read_ascii.pro
; routine.  An ASCII template has to be built on the fly.
; This is normally done interactively with ascii_template.pro
;
; INPUTS:
;  fname       Filename of the data file
;  comment     Comment string.  comment='#' by default.
;  /indef      Replace INDEFs with 999999
;  /allstr     Make all columns strings
;  /noprint    Don't print anything
;  /silent     Same as /noprint.
;  fieldnames  Array of field names.
;  fieldtypes  Array of the IDL field types
;  skipline    Number of lines you want to skip at the beginning
;  /allfloat   All columns are floats (and no INDEFS) and can
;              be read more quickly with READF than READ_ASCII
;  /header     The file has column headers on the first line.
;              Use these for the field names
;  =delimit    The delimiting character, space or tab used by default.
;                Note that there are problems getting the right
;                columns if there are lots of blank entries
;                (e.g. ",,,," with a delimiter of ",").  This can
;                fixed by adding a space for blank entries.
;  =dbformat   Get the field information from this IRAF CHART
;                dbformat file.
;  /preserve_null  Allow for null values separated by the delimiters.
;                    This is set by default for non-space delimiters.
;
; OUTPUTS:
;  arr         The data in an IDL structure
;  =count      The number of elements/lines read in.
;
; USAGE:
;  IDL>arr=importascii('data.txt',fieldnames=['name','ra','dec','mag','err'])
;
; Future improvements:
;  Check for end of file better.
;
; By David Nidever  Feb. 2006
;-

pro importascii_dummy
FORWARD_FUNCTION imp_datatypes
end

;-----------------------

Function imp_datatypes,arr

;+
; This function figures out what data types strings have
;
; INPUT
;  arr   Array of strings
; 
; OUPUT
;  typearr   Array of the data types (3-Long Integer, 4-Float, 5-Double, 7-String)
;            See the documentation for SIZE for a description of data types.
;
; David Nidever   Feb.2006
;-

; No Parmaters input
if n_params() eq 0 then begin
  print,'Syntax - typearr = imp_datatypes(array)'
  return,-1
endif


npar = n_elements(arr)
typearr = lonarr(npar)-1
validnum = valid_num(arr)
validint = valid_num(arr,/integer)
var = strtrim(arr,2)
nvar = strlen(var)

; Figuring out each column's data type
for i=0,npar-1 do begin
  ;var = strtrim(arr[i],2)
  ;bvar = byte(var)
  ;nvar = n_elements(bvar)

  ; String
  if validnum[i] eq 0 then type=7
  ; Float
  if validnum[i] eq 1 and validint[i] eq 0 then begin
    ; Float or Double?
    dec = first_el(strsplit(var[i],'.',/extract),/last)
    ndec = strlen(dec)
    ; What matters is the number of significant digits
    ;  not the decimal places
    ; subtract 1 to get rid of the decimal
    ndig = nvar[i]-1                   ; number of digits
    if strmid(var[i],0,1) eq '-' then ndig-=1   ; don't count the negative sign
    if ndig gt 7 then type=5 else type=4
  endif
  ; Long integer
  if validint[i] eq 1 and nvar[i] le 9 then type=3
  ; Long64 integer
  if validint[i] eq 1 and nvar[i] gt 9 then type=14
  ; NAN's are floats
  if strtrim(strupcase(var[i]),2) eq 'NAN' then type=4     ; float

  typearr[i] = type

;  ; If the string has only 0-9, -, +, ., E, e then it's a float, otherwise a string
;  ; 0-9 is 48-57 (in bytes)
;  ; "-" is 45
;  ; "+" is 43
;  ; "." is 46
;  ; "E" is 69
;  ; "e" is 101
;  ; "D" is 68
;  ; "d" is 100
;  bfloat = [bindgen(10)+48B,43B,45B,46B,69B,101B,68B,100B]
;  bint = [bindgen(10)+48B,43B,45B]
;
;  ; Checking each character in "bvar"
;  ; Are there any non-"number" characters?
;  badfloat = 0           ; float until proven otherwise
;  badint = 0
;  last = 0B
;  for j=0,nvar-1 do begin
;    ; Checking for float characters
;    g = where(bvar[j] eq bfloat,ng)  ; is this character a float characters
;    if ng eq 0 then badfloat=1
;
;    ; Checking for integer characters
;    g = where(bvar[j] eq bint,ng)    ; is this character an integer characters
;    if ng eq 0 then badint=1
;
;    ; Checking for plus or minus, must be at beginning, okay after 'E' or 'e'
;    if (bvar[j] eq 43B) or (bvar[j] eq 45B) then begin
;         if (j ne 0) and (last ne 69B) and (last ne 101B) then badfloat=1
;         badint = 1
;    endif
;
;    ; Checking for period, CAN'T be at beginning
;    if (bvar[j] eq 46B) then begin
;      if (j eq 0) then badfloat=1
;      badint = 1
;    endif
;
;    last = bvar[j]
;  endfor
;
;  ; Use VALID_NUM.PRO as a double-check
;  if valid_num(var) eq 1 then badfloat=0
;  if valid_num(var,/integer) eq 1 then badint=0
;
;  if valid_num(var) eq 0 then badfloat=1
;  if valid_num(var,/integer) eq 0 then badint=1
;
;  ; String
;  if (badfloat eq 1) then type = 7   ; String
;
;  ; Float
;  if (badfloat eq 0 and badint eq 1) then begin
;
;    ; Float or Double?
;    dec = first_el(strsplit(var,'.',/extract),/last)
;    ndec = strlen(dec)
;
;    ;; type = 5, Double
;    ;; type = 4, Float
;    ;if (ndec ge 6) then type=5 else type=4
;
;    ; What matters is the number of significant digits
;    ;  not the decimal places
;    ; subtract 1 to get rid of the decimal
;    ndig = nvar-1                   ; number of digits
;    if var[0] eq '-' then ndig-=1   ; don't count the negative sign
;    if ndig gt 7 then type=5 else type=4
;    ;if (nvar-1) gt 6 then type=5 else type=4
;
;  endif
;
;  ; Long Integer
;  if (badfloat eq 0 and badint eq 0) then type = 3   ; Integer (Long)
;
;  ; Long64 integer
;  if (badfloat eq 0 and badint eq 0 and nvar gt 9) then type = 14   ; Long64
;
;  ; NAN's are floats
;  if strtrim(strupcase(var),2) eq 'NAN' then type = 4     ; float
;
;  typearr[i] = type

endfor

;stop

return,typearr

end


;---------------------------------------------------------------------


Function importascii,fname,indef=indef,allstr=allstr,$
         comment=comment,noprint=noprint,fieldnames=fieldnames0,$
         allfloat=allfloat,skipline=skipline,stp=stp,fieldtypes=fieldtypes,$
         header=header,delimit=delimit,dbformat=dbformat,silent=silent,count=count,$
         preserve_null=preserve_null

; No parameters input
if n_params() eq 0 then begin
  print,'Syntax - arr=importascii(filename,indef=indef,allstr=allstr,comment=comment,noprint=noprint,'
  print,'                         fieldnames=fieldnames,fieldtypes=fieldtypes,allfloat=allfloat,'
  print,'                         skipline=skipline,stp=stp,header=header,delimit=delimit,'
  print,'                         dbformat=dbformat,count=count,preserve_null=preserve_null'
  return,-1
endif

; Skip header line
if keyword_set(header) and not keyword_set(skipline) then skipline=1

if n_elements(comment) eq 0 then comment='#'   ; Comment string
if n_elements(indef) eq 0 then indef=1         ; By default remove indef
if n_elements(allstr) eq 0 then allstr=0       ; By default use proper types
if n_elements(skipline) eq 0 then data_start=0 else data_start=skipline  ; number of lines to skip
if n_elements(silent) gt 0 then noprint=silent
if n_elements(indef) eq 0 then indef=1         ; replace indef's and null's by default

if keyword_set(fieldtypes) then typearr = fieldtypes
if keyword_set(fieldnames0) then fieldnames = fieldnames0

; Getting field information from dbformat file
if keyword_set(dbformat) then begin

  test = file_test(dbformat)
  if test eq 1 then begin

    ; Loading file
    READLINE,dbformat,lines,/noblank
    delim = string(32B)+string(9B)   ; space and tab
    dum = strsplitter(lines,delim,/extract)

    fieldnames = reform(dum[0,*])
    ncol = n_elements(fieldnames)
    type = strtrim(reform(dum[1,*]),2)
    typearr = lonarr(ncol)+7
    gstr = where(type eq 'char',ngstr)
    if ngstr gt 0 then typearr[gstr] = 7
    gint = where(type eq 'int',ngint)
    if ngint gt 0 then typearr[gint] = 3
    greal = where(type eq 'real',ngreal)
    if ngreal gt 0 then typearr[greal] = 4
    gdouble = where(type eq 'double',ngdouble)
    if ngdouble gt 0 then typearr[gdouble] = 5

  endif else begin
    print,dbformat,' NOT FOUND'
  endelse

endif


; Searching for the file
d = file_search(fname)
if (fname[0] eq '') or (d[0] eq '') then begin
  print,'FILE ',fname,' NOT FOUND!!'
  return,-1
endif

; FIGURING out how many columns there are
str1=''
openr,unit,fname,/get_lun

; Skipping lines
if (data_start gt 0) then begin
  dum=''
  for i=0,data_start-1 do readf,unit,dum
endif

; Making sure the line IS NOT commented out
readf,unit,str1
first = strmid(str1,0,1)
while((first eq comment) and (eof(unit) eq 0)) do begin
  readf,unit,str1
  first = strmid(str1,0,1)
endwhile

close,unit
free_lun,unit

; Default delimiter character
if not keyword_set(delimit) then begin
  delimit = 32B  ; space is the default delimiter
 
  ; Split the string
  ;  SPACE, by default do NOT use preserve_null
  if keyword_set(preserve_null) then arr1 = strsplit(str1,' ',/extract,/preserve_null) else $
    arr1 = strsplit(str1,' ',/extract) 
  ;  TAB, by default USE preserve_null
  if n_elements(preserve_null) gt 0 and not keyword_set(preserve_null) then $
    arr2 = strsplit(str1,string(9B),/extract) else $
    arr2 = strsplit(str1,string(9B),/extract,/preserve_null)
 
  ; Use tab as delimiter
  if n_elements(arr2) gt n_elements(arr1) then begin
    delimit = 9B
    arr1 = arr2
    if n_elements(preserve_null) eq 0 then preserve_null=1  ; use preserve_null
  endif

; Using INPUT delimiter
endif else begin
  if delimit ne ' ' and n_elements(preserve_null) eq 0 then preserve_null=1
  arr1 = strsplit(str1,delimit,/extract,/preserve_null)
endelse
ncol = n_elements(arr1)  ; # of columns

; Empty last element
trimlast = 0
if ncol eq n_elements(fieldnames)+1 and keyword_set(preserve_null) and arr1[n_elements(arr1)-1] eq '' then begin
  arr1 = arr1[0:n_elements(arr1)-2]
  ncol = n_elements(arr1)
  trimlast = 1
endif

; Checking for ALL FLOATS
typearr = imp_datatypes(arr1)
bd = where(typearr eq 7,nbd)
if (nbd eq 0) then begin
  if not keyword_set(noprint) then print,'ALL FLOATS?'
endif

; Using the header line
if keyword_set(header) then begin
  headstr=''
  openr,unit,fname,/get_lun
  readf,unit,headstr
  close,unit
  free_lun,unit

  ; Remove comment strings
  headstr = stress(headstr,'D',0,'#')
  headstr = stress(headstr,'D',0,comment)

  ; Split the string
  headarr = strsplit(headstr,delimit,/extract)
  headarr = strtrim(headarr,2)  ; remove blank spaces
  nhead = n_elements(headarr)

  ; Not the same, try both delimiters
  if nhead ne ncol then begin
    headarr = strsplit(headstr,string(9B)+string(32B),/extract)
    nhead = n_elements(headarr)
  endif

  ; Picking fieldnames
  if nhead ge ncol then begin
    fieldnames = headarr[0:(nhead-1) < (ncol-1)]
  endif else begin
    print,'Headers do not match # of columns.  Ncol(header)=',strtrim(nhead,2),$
          ' Ncol(file)=',strtrim(ncol,2)
    return,-1
  endelse

  ; Check that the field names are valid idl tagnames
  nfieldnames = n_elements(fieldnames)
  orig_fieldnames = fieldnames
  for i=0,nfieldnames-1 do begin
    ; Converting to valid IDL tag name
    fieldnames[i] = idl_validname(fieldnames[i],/convert_all)
    ; Messaging the change
    if not keyword_set(noprint) and fieldnames[i] ne orig_fieldnames[i] then $
      print,'Converting to valid IDL tag name: ',orig_fieldnames[i],' -> ',fieldnames[i]
  endfor

  ;stop

endif


; Field Names
if not keyword_set(fieldnames) then $
  fieldnames = 'FIELD'+strtrim(sindgen(ncol),2)
if n_elements(fieldnames) ne ncol then begin
  print,'Input array of field names not of the correct size'
  print,'Input field names = ',strtrim(n_elements(fieldnames),2)
  print,'File columns = ',strtrim(ncol,2)
  fieldnames = 'FIELD'+strtrim(sindgen(ncol),2)
  return,-1
  ;dum=''
  ;read,dum,'Do you want to continue?'
  ;if strlowcase(strmid(strtrim(dum,2),0,1)) ne 'y' then stop
  ;stop
endif 


; Using READ_ASCII to read the data
IF not keyword_set(allfloat) THEN BEGIN

  START:

  ; All strings at the beginning
  typearr = lonarr(ncol)+7

  ; Figuring out field locations.
  fieldlocations = strsplit(str1,string(delimit),preserve_null=preserve_null)   ; this returns the locations
  if keyword_set(trimlast) then fieldlocations=fieldlocations[0:n_elements(fieldlocations)-2] ; trim last

  fieldgroups = lindgen(ncol)

  ; Creating the template
  template = {version:1.0, datastart:0, delimiter:delimit, missingvalue:!values.f_nan,$
              commentsymbol:'#',fieldcount:ncol,fieldtypes:typearr,fieldnames:fieldnames,$
              fieldlocations:fieldlocations,fieldgroups:fieldgroups}

  ; Reading the Data
  arr = READ_ASCII(fname,template=template,data_start=data_start)

  nrow = n_elements(arr.(0))

  ;; Making it into a "normal" structure
  ;cmd = 'dum = {'
  ;for i=0,ncol-1 do begin
  ;  if typearr(i) eq 3 then char = '0L'
  ;  if typearr(i) eq 4 then char = '0.0'
  ;  if typearr(i) eq 5 then char = '0.d0'
  ;  if typearr(i) eq 7 then char = '""'
  ;  if typearr(i) eq 14 then char = '0LL'
  ;  cmd=cmd+fieldnames(i)+':'+char
  ;  if i ne ncol-1 then cmd = cmd+', '
  ;end
  ;cmd = cmd+'}'
  ;
  ;ddd = execute(cmd)
  ;str = replicate(dum,nrow)

  dum = CREATE_STRUCT(fieldnames[0],fix('',type=typearr[0]))
  for i=1,ncol-1 do $
    dum = CREATE_STRUCT(dum,fieldnames[i],fix('',type=typearr[i]))
  str = replicate(dum,nrow)

  ; Transferring the data
  ;STRUCT_ASSIGN does NOT work here
  ; str is a N dimensional structure
  ; arr is a 1 dimensional structure with large arrays
  for i=0,ncol-1 do str.(i) = arr.(i)
  undefine,arr

  ; Converting INDEF's to 999999
  ; Converting INDEF's and null to NAN
  if keyword_set(indef) then begin

    nbad = 0

    ; Looping through the fields
    for i=0,ncol-1 do begin
      col = strtrim(reform(str.(i)),2)
      ;bd = where(col eq 'INDEF' or col eq "'INDEF'",nbd)
      ;if nbd gt 0 then str(bd).(i) = '999999'
      bd = where(col eq 'INDEF' or col eq "'INDEF'" or col eq 'null',nbd)      
      if nbd gt 0 then str[bd].(i) = 'NAN'
      nbad = nbad + nbd
    endfor

    if (not keyword_set(noprint)) then $
    ;  if (nbad gt 0) then print,"INDEF's converted to 999999"
      if (nbad gt 0) then print,"INDEF's and null's converted to NAN"

  endif

  ; Converting them to the proper types
  if not keyword_set(allstr) then begin

    ; Getting the data types
    if not keyword_set(fieldtypes) then begin
      for i=0,ncol-1 do begin

        ; Check the first 100 (or nrow) for the type
        ; Use the maximum type
        ; The higher types encompass all othe lower ones.
        ; 7-string > 5-double > 4-float > 3-long > 2-int
        ;type = datatypes(str(0).(i))
        type100 = imp_datatypes(str[0:99<(nrow-1)].(i))    ; using first 100 rows 
        type = max(type100)                            ; use the maximum 
        typearr[i] = type
      endfor
    endif else begin
      typearr = fieldtypes
    endelse


    ;; Making a new structure with the proper types
    ;cmd = 'dum = {'
    ;for i=0,ncol-1 do begin
    ;  if typearr(i) eq 3 then char = '0L'
    ;  if typearr(i) eq 4 then char = '0.0'
    ;  if typearr(i) eq 5 then char = '0.d0'
    ;  if typearr(i) eq 7 then char = '""'
    ;  if typearr(i) eq 14 then char = '0LL'
    ;  cmd=cmd+fieldnames(i)+':'+char
    ;  if i ne ncol-1 then cmd = cmd+', '
    ;end
    ;cmd = cmd+'}'
    ;
    ;ddd = execute(cmd)
    ;arr = replicate(dum,nrow)

    dum = CREATE_STRUCT(fieldnames[0],fix('',type=typearr[0]))
    for i=1,ncol-1 do $
      dum = CREATE_STRUCT(dum,fieldnames[i],fix('',type=typearr[i]))
    arr = replicate(dum,nrow)

    ; Check for missing data in floating point columns here
    ;   harder to do that above with the INDEFs because
    ;   we might not know the data types yet.
    nbad2 = 0
    for i=0,ncol-1 do begin
      if typearr[i] eq 4 or typearr[i] eq 5 then begin
        col = strtrim(reform(str.(i)),2)  ; STR is all strings
        bd = where(col eq '',nbd)
        if nbd gt 0 then str[bd].(i) = 'NAN'
        nbad2 += nbd
      endif
    endfor
    if not keyword_set(noprint) and nbad2 gt 0 then print,'Missing data converted to NAN'

    ; Transferring the data
    STRUCT_ASSIGN,str,arr          ; source, destination
    ;for i=0,ncol-1 do arr.(i) = str.(i)
    ; Rename
    str = arr
    undefine,arr
  endif


; /ALLFLOAT, Using READF to read the data
ENDIF ELSE BEGIN

  ; Checking the data types
  typearr = imp_datatypes(arr1)
  bd = where(typearr eq 7,nbd)

  ; Some of the columns are strings
  if nbd gt 0 then begin
    print,'Some of the columns are strings'
    print,'Using READ_ASCII to read the data'
    goto,START
  endif

  nrow = file_lines(fname)
  if keyword_set(header) then nrow=nrow-1
  arr = fltarr(ncol,nrow)

  ; Reading the data
  openr,unit,fname,/get_lun

  ; Skipping lines
  if (data_start gt 0) then begin
    dum=''
    for i=0,data_start-1 do readf,unit,dum
  endif

  readf,unit,arr
  close,unit
  free_lun,unit

  ;; Creating the normal structure
  ;cmd = 'dum = {'
  ;for i=0,ncol-1 do begin
  ;  cmd=cmd+fieldnames(i)+':0.0'
  ;  if i ne ncol-1 then cmd = cmd+', '
  ;end
  ;cmd = cmd+'}'
  ;
  ;ddd = execute(cmd)
  ;str = replicate(dum,nrow)

  dum = CREATE_STRUCT(fieldnames[0],0.0)
  for i=1,ncol-1 do $
    dum = CREATE_STRUCT(dum,fieldnames[i],0.0)
  str = replicate(dum,nrow)

  ; Transferring the data
  for i=0,ncol-1 do str.(i) = reform(arr(i,*))

ENDELSE


; Report on the data
if not keyword_set(noprint) then $
  print,strtrim(ncol,2),' columns x ',strtrim(nrow,2),' rows'


if keyword_set(stp) then stop

count = n_elements(str)

return,str

end
