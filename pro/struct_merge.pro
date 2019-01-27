;+
;
; STRUCT_MERGE
;
; This program merged structures with the same number of elements.
;
; INPUTS:
;  str1   The first structure to merge.
;  str2   The second structure to merge.
;  str3   The third structure to merge (optional).
;
; OUTPUTS:
;  final   The final merged structure.
;
; USAGE:
;  IDL>final = struct_merge(str1,str2,str3)
;
; By D. Nidever  Sep2018
;-

function struct_merge,str1,str2,str3

  ;; Merge structures

  nstr1 = n_elements(str1)
  nstr2 = n_elements(str2)
  nstr3 = n_elements(str3)  
  ;; Not enough inputs
  if nstr1 eq 0 or nstr2 eq 0 then begin
    print,'Syntax - struct_merge,str1,str2,str3'
    return,-1
  endif
  if nstr1 ne nstr2 then begin
    print,'The number of elements in STR1 and STR2 are not the same'
    return,-1
  endif
  if nstr3 gt 0 and nstr3 ne nstr1 then begin
    print,'The number of elements in STR1 and STR3 are not the same'
    return,-1
  endif
  if size(str1,/type) ne 8 or size(str2,/type) ne 8 then begin
    print,'Inputs must be structures'
    return,-1
  endif

  tags1 = tag_names(str1)
  tags2 = tag_names(str2)
  tags3 = tag_names(str3)  
  ntags1 = n_tags(str1)
  ntags2 = n_tags(str2)
  ntags3 = n_tags(str3)
  tags = [tags1,tags2,tags3]
  ntags = n_elements(tags)
  ui = uniq(tags,sort(tags))
  ;if n_elements(ui) ne ntags then begin
  ;  print,'Some columns are duplicated.  Adding unique characters'
  ;  stop 
  ;endif
  
  ;; Get the schema
  schema = str1[0]
  struct_assign,{dumdum:''},schema
  schema2 = str2[0]
  struct_assign,{dumdum:''},schema2
  for i=0,ntags2-1 do begin
     colname = tags2[i]
     ;; Duplicate column name
     cnt = 1
     while total(tag_names(schema) eq colname) gt 0 do begin
       colname = tags2[i]+strtrim(cnt,2)
       cnt++
     endwhile
     schema = create_struct(schema,colname,schema2.(i))
  endfor
  ;schema = create_struct(schema,schema2)
  if nstr3 gt 0 then begin
    schema3 = str3[0]
    struct_assign,{dumdum:''},schema3
    for i=0,ntags3-1 do begin
      colname = tags3[i]
      ;; Duplicate column name
      cnt = 1
      while total(tag_names(schema) eq colname) gt 0 do begin
        colname = tags3[i]+strtrim(cnt,2)
        cnt++
      endwhile
      schema = create_struct(schema,colname,schema3.(i))
    endfor
    ;schema = create_struct(schema,schema3)
  endif
  final = replicate(schema,nstr1)

  ;; Start loading in the data
  for i=0,ntags1-1 do final.(i)=str1.(i)
  for i=0,ntags2-1 do final.(i+ntags1)=str2.(i)
  for i=0,ntags3-1 do final.(i+ntags1+ntags2)=str3.(i)  

  return,final
  
  end
