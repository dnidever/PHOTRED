function strtrim0,string,back=back,both=both

; gets rid of zeros at the front of a string

string = strtrim(string,2)
old_string = string
len = strlen(string)

flag = 1

if not keyword_set(back) then begin
  while flag eq 1 do begin
    len = strlen(string)
    dig = strmid(string,0,1)
    if dig eq '0' then string = strmid(string,1,len-1)
    if dig ne '0' then flag = 0
  end
endif

if keyword_set(back) or keyword_set(both) then begin
  len = strlen(string)

  flag = 1

  while flag eq 1 do begin
    len = strlen(string)
    dig = strmid(string,len-1,1)
    if dig eq '0' then string = strmid(string,0,len-1)
    if dig ne '0' then flag = 0
  end
endif

;stop

return,string

end
