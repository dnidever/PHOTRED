function strmult,string,n

; multiplies strings together

fstring = ''

if n eq 0 then fstring = ''

if n eq 1 then fstring = string

if n gt 1 then begin
  for i=0,n-1 do fstring = fstring+string
endif

return,fstring

end