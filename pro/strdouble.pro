function strdouble,num,nsig
; Return a high-precision string can be used
; to recreate the number in IDL
if n_elements(nsig) eq 0 then nsig=10
str = strtrim(string(num,format='(G'+strtrim(2*nsig>20,2)+'.'+strtrim(nsig,2)+')'),2)
if strpos(str,'E') ne -1 then str=repstr(str,'E','D') else str+='D'
return,str
end
