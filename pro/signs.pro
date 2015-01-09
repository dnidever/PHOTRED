function signs,input

in = input
s=sign(in)
bad = where(s eq 0,nbad)
if nbad gt 0 then s(bad) = 1

return,s
end
