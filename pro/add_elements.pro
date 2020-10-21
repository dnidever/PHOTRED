function add_elements,str,nnew,count=nstr
orig = str
norig = n_elements(orig)
str = make_structure(orig[0],norig+nnew)
str[0:norig-1] = orig
undefine,orig
nstr = n_elements(str)
return,str
end
