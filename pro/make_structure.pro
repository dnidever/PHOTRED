function make_structure,str1,num,count=count
  
schema = str1[0]
struct_assign,{dum:''},schema
count = num
return,replicate(schema,num)

end
