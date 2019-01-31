function create_index,arr,index

; Create an index of array values like reverse indices.
if n_elements(arr) eq 0 then begin
  print,'Syntax - index = create_index(arr)'
  return,-1
endif

narr = n_elements(arr)
si = sort(arr)
sarr = arr[si]
brklo = where(sarr ne shift(sarr,1),nbrk)
if nbrk gt 0 then begin
  brkhi = [brklo[1:nbrk-1]-1,narr-1]
  num = brkhi-brklo+1
  index = {index:si,value:sarr[brklo],num:num,lo:brklo,hi:brkhi}
endif else begin
  index = {index:si,value:arr[0],num:narr,lo:0L,hi:narr-1}
endelse
return,index

end
