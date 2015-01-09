function array_indices2,im,ind

; This function is the opposite of IDL's array_indices.pro that converts
; 1D array indices to 2D array indices.

nim = n_elements(im)
nind = n_elements(ind)

if nim eq 0 or nind eq 0 then begin
  print,'Syntax - ind = array_indices2(im,ind2d)
  return,-1
end

b = where(ind eq -1,nb)
if nb gt 0 then begin
  print,'Indices cannot be negative!'
  return,-1
end

index = ind

sz = size(im)
dim0 = sz(1)

sz2 = size(ind)

; only one point input
if (sz2(0) eq 1 and sz2(1) eq sz(0)) then begin
  index = fltarr(sz(0),1)
  index(*,0) = ind
endif

if (sz2(1) ne sz(0)) then begin
  print,'Index must be a '+strtrim(sz(0),2)+'D array'
  return,-1
endif

; If im is NxM then:
; 0=[0,0], 1=[1,0], 2=[2,0], ..., N+5=[5,1], 2N+6=[6,2], etc.
index2 = reform(index(0,*))
for i=1,sz(0)-1 do begin
  mult = 1L
  for j=1,i do mult = mult*sz(j)
  index2=index2+reform(index(i,*))*mult
end

return,index2

end