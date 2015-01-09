function minloc,inarr,first=first,last=last
dum=where(inarr eq min(inarr,/nan),ndum) & mud = dum
;  dum=where(inarr eq min(inarr),ndum)  &  mud=dum
if keyword_set(first) then mud=dum(0)
if keyword_set(last) then mud=dum(ndum-1)
if keyword_set(first) and keyword_set(last) then begin
   mud=[dum(0),dum(ndum-1)]  &  if mud(0) eq mud(1) then mud=mud(0)
endif
return,mud
end ;minloc
