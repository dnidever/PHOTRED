;this program returns the indices of maximum points and minimum points
;it also plots the array and max/min points if desires.

;function sign,dum
;  return,dum/abs(dum)
;end

pro dln_maxmin,array,minarr,maxarr,plot=plot,win=win

;returns maxima and minima points

n=n_elements(array)

if n eq 0 then begin
  print,'Syntax - dln_maxmin,array,minarr,maxarr,plot=plot,win=win'
  return
endif

minarr=-1L
maxarr=-1L

; Only one element
if n eq 1 then return  

;finding first derivative
deriv=array(1:n-1)-array(0:n-2)

nel=n_elements(deriv)

for i=0.,nel-2 do begin
;  if deriv(0) eq 0.0 then begin		;zeros
;    if i eq 0 then begin		;at the beginning
;      if deriv(i+1) gt 0.0 then minarr=[minarr,i+1]
;      if deriv(i+1) lt 0.0 then maxarr=[maxarr,i+1]
;    end
;			
;    if i gt 0 then begin		;in the middle
;      lft=deriv(i-1)
;      rt=deriv(i+1)
;
;      if lft gt 0.0 and rt lt 0.0 then maxarr=[maxarr,i+1]
;      if lft lt 0.0 and rt gt 0.0 then minarr=[minarr,i+1]
;    end
;    stop
;  end
  if deriv[i] gt 0.0 and deriv[i+1] le 0.0 then maxarr=[maxarr,i+1]	;maximum
  if deriv[i] lt 0.0 and deriv[i+1] ge 0.0 then minarr=[minarr,i+1]	;minimum
end

nmin = n_elements(minarr)
if nmin gt 1 then minarr=minarr[1:nmin-1]
nmax = n_elements(maxarr)
if nmax gt 1 then maxarr=maxarr[1:nmax-1]

;deriv1=deriv/abs(deriv)
;n1=n_elements(deriv1)
;deriv2=deriv1(1:n1-1)-deriv1(0:n1-2)
;minarr=where(deriv2 eq 2)+1
;maxarr=where(deriv2 eq -2)+1

if keyword_set(plot) then begin
  if keyword_set(win) then begin
    co1=200 & co2=200
  endif
  if not keyword_set(win) then begin
    ;co1=300 & co2=800
    co1=100 & co2=100
  endif

  x=findgen(n)
  plot,array,co=co1
  if nmax gt 1 and nmin gt 1 then begin
    oplot,x[minarr],array[minarr],ps=1,co=co2
    oplot,x[maxarr],array[maxarr],ps=2,co=co2
  end
end

end


