function first_el,expre,nth_el=nth_el,last=last
;
if keyword_set(nth_el) then nth_el=fix(abs(nth_el)) else nth_el=0
if keyword_set(last) then nth_el=n_elements(expre)-1
return,expre(nth_el)
end
