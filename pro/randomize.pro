;+
;
;  This program randomly picks a given number of elements
;  out of a given array.  It returns the values picked.
;
; INPUTS:
;  inarray	The input array of values to use
;  num		How many numbers to output
;  =seed        The seed value to use.
;
; OUTPUTS:
;  outarray     The outputs array of randomly selected values
;  indx=indx    Indices of used values
;
; USAGE:
;  IDL>randomize,y,50,out,indx=indx
;
; by D.Nidever   June 2007
;-

pro randomize,inarray,num,outarray,indx=indx,seed=seed

; Not enough inputs
if n_elements(inarray) eq 0 or n_elements(num) eq 0 then begin
  print,'Syntax - randomize,inarray,num,outarray,indx=indx'
  return
endif 

indx=lonarr(num)-1
n=n_elements(inarray)
if n lt num then begin
  print,'Input array not large enough!'
  return
end
good=0

; Initialize seed
;  Use the first value of INARRAY to set seed
;  so it's always the same and reproducible for
;  the same data set.
if n_elements(seed) eq 0 then begin
  val = inarray[0]
  type = size(val,/type)
  ; Check if we have a number that we can use
  dum = where(type eq [1,2,3,4,5,6,9,12,13,14,15],isnum)
  if type eq 7 then isnum=valid_num(val,val)
  ; Number
  if isnum eq 1 then begin
    ; Convert all types to float
    val = fix(val,type=4)
    ; It's a finite number
    if finite(val) eq 1 then begin
       frac = abs(alog10(val))-floor(abs(alog10(val)))
       seed = round(frac*1000) > 1  ; number between 1 and 1000
    endif else seed=1
  ; Non-number
  endif else begin
    seed = 1
  endelse
endif

rnd = randomu(seed,n)   ; uniformly distributed random numbers
rint = sort(rnd)          ; the sorted indices will be random
indx = rint[0:num-1]

;; Now sort the indices
;si = sort(ind1)
;indx = ind1[si]

; Randomly picked values
outarray = inarray[indx]


; Old code
;for i=0,num-1 do begin
;
;  while good eq 0 do begin
;    ind=first_el(round((n-1)*randomu(seed,1)))		;random index
;    dum=where(indx eq ind,nbd)
;    if nbd eq 0 then begin
;      indx(i)=ind
;      good=1
;    end
;  end
;
;  good=0 		;restart the clock
;
;end	;for i
;
;;sorting final results
;si=sort(indx)
;indx=indx(si)
;
;outarray=inarray(indx)

;stop

end
