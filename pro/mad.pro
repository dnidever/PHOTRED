function mad,array,zero=zero,dim=dim

;+
;
; MAD - Median Absolute Deviation
; A robust method of variability.
; This program outputs 1.4826*MAD
; which is approximately equal to the
; Standard Deviation.
; 
; INPUTS:
;  array   An array to find the MAD for
;  /zero   Calculate the dispersion w.r.t to zero.
;            Use this if array is a vector of
;            residuals.  Similar to robust_sigma.pro
;  =dim    If array is 2D which dimension to
;            find the sigma over (starting with 1)
;
; OUTPUTS:
;  1.4826*MAD
;
; USAGE:
;  IDL>mad = mad(array)
;
; By D.Nidever  July 2006
;-

; Not enough inputs
if n_elements(array) eq 0 then begin
  print,'syntax - mad = mad(array,/zero,dim=dim)'
  return,-1
end

if n_elements(array) eq 1 then begin
  return,0.0
  ; return the single element value
  ;return,array[0]
end

if n_elements(dim) gt 0 and size(array,/n_dim) ne 2 then begin
  print,'ERROR: DIM can ONLY be used if ARRAY is 2D'
  return,-1
endif

;  This is what stat returns:
;  Median Absolute Deviation estimate of St.Dev.

sz = size(array)

; Array is greater than 1D, want sigma long one dimension
;   e.g. like TOTAL and MEDIAN
;----------------------------------------------------------
if size(array,/n_dim) eq 2 and n_elements(dim) gt 0 then begin

  ; Check that "dim" makes sense
  if dim lt 1 or dim gt size(array,/n_dim) then begin
    print,'ERROR: DIM must be between 1 and N_DIM(ARRAY)'
    return,-1
  endif

  if not keyword_set(zero) then begin

    med = median(array,dim=dim,/even)  ; 1D median
    if dim eq 1 then med2d = (fltarr(sz[1])+1.0)#med    else $  ; make 2D median
       med2d = med#(fltarr(sz[2])+1.0)

    mad = 1.4826*median(abs(array-med2d),dim=dim,/even)

  ; Dispersion w.r.t to zero
  endif else begin
    mad = 1.4826*median(abs(array),dim=dim,/even)
  endelse


; Array is 1D or want sigma of all elements
;--------------------------------------------
endif else begin

  if not keyword_set(zero) then begin
    mad = 1.4826*median(abs(array-median(array,/even)),/even)

  ; Dispersion w.r.t to zero
  endif else begin
    mad = 1.4826*median(abs(array),/even)
  endelse

endelse

return,mad

end
