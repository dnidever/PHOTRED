;------------------------------------------------------------------------------
;+
; NAME:
;   djs_int2bin
;
; PURPOSE:
;   Convert integer number(s) to binary numbers.
;
; CALLING SEQUENCE:
;   binval = djs_int2bin(intval, [ndigit=ndigit])
;
; INPUTS:
;   intval:   Integer number(s)
;
; OPTIONAL INPUTS:
;   ndigit:   Number of binary digits in output; if not supplied, then the
;             minimum number of digits are used
;
; OUTPUTS:
;   binval:   Byte array(s) of binary values
;
; PROCEDURES CALLED:
;   djs_ceil()
;
; REVISION HISTORY:
;   Written D. Schlegel, 30 June 1997, Durham
;   31-Jul-1998  DJS - Subscripts modified to IDL 5 convention.
;-
;------------------------------------------------------------------------------
function djs_int2bin, intval, ndigit=ndigit
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - binval = djs_int2bin( intval, [ndigit=ndigit] )'
      return, -1
   endif

;   if (NOT keyword_set(ndigit)) then $
;    ndigit = djs_ceil( alog(max(intval)) / alog(2.0) )
   if (max(intval) LE 0) then maxdigit = 1 $
   else maxdigit = djs_ceil( alog(max(intval)) / alog(2.0) )

   mfac = 2L^indgen(maxdigit)

   ; Case where input value is a scalar
   if ((size(intval))[0] EQ 0) then begin
      nnum = 1
      binval = bytarr(maxdigit)
      for idig = maxdigit-1, 0, -1 do begin
         binval[idig] = intval / mfac[idig]
         intval = intval - binval[idig] * mfac[idig]
      endfor

   ; Case where input value is a vector
   endif else begin
      nnum = (size(intval))[1]
      binval = bytarr(maxdigit, nnum)
      for inum=0L, nnum-1 do begin
         for idig = maxdigit-1, 0, -1 do begin
            binval[idig,inum] = intval[inum] / mfac[idig]
            intval[inum] = intval[inum] - binval[idig,inum] * mfac[idig]
         endfor
      endfor
   endelse

   if (keyword_set(ndigit)) then begin
      if (ndigit LT maxdigit) then $
       binval = binval[0:ndigit-1,*]
      if (ndigit GT maxdigit) then $
       binval = [binval, bytarr(ndigit-maxdigit,nnum)]
   endif

   return, binval
end
;------------------------------------------------------------------------------
