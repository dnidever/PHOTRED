;+
; NAME:
;	wmeanerr
; PURPOSE: (one line)
;	Calculate the mean and estimated error for a set of weighted data points
; DESCRIPTION:
;	This routine is adapted from Program 5-1, XFIT, from "Data Reduction
;	and Error Analysis for the Physical Sciences", p. 76, by Philip R.
;	Bevington, McGraw Hill.  This routine computes the weighted mean using
;	Instrumental weights (w=1/sigma^2).  The final uncertainty is
;	insensitive to a multiplicative constant on the weights.
; CATEGORY:
;	Statistics
; CALLING SEQUENCE:
;	wmeanerr,x,sigmax,xmean,xsigma
; INPUTS:
;	x      - Array of data points
;	sigmax - Array of errors in x
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;	None.
; OUTPUTS:
;	xmean   - weighted mean
;	xsigma  - uncertainty of xmean.  Weighted St.Dev.
;                   NOT Stdev. of the Mean.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
; MODIFICATION HISTORY:
;	Written by Marc W. Buie, Lowell Observatory, 2000 Oct 1
;  2000/10/3, MWB, mostly WMG, it now works.
;  Modified by D.Nidever
;-
pro wmeanerr,x,sigmax,xmean,xsigma

   undefine,xmean,xsigma

   ; Not enough inputs
   nx = n_elements(x)
   nsigmax = n_elements(sigmax)
   if (nx eq 0) then begin
     print,'Syntax - wmeanerr,x,sigmax,xmean,xsigma'
     return
   endif

   ; X and SIGMAX not the same size
   if nx ne nsigmax then begin
     print,'X and SIGMAX not the same size'
     return
   endif

   ; Calculate
   if n_elements(x) eq 1 then begin
      xmean  = x[0]
      xsigma = sigmax(0)
   endif else begin
      weight = 1.0/sigmax^2.
      sum    = total(weight)
      if sum eq 0.0 then begin
        print,'WMEANERR: sum is zero.  Assuming UNWEIGHTED'
        sum = 1.0
        weight = weight*0.+1.0
      endif
      sumx   = total(weight*x)
      xmean  = sumx/sum
      nx = float(n_elements(x))
      xsigma = sqrt( total( ((x-xmean)^2.)*weight ) * nx / $
                     (( nx-1.) * sum) )
   endelse

end
