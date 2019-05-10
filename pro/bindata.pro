;+
;
; BINDATA
;
; Use this to bin X/Y data.  Uses Histogram.pro.  NANs in Y are automatically ignored.
;
; INPUTS:
;  x          Array of X values to be binned
;  y          Array of Y values to be binned
;  =binsize   Binsize for X dimension
;  =weights   Weights for the Y points. Only works for mean (the
;                default) and /std.
;  =min       Minimum X value   
;  =max       Maximum X value
;  /mean      Calculate the mean of Y for each bin (the default).
;  /med       Calculate the median of Y for each bin instead of the mean
;  /std       Return the standard deviation of Y for each bin instead
;             of the mean
;  /mad       Calculate a robust std.dev. using the Median Absolute Deviation
;  =perc      Calculate the Nth percentile.  This can also be used to
;               obtain the MINIMUM (perc=0) or MAXIMUM (perc=1)
;               values in each bin.
;  /hist      Histogram, just the number of points.
;  /sum       Sum of Y values.
;  /mode      Mode of Y values
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  xbin       The X-value for each bin (center).
;  ybin       The mean of the binned Y values per bin.  If /std is
;             set the standard deviation is used.  Empty bins that have no Y values will be NAN.
;  =gdind     The indices of ybin that are not NAN.
;  =xmnbin    The mean X-value for each bin.
;  =xstatbin  The X-value for each bin using the applied statistic.
;               Only for /MEAN and =PERC.
;
; USAGE:
;  IDL>bindata,x,y,xbin,ybin,binsize=5,min=-1,max=105
;
;
; By D.Nidever  April 2007
;-

pro bindata,x0,y0,xbin,ybin,binsize=binsize,weights=weights,$
            min=min,max=max,stp=stp,std=std,mean=meanset,med=med,mad=madset,$
            gdind=gdind,xmnbin=xmnbin,perc=perc,hist=hist,sum=sum,mode=mode,$  ;,nan=nan
            xstatbin=xstatbin
undefine,xbin,ybin,gdind,xmnbin,xstatbin

; Not enough inputs
if n_elements(x0) eq 0 or n_elements(y0) eq 0 then begin
  print,'Syntax - bindata,x,y,xbin,ybin,binsize=binsize,min=min,max=max,mean=mean,med=med,'
  print,'                 std=std,mad=mad,perc=perc,hist=hist,sum=sum,mode=mode,gdind=gdind,'
  print,'                 xmnbin=xmnbin,stp=stp'
  return
endif

; They need to be arrays
x = [x0]
y = [y0]

; Bin the X data and get the indices
;H = HISTOGRAM(X, binsize=binsize, min=min, max=max, nan=nan, REVERSE_INDICES = R, locations=locations)
H = HISTOGRAM(X, binsize=binsize, min=min, max=max, /nan, REVERSE_INDICES = R, locations=locations)

; How the REVERSE_INDICES WORK
;Set all elements of A that are in the ith bin of H to 0.  
;IF R[i] NE R[i+1] THEN A[R[R[I] : R[i+1]-1]] = 0  
;The above is usually more efficient than the following: 
;bini = WHERE(A EQ i, count)  
;IF count NE 0 THEN A[bini] = 0  

; What operation are we performing
oper = 1   ; Mean by default
if n_elements(meanset) gt 0 then oper=1  ; Mean, this is redundant
if n_elements(med) gt 0 then oper=2   ; Median
if n_elements(std) gt 0 then oper=3   ; Standard Deviation
if n_elements(madset) gt 0 then oper=4   ; MAD
if n_elements(perc) gt 0 then oper=5  ; Percentile
if n_elements(hist) gt 0 then oper=6  ; Histogram, i.e. number of points
if n_elements(sum) gt 0 then oper=7  ; Sum of points
if n_elements(mode) gt 0 then oper=8  ; Mode

; Check PERC
if oper eq 5 then begin
  if perc[0] lt 0.0 or perc[0] gt 1.0 then begin
    print,'PERC must be between 0.0 and 1.0'
    return
  endif
endif

nh = n_elements(h)
ybin = h*0.0+!values.f_nan   ; all bad by default
xmnbin = H*0.0
xstatbin = h*0.0+!values.f_nan   ; all bad by default  
if n_elements(weights) eq 0 then weights=x*0.0+1.0

; Loop through the bins
; and bin the Y values
for i=0,nh-1 do begin
  IF R[i] NE R[i+1] THEN begin
    x2 =  X[R[R[I] : R[i+1]-1]]
    y2 =  Y[R[R[I] : R[i+1]-1]]
    weights2 = WEIGHTS[R[R[I] : R[i+1]-1]]
    ny2 = n_elements(y2)
    gg = where(finite(y2) eq 1,ngdy2)
    ;ngdy2 = total(finite(y2))

    ; Seven operations: mean, median, std.dev, mad, percentile, histogram, sum
    CASE oper of
    ; Mean
    ;1: ybin[i] = TOTAL(y2,/nan)/ny2
    ;1: ybin[i] = TOTAL(y2,/nan)/ngdy2
    1: begin
         ybin[i] = TOTAL(y2[gg]*weights2[gg])/total(weights2[gg])
         xstatbin[i] = TOTAL(x2[gg]*weights2[gg])/total(weights2[gg])
       end
    ; Median
    2: ybin[i] = MEDIAN([y2],/even)
    ; Standard Deviation
    3: if ngdy2 eq 1 then ybin[i]=0.0 else begin
          ;ybin[i] = STDDEV(y2,/nan)
          ;weight = 1.0/sigmax^2.
          ; This is from WMEANERR.PRO
          sum    = total(weights[gg])
          ymean  = TOTAL(y2[gg]*weights2[gg])/sum
          ybin[i] = sqrt( total( ((y2-ymean)^2.)*weights2[gg] ) * ngdy2 / $
                         (( ngdy2-1.) * sum) )
       end
    ; MAD
    4: ybin[i] = MAD([y2])
    ; Percentile
    5:  begin
          gg = where(finite(y2) eq 1,ngg)
          if ngg gt 0 then begin
            si = sort(y2[gg])
            y3 = y2[gg[si]]
            x3 = x2[gg[si]]
            ybin[i] = y3[ round(ngg*perc[0]) < (ngg-1) ]
            xstatbin[i] = x3[ round(ngg*perc[0]) < (ngg-1) ]
          endif else begin
            ybin[i] = !values.f_nan
            xstatbin[i] = !values.f_nan             
          endelse
        end
    ; Histogram, number of points
    6:  ybin[i] = ny2
    ; Sum of points
    7:  ybin[i] = total(y2[gg])
    ; Mode
    8: begin
        ; Make a histogram and take the maximum point
        sig = MAD([y2])
        med = MEDIAN([y2])
        binsize2 = sig/5.0  ;10.0
        h2 = HISTOGRAM(y2[gg], binsize=binsize2, locations=xbin2,min=med-5.0*sig,max=med+5.0*sig)
        xbin2 = xbin2+0.5*binsize2
        maxind = first_el(maxloc(h2))
        ybin[i] = xbin2[maxind]
      end
    else: begin
       print,'Not a supported option'
       return
      end
    ENDCASE

    ; Mean x-values
    ;xmnbin[i] = TOTAL(x2,/nan)/ngdy2
    xmnbin[i] = TOTAL(x2*finite(y2),/nan)/ngdy2
    if finite(xmnbin[i]) eq 0 then xmnbin[i]=locations[i]+0.5*binsize

  endif
end

; X-values for each bin
xbin = locations+0.5*binsize

; Good indices
gdind = where(finite(ybin) eq 1,ngdind)

if keyword_set(stp) then stop

end
