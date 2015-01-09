;+
;
; ALLFRAME_CALCWEIGHTS
;
;-

pro allframe_calcweights,mag,err,fwhm,rdnoise,medsky,actweight,scales

sz = size(mag)
ngd = sz[2]
nfiles = sz[1]

; Make 2D arrays for fwhm, rdnoise and medsky
fwhm2 = fwhm#(fltarr(ngd)+1.0)
rdnoise2 = rdnoise#(fltarr(ngd)+1.0)
medsky2 = medsky#(fltarr(ngd)+1.0)

; Computing intensity and weight for each star
;C Compute the intensity from the magnitude, easy enough
;            intensity=10**( (Mag(i,num(i))-25.0)/(-2.5))
;C Now compute the S/N (called the weight here)
;            weight(i,n)=(intensity/(Pi*FWHM(i)**2)) /
;     &      (sqrt(RDNOISE(i)**2
;     &      + MEDSKY(i) + intensity/(Pi*FWHM(i)**2) ) )
;C            print*, id(i,num(i))
; magnitude zero-point: 1 star ADU == 25.0 mag
intensity = 10.0^( (mag-25.0)/(-2.5) )
weight = (intensity/(!dpi*fwhm2^2.0))/ (sqrt(rdnoise2^2.0 + medsky2 + intensity/(!dpi*fwhm2^2.0) ) )
;weight(i,n)=(intensity/(Pi*FWHM(i)**2)) / (sqrt(RDNOISE(i)**2+ MEDSKY(i) + intensity/(Pi*FWHM(i)**2) ) )

; You get similar results if you use: weight = 1.0/err
; since magnitude errors are basically N/S

;C Now lets normalize the S/N (weight) for each star, take the maximimum S/N
;C and normalize so that is 1 S/N(i)/maximum S/N
;       do j=1,n-1
;          maxweight=-99999999.98
;          do i=1,FMAX
;           if(weight(i,j).gt.maxweight) maxweight=weight(i,j)
;          end do
;           weight(1:FMAX,j)=weight(1:FMAX,j)/maxweight
;C           print*,weight(1:FMAX,j),maxweight
;       end do

maxw = max(weight,dim=1)
maxw2 = (fltarr(nfiles)+1.0)#maxw
nweight = weight/maxw2

  
;C Finally computed the actual weight, by summing up the normalized weights,
;C and dividing thru by the sum
;       do j=1,FMAX
;        avgweight(j)=sum(weight(j,1:n))
;C/dble(n-1)
;C        print*, avgweight(j),n
;       end do
;        actweight(1:FMAX)=avgweight(1:FMAX)/sum(avgweight(1:FMAX))
;C Print them out, so we can put them somewhere useful
;       do j=1,FMAX
;          print*,actweight(j)
;       end do
;       stop
;       end

;avgweight = total(weight2,2)/ngd
avgweight = total(nweight,2)
actweight = avgweight/total(avgweight)
; They sum to 1.0


; Compute the scaling for each frame.  scale=1 for the first frame
scales = fltarr(nfiles)
for i=0,nfiles-1 do begin
  ratio = reform(intensity[i,*]/intensity[0,*])
  ratio_error = reform(sqrt(err[i,*]^2 + err[0,*]^2))  ; mag errors are ~fractional
  med_ratio = median(ratio)
  ;wmeanerr,ratio,ratio_error,xmean,xsigma
  scales[i] = med_ratio
end


;print,'Files: ',files
;print,'Weights: ',actweight
;print,'Scales: ',scales
;print,'Sky: ',medsky

;stop

if keyword_set(stp) then stop

end
