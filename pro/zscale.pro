pro zscale,im,z1,z2,contrast=contrast,nsample=nsample,stp=stp

;+
;
; ZSCALE.PRO
;
; PURPOSE:
;  This finds the IRAF display z1,z2 scalings
;  for an image
;
; INPUTS:
;  im         2D image
;  =contrast  Contrast to use (contrast=0.25 by default)
;  =nsample   The number of points of the image to use for the sample
;              (nsample=5e4 by default).
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  z1         The minimum value to use for display scaling.
;  z2         The maximum value to use for display scaling.
;
; USAGE:
;  IDL>zscale,im,z1,z2
;  IDL>display,im,min=z1,max=z2
;
;
;  From IRAF display help
;
;    If  the  contrast  is  not  zero  the  sample  pixels  are ranked in
;    brightness to form the function I(i) where i  is  the  rank  of  the
;    pixel  and  I is its value.  Generally the midpoint of this function
;    (the median) is very near the peak of the image histogram and  there
;    is  a  well defined slope about the midpoint which is related to the
;    width of the histogram.  At the ends of the I(i) function there  are
;    a  few very bright and dark pixels due to objects and defects in the
;    field.  To determine  the  slope  a  linear  function  is  fit  with
;    iterative rejection;
;    
;            I(i) = intercept + slope * (i - midpoint)
;    
;    If  more  than half of the points are rejected then there is no well
;    defined slope and the full range of the sample defines  z1  and  z2.
;    Otherwise  the  endpoints  of the linear function are used (provided
;    they are within the original range of the sample):
;    
;            z1 = I(midpoint) + (slope / contrast) * (1 - midpoint)
;            z2 = I(midpoint) + (slope / contrast) * (npoints - midpoint)
;
;  The actual IRAF program is called zsc_zlimits and is in the file
;  /net/astro/iraf2.12/iraf/pkg/images/tv/display/zscale.x
;
;
; By D.Nidever  Oct 2007  (using IRAF display zscale algorithm)    
;-


; No image input
if n_elements(im) eq 0 then begin
  print,'Syntax - zscale,im,z1,z2,contrast=contrast,stp=stp'
  return
endif

sz = size(im)
ndim = sz[0]
nx = sz[1]
ny = sz[2]
n = n_elements(im)

; Not an image
if (ndim ne 2) then begin
  print,'The input must be 2 dimensional'
  return
endif

if n_elements(contrast) eq 0 then contrast = 0.25

if n_elements(nsample) eq 0 then nsample = 5e4
nsample = (nsample < n)

xind = randomu(seed,nsample)*(nx-1L)+1L
yind = randomu(seed,nsample)*(ny-1L)+1L

f = im[xind,yind]
si = sort(f)
f2 = f[si]
x = findgen(nsample)
midpoint = round(nsample*0.5)
med = median(f)

zmin = min(im)
zmax = max(im)
zmed = median(im)

frac = 0.3 ;0.4
x3 = x[midpoint-frac*nsample:midpoint+frac*nsample]
f3 = f2[midpoint-frac*nsample:midpoint+frac*nsample]
;coef = poly_fit(x3,f3,1)
; Robust fitting program
coef = goodpoly(x,f2,1,2.0)
;print,coef

; y = m*x + b
; I = intercept + slope * (i-midpoint)
; I = intercept + slope * i - slope*midpoint
; I = slope * i + (intercept - slope*midpoint)
; b = intercept - slope*midpoint
; intercept = b + slope*midpoint

slope = coef[0,1]
intercept = coef[0,0] + slope*midpoint

;    z1 = I(midpoint) + (slope / contrast) * (1 - midpoint)
;    z2 = I(midpoint) + (slope / contrast) * (npoints - midpoint)
;z1 = f2[midpoint] + (slope/contrast) * (1L - midpoint)
;z2 = f2[midpoint] + (slope/contrast) * (nsample - midpoint)
z1 = zmed + (slope/contrast) * (1L - midpoint)
z2 = zmed + (slope/contrast) * (nsample - midpoint)

z1 = (z1 > zmin)
z2 = (z2 < zmax)

if keyword_set(stp) then stop

end
