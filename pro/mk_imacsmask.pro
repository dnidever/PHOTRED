pro mk_imacsmask,chip,mask,stp=stp

; Make a mask for IMACS images
; mask out vignetted areas.

undefine,mask

nchip = n_elements(chip)

if nchip eq 0 then begin
  print,'Syntax - mk_imacsmask,chip,mask,stp=stp'
  return
endif

nx = 2048.0
ny = 4096.0
xx = findgen(nx)#(fltarr(ny)+1.0)
yy = (fltarr(nx)+1.0)#findgen(ny)

; Remove regions that are vignetted
CASE chip of
1: coef = [ 2411.5237d0,  -2.1598841d0, 0.0010141577d0, -3.6877898d-07, 6.3852884d-11 ]
2: coef = [ 117.57659d0,  -0.52246527d0,   0.00019146285d0]
3: coef = [ -1353.9826d0,  0.99960654d0,  -0.00012923529d0]
4: coef = [ 180.94498d0, 0.63072300d0, 0.00016897444d0,  -9.6062773d-08, 7.1764188d-11]
5: coef = [ -391.74482d0, 0.026080168d0, 0.00011564051d0]
6: coef = [177.27057d0, 0.55727058d0, 0.00016613567d0, -1.4690285d-08, 3.4228032d-11]
7: coef = [ 2619.7328d0, -2.3845349d0, 0.0012565381d0, -4.9410640d-07, 8.7625961d-11]
8: coef = [ 170.84357d0, -0.51338861d0, 0.00013863439d0 ]
else: stop,'WRONG CHIP NUMBER ',strtrim(i,2)

ENDCASE

mask = float(yy gt poly(xx,coef))

if keyword_set(stp) then stop

end
