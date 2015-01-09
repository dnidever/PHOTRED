; This puts the TPV WCS back into the header
;
; This does NOT put the non-standard CRVALOFFSET
; or ASTCD values into the header.

pro wcstpv2hdr, head, wcs,maxpv=maxpv

    ; Not enough inputs
    nhead = n_elements(head)
    nwcs = n_elements(wcs)
    if (nhead eq 0 or nwcs eq 0) then begin
      print,'Syntax - wcstpv2hdr, head, wcs'
      return
    endif

    ; Put the "normal" WCS back in the header
    ast = wcs.ast
    ;PUTAST, head, ast

    ; Update all the standard WCS values
    ; This does NOT update LONGPOLE, LATPOLE, PV2_1, PV2_2
    SXADDPAR,head,'CTYPE1',ast.ctype[0]
    SXADDPAR,head,'CTYPE2',ast.ctype[1]
    SXADDPAR,head,'CRVAL1',ast.crval[0],format='G19.13'
    SXADDPAR,head,'CRVAL2',ast.crval[1],format='G19.13'
    SXADDPAR,head,'CRPIX1',double(ast.crpix[0]),format='G19.13'
    SXADDPAR,head,'CRPIX2',double(ast.crpix[1]),format='G19.13'
    SXADDPAR,head,'CD1_1',ast.cd[0,0],format='G19.13'
    SXADDPAR,head,'CD2_1',ast.cd[1,0],format='G19.13'
    SXADDPAR,head,'CD1_2',ast.cd[0,1],format='G19.13'
    SXADDPAR,head,'CD2_2',ast.cd[1,1],format='G19.13'

    ; Add PVi_X values
    if n_elements(maxpv) eq 0 then maxpv=39
    for i=0,maxpv do begin
      sxaddpar,head,'PV1_'+strtrim(i,2),wcs.pv1[i]
      sxaddpar,head,'PV2_'+strtrim(i,2),wcs.pv2[i]
    endfor

  if keyword_set(stp) then stop

end
