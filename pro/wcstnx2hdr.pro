; This puts the TNX WCS back into the header
;
; This does NOT put the non-standard CRVALOFFSET
; or ASTCD values into the header.

pro wcstnx2hdr, head, wcs

    ; Not enough inputs
    nhead = n_elements(head)
    nwcs = n_elements(wcs)
    if (nhead eq 0 or nwcs eq 0) then begin
      print,'Syntax - wcstnx2hdr, head, wcs'
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



    ; Put the TNX WAT parameters in the header

    ; WAT1
    ; Maximum length of string is 68 characters
    ctype1 = wcs.ast.ctype[0]
    axtype1 = first_el(strsplit(ctype1,'-',/extract))
    axtype1 = strlowcase(axtype1)

    ss1 = 'wtype=tnx axtype='+axtype1+' lngcor = "'
    ss1 = ss1+string(wcs.tnx1.fun_type,format='(F2.0)')
    ss1 = ss1+string(wcs.tnx1.xiorder,format='(F3.0)')
    ss1 = ss1+string(wcs.tnx1.etaorder,format='(F3.0)')
    ss1 = ss1+string(wcs.tnx1.cross_type,format='(F3.0)')
    ss1 = ss1+' '+strtrim(string(wcs.tnx1.ximin,format='(G25.16)'),2)
    ss1 = ss1+' '+strtrim(string(wcs.tnx1.ximax,format='(G25.16)'),2)
    ss1 = ss1+' '+strtrim(string(wcs.tnx1.etamin,format='(G25.16)'),2)
    ss1 = ss1+' '+strtrim(string(wcs.tnx1.etamax,format='(G25.16)'),2)
    for i=0,9 do $
      ss1 = ss1+' '+strtrim(string(wcs.tnx1.c[i],format='(G25.16)'),2)
    ss1 = ss1+' "'

    ; Add to header
    len1 = strlen(ss1)
    nwat1 = ceil(len1/68.)
    for i=0,nwat1-1 do begin
      num = strtrim(i+1,2)
      if i lt 10 then num='0'+num
      if i lt 100 then num='0'+num
      str = strmid(ss1,i*68,68)
      SXADDPAR,head,'WAT1_'+num,str
    end
 

    ; WAT2
    ctype2 = wcs.ast.ctype[1]
    axtype2 = first_el(strsplit(ctype2,'-',/extract))
    axtype2 = strlowcase(axtype2)

    ss2 = 'wtype=tnx axtype='+axtype2+' latcor = "'
    ss2 = ss2+string(wcs.tnx2.fun_type,format='(F2.0)')
    ss2 = ss2+string(wcs.tnx2.xiorder,format='(F3.0)')
    ss2 = ss2+string(wcs.tnx2.etaorder,format='(F3.0)')
    ss2 = ss2+string(wcs.tnx2.cross_type,format='(F3.0)')
    ss2 = ss2+' '+strtrim(string(wcs.tnx2.ximin,format='(G25.16)'),2)
    ss2 = ss2+' '+strtrim(string(wcs.tnx2.ximax,format='(G25.16)'),2)
    ss2 = ss2+' '+strtrim(string(wcs.tnx2.etamin,format='(G25.16)'),2)
    ss2 = ss2+' '+strtrim(string(wcs.tnx2.etamax,format='(G25.16)'),2)
    for i=0,9 do $
      ss2 = ss2+' '+strtrim(string(wcs.tnx2.c[i],format='(G25.16)'),2)
    ss2 = ss2+' "'

    ; Add to header
    len2 = strlen(ss2)
    nwat2 = ceil(len2/68.)
    for i=0,nwat2-1 do begin
      num = strtrim(i+1,2)
      if i lt 10 then num='0'+num
      if i lt 100 then num='0'+num
      str = strmid(ss2,i*68,68)
      SXADDPAR,head,'WAT2_'+num,str
    end

  if keyword_set(stp) then stop

end
