; extract the tpv WCS info from the extention

function hdr2wcstpv, h, $
                  crvaloffset_deg=crvaloffset_deg, $
                  astscale = astscale, astrota_deg = astrota_deg, $
                  astcd = astcd

    extast, h, ast 

    ; http://iraf.noao.edu/projects/ccdmosaic/tpv.html

    ; Get all of the PV1_X and PV2_X values from header
    pv1 = fltarr(40)  ; PV1_0 to PV1_39
    pv1[1] = 1
    pv2 = fltarr(40)  ; PV2_0 to PV2_39
    pv2[1] = 1
    ; all parameters default to 0 except PV1_1 and PV2_1
    ; which default to 1
    for i=0,39 do begin
      pv1num = sxpar(h,'PV1_'+strtrim(i,2),count=npv1num,/silent)
      if npv1num gt 0 then pv1[i]=pv1num
      pv2num = sxpar(h,'PV2_'+strtrim(i,2),count=npv2num,/silent)
      if npv2num gt 0 then pv2[i]=pv2num
    endfor

    ;ss1 = sxpar(h,'WAT1_*')
    ;s1 = ''
    ;for i=0,n_elements(ss1)-1 do s1 = s1 + ss1[i]
    ;;for i=0,4 do s1 = s1 + ss1[i]
    ;tnx1 = parsetnx(s1)
    ;ss2 = sxpar(h,'WAT2_*')
    ;s2 = ''
    ;for i=0,n_elements(ss2)-1 do s2 = s2 + ss2[i]
    ;;for i=0,4 do s2 = s2 + ss2[i]
    ;tnx2 = parsetnx(s2)
    ccdsec  = fix(strsplit(sxpar(h,'CCDSEC',/silent), '[:,]',/ex))
    if n_elements(ccdsec) eq 1 and ccdsec[0] eq 0 then $
      ccdsec = [1,sxpar(h,'naxis1'),1,sxpar(h,'naxis2')]
    datasec = fix(strsplit(sxpar(h,'CCDSEC',/silent),'[:,]',/ex))
    ;datasec = fix(strsplit(sxpar(h,'DATASEC'),'[:,]',/ex))
    if n_elements(datasec) eq 1 and datasec[0] eq 0 then $
      datasec = [1,sxpar(h,'naxis1'),1,sxpar(h,'naxis2')]
    ;datasec = fix(strsplit(sxpar(h,'DATASEC'),'[:,]',/ex))
    ;et = utc2et(sxpar(h,'DATE-OBS'))

    ; astcd - transforms nominal XY->xi,eta (coordinate xi,eta)
    ; to astrometric xi,eta, thus incorporating a camera
    ; rotation and other effects such as stellar aberration
    if not keyword_set(astcd) then begin
        if not keyword_set(astscale) then astscale = 1.d
        if not keyword_set(astrota_deg) then astrota_deg = 0.d
        deg = !dpi / 180.d
        sc = astscale * cos(astrota_deg * deg)
        ss = astscale * sin(astrota_deg * deg)
        astcd = astscale * [[sc,-ss],[ss,sc] ]
    endif

    ; crvaloffset - crpix_true = crpix_head + crvaloffset
    ; accounts for "zeroing the offset's"
    if not keyword_set(crvaloffset_deg) then crvaloffset_deg = [0.d, 0.d] 

    return, {ast:ast, pv1:pv1, pv2:pv2, $
             ccdsec:ccdsec, datasec:datasec, $
             crvaloffset_deg:crvaloffset_deg, $
             astcd:astcd}

    ;return, {ast:ast, tnx1:tnx1, tnx2:tnx2, $
    ;         ccdsec:ccdsec, datasec:datasec, $
    ;         crvaloffset_deg:crvaloffset_deg, $
    ;         et:et, $
    ;         astcd:astcd}

end
