; extract the (somewhat strange) tnx WCS info from the extention

function hdr2wcstnx, h, $
                  crvaloffset_deg=crvaloffset_deg, $
                  astscale = astscale, astrota_deg = astrota_deg, $
                  astcd = astcd

    extast, h, ast 

    ss1 = sxpar(h,'WAT1_*')
    s1 = ''
    for i=0,n_elements(ss1)-1 do s1 = s1 + ss1[i]
    ;for i=0,4 do s1 = s1 + ss1[i]
    tnx1 = parsetnx(s1)
    ss2 = sxpar(h,'WAT2_*')
    s2 = ''
    for i=0,n_elements(ss2)-1 do s2 = s2 + ss2[i]
    ;for i=0,4 do s2 = s2 + ss2[i]
    tnx2 = parsetnx(s2)
    ccdsec  = fix(strsplit(sxpar(h,'CCDSEC'), '[:,]',/ex))
    datasec = fix(strsplit(sxpar(h,'CCDSEC'),'[:,]',/ex))
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

    return, {ast:ast, tnx1:tnx1, tnx2:tnx2, $
             ccdsec:ccdsec, datasec:datasec, $
             crvaloffset_deg:crvaloffset_deg, $
             astcd:astcd}

    ;return, {ast:ast, tnx1:tnx1, tnx2:tnx2, $
    ;         ccdsec:ccdsec, datasec:datasec, $
    ;         crvaloffset_deg:crvaloffset_deg, $
    ;         et:et, $
    ;         astcd:astcd}

end
