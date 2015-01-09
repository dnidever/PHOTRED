; given x and y (in IDL coordinates - 1st pixel centered on 0)
; calculate ra and dec in radians, unless /degree set

pro wcstnx_xy2rd, xin, yin, wcs, ra, dec, degree=degree

    ; decompose the wcs
    CRPIX1 = wcs.ast.crpix[0]
    CRPIX2 = wcs.ast.crpix[1]
    CD1_1 = wcs.ast.cd[0,0]
    CD1_2 = wcs.ast.cd[0,1]
    CD2_1 = wcs.ast.cd[1,0]
    CD2_2 = wcs.ast.cd[1,1]
    CCDSEC1 = wcs.ccdsec[0]
    CCDSEC2 = wcs.ccdsec[2]
    DATASEC1 = wcs.datasec[0]
    DATASEC2 = wcs.datasec[2]
    deg = !dpi/180.d


    ; convert to FITS indexing of original data coordinates
    x = xin + 1 + datasec1 - ccdsec1
    y = yin + 1 + datasec2 - ccdsec2


    ; Compute the first order standard
    ; coordinates xi and eta from the
    ; linear part of the solution stored
    ; in CRPIX and the CD matrix
    xi  = CD1_1 * (x - CRPIX1) + CD1_2 * (y - CRPIX2)
    eta = CD2_1 * (x - CRPIX1) + CD2_2 * (y - CRPIX2)

    ;Add the non-linear part of the projection using the coefficients 
    ;in the WAT keywords as described below
    xip  =  xi + wcstnxcor (xi, eta, wcs.tnx1)
    etap = eta + wcstnxcor (xi, eta, wcs.tnx2)
    ;xip = xi
    ;etap = eta


    ;Apply additional tranformation to go from nominal
    ;XY->xi,eta mapping (to account for rotation of the camera,
    ;changes in focus, or even stellar aberration
    xiast  = (wcs.astcd[0,0] * xip + wcs.astcd[0,1] * etap) * deg
    etaast = (wcs.astcd[1,0] * xip + wcs.astcd[1,1] * etap) * deg


    ;Apply the standard tangent plane projection to xi' amd eta' 
    ;using the CRVAL values as the tangent point to get the 
    ;RA and DEC in degrees. 
    A = (wcs.ast.crval[0] + wcs.crvaloffset_deg[0])*deg
    D = (wcs.ast.crval[1] + wcs.crvaloffset_deg[1])*deg
    xieta2rd, A, D, xiast, etaast, ra, dec

    ;; Convert to hours and degrees
    ;ra = ra*!radeg/15.0
    ;dec = dec*!radeg

    ; Convert RA/DEC to degrees
    if keyword_set(degree) then begin
      ra = ra / deg
      dec = dec / deg
    endif

;    print, 'wcstnx_xy2rd: x', x, fo='(A20,2F15.9)'
;    print, 'wcstnx_xy2rd: xi', xi, fo='(A20,2F15.9)'
;    print, 'wcstnx_xy2rd: xip', xip, fo='(A20,2F15.9)'
;    print, 'wcstnx_xy2rd: xiast', xiast, fo='(A20,2F15.9)'
;    print, 'wcstnx_xy2rd: ra', ra, fo='(A20,2F15.9)'

end

