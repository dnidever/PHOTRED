;+
; NAME:
;  wcstnx_rd2xy
; PURPOSE: (one line)
;  World Coordinate System transformation
; DESCRIPTION:
; given ra and dec in radians
; calculate x and y (in IDL coordinates - 1st pixel centered on 0)
; CATEGORY:
;  Astrometry
; CALLING SEQUENCE:
;  wcstnx_rd2xy, ra, dec, wcs, xout, yout
; INPUTS:
;  ra - right ascention in radians
;  dec - declination in radians
;  wcs - a structure containing astrometric info (see hdr2wcs)
; OPTIONAL INPUT PARAMETERS:
;  /degree - RA/DEC inputs are in degrees
; KEYWORD INPUT PARAMETERS:
; KEYWORD OUTPUT PARAMETERS:
;  xiast, etaast - 'astrometric' xi, eta
;  xip, etap - the non-linear xi,eta
; OUTPUTS:
;  x, y - IDL coordinates - 1st pixel centered on 0
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;  Written 2006 Feb, by Leslie Young, SwRI
;-


pro wcstnx_rd2xy, ra0, dec0, wcs, xout, yout, degree=degree, $
                  xip=xip, etap=etap, xia=xia, etaa=etaa
    ; decompose the wcs
    CRPIX1 = wcs.ast.crpix[0]
    CRPIX2 = wcs.ast.crpix[1]
    invcd = invert(wcs.ast.cd)
    INVCD1_1 = invcd[0,0]
    INVCD1_2 = invcd[0,1]
    INVCD2_1 = invcd[1,0]
    INVCD2_2 = invcd[1,1]
    CCDSEC1 = wcs.ccdsec[0]
    CCDSEC2 = wcs.ccdsec[2]
    DATASEC1 = wcs.datasec[0]
    DATASEC2 = wcs.datasec[2]
    deg = !dpi/180.d

    ra = ra0
    dec = dec0

    ; Convert RA/DEC to radians if /degree set
    if keyword_set(degree) then begin
      ra = ra * deg
      dec = dec *deg
    endif

    ;Apply the standard tangent plane projection to xi' amd eta' 
    ;using the CRVAL values as the tangent point to get the 
    ;RA and DEC in degrees. 
    A = (wcs.ast.crval[0] + wcs.crvaloffset_deg[0])*deg
    D = (wcs.ast.crval[1] + wcs.crvaloffset_deg[1])*deg
    rd2xieta, A, D, ra, dec, xiast, etaast

    ;Apply additional tranformation to go from 
    ;astrometric xi,eta
    ;nominal XY->xi,eta mapping 
    ;(to account for rotation of the camera,
    ;changes in focus, or even stellar aberration
    invastcd = invert(wcs.astcd)
    xip  = (invastcd[0,0] * xiast + invastcd[0,1] * etaast) / deg
    etap = (invastcd[1,0] * xiast + invastcd[1,1] * etaast) / deg

    ;Subtract the non-linear part of the projection using the coefficients 
    ;in the WAT keywords as described below
    xi  = xip
    eta = etap
    for i=0,3 do begin
        dxi  =  wcstnxcor (xi, eta, wcs.tnx1)
        deta =  wcstnxcor (xi, eta, wcs.tnx2)
        xi  =  xip - dxi
        eta = etap - deta
    end
    
    ; Compute the first order standard
    ; coordinates xi and eta from the
    ; linear part of the solution stored
    ; in CRPIX and the INVCD matrix
    x  = INVCD1_1 * xi + INVCD1_2 * eta + CRPIX1
    y  = INVCD2_1 * xi + INVCD2_2 * eta + CRPIX2


    ; convert from FITS indexing of original data coordinates
    xout = x - (1 + datasec1 - ccdsec1)
    yout = y - (1 + datasec2 - ccdsec2)

    if not isarray(ra) then begin
        xout = xout[0]
        yout = yout[0]
    end

;    print, 'wcstnx_rd2xy: x', x, fo='(A20,2F15.9)'
;    print, 'wcstnx_rd2xy: xi', xi, fo='(A20,2F15.9)'
;    print, 'wcstnx_rd2xy: xip', xip, fo='(A20,2F15.9)'
;    print, 'wcstnx_rd2xy: xiast', xiast, fo='(A20,2F15.9)'
;    print, 'wcstnx_rd2xy: ra', ra, fo='(A20,2F15.9)'

end

;pro wcstnx_rd2xyTEST
;
;    ; cal046_CCD7C313.2(1937,105)
;    ; cal046_CCD7Pluto(1945,149.2)
;    ; cal046_CCD7SOAR std(1912,152)
;    
;    xin = [1937, 1912]
;    yin = [ 105,  152]
;    d = rd_ccd7(46, hprime, hext)
;    wcs = hdr2wcs(hext,crvalof=[-0.01102,-0.15312])
;    wcstnx_xy2rd, xin, yin, wcs, ra, dec
;    wcstnx_rd2xy, ra,dec, wcs, xout, yout
;    print, xout, yout
;

;;  WHAT WE EXPECT
;; UCAC2
;;C313.2    17:28:55.0174 -15:00:54.749   
;;comp      17:28:55.8893 -15:01:01.395
;
;
;    ; cal130_CCD7C313.2/Pluto(1969,94)
;    ; cal130_CCD7SOAR Std.(1941,148)
;
;
;end
