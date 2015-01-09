pro wcstpvcor, xi, eta, wcs, xi1, eta1


    ; http://iraf.noao.edu/projects/ccdmosaic/tpv.html

    ; radius
    r = sqrt(xi^2 + eta^2)


     ;xi' = PV1_0 + PV1_1 * xi + PV1_2 * eta + PV1_3 * r +
     ;      PV1_4 * xi^2 + PV1_5 * xi * eta + PV1_6 * eta^2 +
     ;      PV1_7 * xi^3 + PV1_8 * xi^2 * eta + PV1_9 * xi * eta^2 + PV1_10 * eta^3 + PV1_11 * r^3 +
     ;      PV1_12 * xi^4 + PV1_13 * xi^3 * eta + PV1_14 * xi^2 * eta^2 + PV1_15 * xi * eta^3 + PV1_16 * eta^4 +
     ;      PV1_17 * xi^5 + PV1_18 * xi^4 * eta + PV1_19 * xi^3 * eta^2 +
     ;      PV1_20 * xi^2 * eta^3 + PV1_21 * xi * eta^4 + PV1_22 * eta^5 + PV1_23 * r^5 +
     ;      PV1_24 * xi^6 + PV1_25 * xi^5 * eta + PV1_26 * xi^4 * eta^2 + PV1_27 * xi^3 * eta^3 +
     ;      PV1_28 * xi^2 * eta^4 + PV1_29 * xi * eta^5 + PV1_30 * eta^6 +
     ;      PV1_31 * xi^7 + PV1_32 * xi^6 * eta + PV1_33 * xi^5 * eta^2 + PV1_34 * xi^4 * eta^3 +
     ;      PV1_35 * xi^3 * eta^4 + PV1_36 * xi^2 * eta^5 + PV1_37 * xi * eta^6 + PV1_38 * eta^7 + PV1_39 * r^7
     ;eta' = PV2_0 + PV2_1 * eta + PV2_2 * xi + PV2_3 * r +
     ;      PV2_4 * eta^2 + PV2_5 * eta * xi + PV2_6 * xi^2 +
     ;      PV2_7 * eta^3 + PV2_8 * eta^2 * xi + PV2_9 * eta * xi^2 + PV2_10 * xi^3 + PV2_11 * r^3 +
     ;      PV2_12 * eta^4 + PV2_13 * eta^3 * xi + PV2_14 * eta^2 * xi^2 + PV2_15 * eta * xi^3 + PV2_16 * xi^4 +
     ;      PV2_17 * eta^5 + PV2_18 * eta^4 * xi + PV2_19 * eta^3 * xi^2 +
     ;      PV2_20 * eta^2 * xi^3 + PV2_21 * eta * xi^4 + PV2_22 * xi^5 + PV2_23 * r^5 +
     ;      PV2_24 * eta^6 + PV2_25 * eta^5 * xi + PV2_26 * eta^4 * xi^2 + PV2_27 * eta^3 * xi^3 +
     ;      PV2_28 * eta^2 * xi^4 + PV2_29 * eta * xi^5 + PV2_30 * xi^6 +
     ;      PV2_31 * eta^7 + PV2_32 * eta^6 * xi + PV2_33 * eta^5 * xi^2 + PV2_34 * eta^4 * xi^3 +
     ;      PV2_35 * eta^3 * xi^4 + PV2_36 * eta^2 * xi^5 + PV2_37 * eta * xi^6 + PV2_38 * xi^7 + PV2_39 * r^7

     PV1 = wcs.pv1
     PV2 = wcs.pv2

     xi1 = PV1[0] + PV1[1] * xi + PV1[2] * eta + PV1[3] * r + $
           PV1[4] * xi^2 + PV1[5] * xi * eta + PV1[6] * eta^2 + $
           PV1[7] * xi^3 + PV1[8] * xi^2 * eta + PV1[9] * xi * eta^2 + PV1[10] * eta^3 + PV1[11] * r^3 + $
           PV1[12] * xi^4 + PV1[13] * xi^3 * eta + PV1[14] * xi^2 * eta^2 + PV1[15] * xi * eta^3 + PV1[16] * eta^4 + $
           PV1[17] * xi^5 + PV1[18] * xi^4 * eta + PV1[19] * xi^3 * eta^2 + $
           PV1[20] * xi^2 * eta^3 + PV1[21] * xi * eta^4 + PV1[22] * eta^5 + PV1[23] * r^5 + $
           PV1[24] * xi^6 + PV1[25] * xi^5 * eta + PV1[26] * xi^4 * eta^2 + PV1[27] * xi^3 * eta^3 + $
           PV1[28] * xi^2 * eta^4 + PV1[29] * xi * eta^5 + PV1[30] * eta^6 + $
           PV1[31] * xi^7 + PV1[32] * xi^6 * eta + PV1[33] * xi^5 * eta^2 + PV1[34] * xi^4 * eta^3 + $
           PV1[35] * xi^3 * eta^4 + PV1[36] * xi^2 * eta^5 + PV1[37] * xi * eta^6 + PV1[38] * eta^7 + PV1[39] * r^7

     eta1 = PV2[0] + PV2[1] * eta + PV2[2] * xi + PV2[3] * r + $
           PV2[4] * eta^2 + PV2[5] * eta * xi + PV2[6] * xi^2 + $
           PV2[7] * eta^3 + PV2[8] * eta^2 * xi + PV2[9] * eta * xi^2 + PV2[10] * xi^3 + PV2[11] * r^3 + $
           PV2[12] * eta^4 + PV2[13] * eta^3 * xi + PV2[14] * eta^2 * xi^2 + PV2[15] * eta * xi^3 + PV2[16] * xi^4 + $
           PV2[17] * eta^5 + PV2[18] * eta^4 * xi + PV2[19] * eta^3 * xi^2 + $
           PV2[20] * eta^2 * xi^3 + PV2[21] * eta * xi^4 + PV2[22] * xi^5 + PV2[23] * r^5 + $
           PV2[24] * eta^6 + PV2[25] * eta^5 * xi + PV2[26] * eta^4 * xi^2 + PV2[27] * eta^3 * xi^3 + $
           PV2[28] * eta^2 * xi^4 + PV2[29] * eta * xi^5 + PV2[30] * xi^6 + $
           PV2[31] * eta^7 + PV2[32] * eta^6 * xi + PV2[33] * eta^5 * xi^2 + PV2[34] * eta^4 * xi^3 + $
           PV2[35] * eta^3 * xi^4 + PV2[36] * eta^2 * xi^5 + PV2[37] * eta * xi^6 + PV2[38] * xi^7 + PV2[39] * r^7

end

