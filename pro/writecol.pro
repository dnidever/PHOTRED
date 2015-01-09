;+ 
; NAME:
; writecol
;  Version 1.1
;
; PURPOSE:
;    Prints a series of arrays to a file in ASCII format
;
; CALLING SEQUENCE:
;   
;   writecol, filename, v1, v2, FMT=''
;
; INPUTS:
;   file     - Name of the ASCII file
;   v1       - Vector 1
;   v2       - Vector 2
;   [v3-v40]       - Vectors 3-40
;
; RETURNS:
;
; OUTPUTS:
;   Prints v1, v2 to screen
;
; OPTIONAL KEYWORDS:
;   FMT -  FORTRAN formatting
;   FILNUM - File number (as opposed to file)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The program keys off the number of elements in v1
;
; EXAMPLES:
;   writecol, 'arrays.dat', array1, array2
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-June-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro writecol, file, v1, v2, v3, v4, v5, v6, v7, v8, v9, $
              v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, $
              v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, $
              v31, v32, v33, v34, v35, v36, v37, v38, v39, v40, $
              FMT=fmt, FILNUM=filnum


; writecol -- Writes a 2 column ascii file

  if (N_params() LT 2) or n_elements(file) eq 0 then begin 
    print,'Syntax - ' + $
             'writecol, file, v1, v2, [v3-v19] FMT=, FILNUM= '
    return
  endif 

;

  flgvn = N_params()-1
  if not keyword_set( FMT ) then    flgfmt    = 0 else begin
      flgfmt = 1 
      fmt = fmt[0]
  endelse

;

  ; Print to screen
  if file eq '' or strtrim(file,2) eq '-1' then begin
    filnum = -1
  endif 

  if not keyword_set(FILNUM) then begin
      filnum = 91
      close, filnum
      openw, filnum, file
      flg_fil = 91
  endif

  for i=0LL,n_elements(v1)-1 do begin
      case flgvn of 

          40: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i],v33[i],v34[i],v35[i],v36[i],v37[i],v38[i],v39[i],v40[i]
          39: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i],v33[i],v34[i],v35[i],v36[i],v37[i],v38[i],v39[i]
          38: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i],v33[i],v34[i],v35[i],v36[i],v37[i],v38[i]
          37: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i],v33[i],v34[i],v35[i],v36[i],v37[i]
          36: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i],v33[i],v34[i],v35[i],v36[i]
          35: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i],v33[i],v34[i],v35[i]
          34: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i],v33[i],v34[i]
          33: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i],v33[i]
          32: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i],v32[i]
          31: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],$
            v30[i],v31[i]
          30: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i],v30[i]
          29: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i],v29[i]
          28: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i],v28[i]
          27: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i],v27[i]
          26: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i],v26[i]
          25: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i],v25[i]
          24: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i],v24[i]
          23: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i],v23[i]
          22: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i],v22[i]
          21: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i],$
            v19[i],v20[i],v21[i]
          20: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
            v16[i],v17[i], v18[i], v19[i],v20[i]
          19: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
            v16[i],v17[i], v18[i], v19[i]
          18: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
            v16[i],v17[i], v18[i]
          17: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i]
          16: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i]
          15: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i]
          14: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i]
          13: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i]
          12: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i]
          11: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i]
          10: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i]
          9: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i]
          8: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i]
          7: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i]
          6: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i]
          5: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i], v4[i], v5[i]
          4: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i], v4[i] 
          3: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i]
          2: printf, filnum, FORMAT=fmt, v1[i], v2[i]
          1: printf, filnum, FORMAT=fmt, v1[i]
          else: stop
      endcase
  endfor
  if keyword_set(FLG_FIL) then close, filnum


return
end
