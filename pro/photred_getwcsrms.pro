;+
;
; PHOTRED_GETWCSRMS
;
; Return the WCS RMS and NMATCH information from WCSFIT for the input files.
;
; INPUTS:
;  files   Array of file names.
;
; OUTPUTS:
;  str     Structure with WCS RMS and NMATCH information.
;
; USAGE:
;  IDL>str = photred_getwcsrms(files)
;
; By D. Nidever, Dec 2019
;-

function photred_getwcsrms,files

nfiles = n_elements(files)
if nfiles eq 0 then begin
  print,'Syntax - str=photred_getwcsrms(files)'
  return,-1
endif

str = replicate({file:'',rms:999999.0,nmatch:-1L},nfiles)
str.file = files
for i=0,nfiles-1 do begin
  ;;print,i,files[i]
  if file_test(files[i]) eq 1 then begin
    if strmid(files[i],6,7,/reverse_offset) eq 'fits.fz' then $
     head = PHOTRED_READFILE(files[i],exten=1,/header) else $
     head = PHOTRED_READFILE(files[i],exten=0,/header)
    rmsind = where(stregex(head,'WCSFIT: RMS',/boolean) eq 1,nrmsind)
    if nrmsind gt 0 then begin
      line = head[first_el(rmsind,/last)]
      dum1 = strsplitter(line,' ',/extract)
      dum2 = strsplit(dum1[2],'=',/extract)
      str[i].rms = float(dum2[1])
    endif
    nmatchind = where(stregex(head,'WCSFIT: NMATCH',/boolean) eq 1,nnmatchind)
    if nnmatchind gt 0 then begin
      line = head[first_el(nmatchind,/last)]
      dum1 = strsplitter(line,' ',/extract)
      dum2 = strsplit(dum1[2],'=',/extract)
      str[i].nmatch = float(dum2[1])
    endif
  endif
endfor

return,str

end
