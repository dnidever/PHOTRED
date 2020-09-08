;+
;
; LOADTFR
;
; This loads a DAOMATCH/DAOMASTER tfr file.
;
; INPUTS:
;  tfrfile  The TFR filename
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  files    The list of files in the TFR file
;  str      The structure with final ID, X, Y and
;             the index array.  The indices are
;             FORTRAN/IRAF 1-index based.
;
; USAGE:
;  IDL>loadtfr,'ccd1001.tr',files,str
;
; By D.Nidever   August 2016
;-

pro loadtfr,tfrfile,files,str,stp=stp

undefine,files
undefine,str

; Enough inputs
ntfrfile = n_elements(tfrfile)
if ntfrfile eq 0 then begin
  print,'Syntax - loadtfr,tfrfile,files,str,stp=stp'
  return
endif

; Test the file
test = file_test(tfrfile)
if test eq 0 then begin
  print,file,' NOT FOUND'
  return
endif

; A TFR file looks like this:
; F1-00423034_01.als              99.9999   9.9999
; F1-00423033_01.als              99.9999   9.9999
; F1-00423035_01.als              99.9999   9.9999
; F1-00423036_01.als              99.9999   9.9999
; F1-00423037_01.als              99.9999   9.9999
; F1-00423038_01.als              99.9999   9.9999
; F1-00423039_01.als              99.9999   9.9999
; F1-00423040_01.als              99.9999   9.9999
; F1-00423041_01.als              99.9999   9.9999
; F1-00423042_01.als              99.9999   9.9999
; F1-00423043_01.als              99.9999   9.9999
; F1-00423044_01.als              99.9999   9.9999
; F1-00423045_01.als              99.9999   9.9999
; F1-00423046_01.als              99.9999   9.9999
; F1-00423047_01.als              99.9999   9.9999
; F1-00423048_01.als              99.9999   9.9999
; ==============================
;      1 -1037.98 2452.949      0      0      0      0      0      0    340      0      0      0      0      0      0      0      0      0
;      2 -1037.67 3505.380      0      0      0      0      0      0   1222      0      0      0      0      0      0      0      0      0
;      3 -1036.54 3174.594      0      0      0      0      0      0    448      0      0      0      0      0      0      0      0      0
;      4 -1035.85 4321.116      0      0      0      0      0      0   1263      0      0      0      0      0      0      0      0      0
;      5 -1035.28 5458.115      0      0      0      0      0      0    729      0      0      0      0      0      0      0      0      0
;      6 -1033.22 2134.540      0      0      0      0      0      0    838      0      0      0      0      0      0      0      0      0
;      7 -1032.40 3823.881      0      0      0      0      0      0   1126      0      0      0      0      0      0      0      0      0
;      8 -1031.18 5946.214      0      0      0      0      0      0   1075      0      0      0      0      0      0      0      0      0
;      9 -1030.97 5823.931      0      0      0      0      0      0      0      0   1773      0      0      0      0      0      0      0
;     10 -1030.16 5403.574      0      0      0      0      0      0    725      0   2157      0      0      0      0      0      0      0
;     11 -1029.83 4989.178      0      0      0      0      0      0      0      0   2110      0      0      0      0      0      0      0
;     12 -1029.31 5322.905      0      0      0      0      0      0      0      0    700      0      0      0      0      0      0      0
;     13 -1029.17 3798.451      0      0      0      0      0      0      0      0    377      0      0      0      0      0      0      0

; Read in the information
READLINE,tfrfile,tfrlines

; Find the break in in the list
brkind = where(stregex(tfrlines,'====',/boolean) eq 1,nbrkind)
if nbrkind eq 0 then begin
  print,'ERROR. No ===== break line'
  return
endif

; Filenames
filelines = tfrlines[0:brkind[0]-1]
files = reform( (strsplitter(filelines,' ',/extract))[0,*] )
nfiles = n_elements(files)

; Transformation info
tlines = tfrlines[brkind[0]+1:*]
ntlines = n_elements(tlines)
tarr = strsplitter(tlines,' ',/extract)
allid = long(reform(tarr[0,*]))
allx = float(reform(tarr[1,*]))
ally = float(reform(tarr[2,*]))
index = tarr[3:*,*]
nindexcol = (size(index,/dimensions))[0]

; Get unique star IDs and X/Y values
ui = uniq(allid)
uid = allid[ui]
ux = allx[ui]
uy = ally[ui]
nstars = n_elements(uid)
; Create structure and stuff in ID, X, Y
str = replicate({id:0L,x:0.0,y:0.0,index:lonarr(nfiles)},nstars)
str.id = uid
str.x = ux
str.y = uy

; Load the INDEX information
nlinesperstar = ntlines / nstars

; Load the index information
uindex = lonarr(nfiles,nstars)
For i=0,nstars-1 do begin
  ; Pull out the Nlinesperstar from INDEX for this star
  index1 = index[*,i*nlinesperstar:(i+1)*nlinesperstar-1]
  ; Reformat from 2D to 1D
  index2 = reform(index1,nlinesperstar*nindexcol)
  ; Take only the lines we need
  ;   there might be extra blank elements
  uindex[*,i] = long(index2[0:nfiles-1])
Endfor
str.index = uindex

;nstars = n_elements(tlines)
;str = replicate({id:0L,x:0.0,y:0.0,index:lonarr(nfiles)},nstars)
;str.id = reform(tarr[0,*])
;str.x = reform(tarr[1,*])
;str.y = reform(tarr[2,*])
;str.index = tarr[3:*,*]

if keyword_set(stp) then stop

end
