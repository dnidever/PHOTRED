pro ps2gif,filename,rot=rot,eps=eps

;+
;
; PS2GIF
;
; This program converts a postscript image to gif
; using the CONVERT function, and also allows for
; rotations.
;
; INPUTS:
;  filename  Filename with ending
;  =rot      Angle to rotate by
;  /eps      This is an EPS file.
;
; OUTPUTS:
;  GIF version of a PS or EPS file
;
; USAGE:
;  IDL>ps2gif,'plot.ps',rot=-90
;
; By D.Nidever  2006
;-

if n_params(0) eq 0 then begin
  print,'Syntax - ps2gif,filename,rot=rot,eps=eps'
  return
endif

files = file_search(filename)
if (files(0) eq '') then begin
  print,'FILE '+filename+' DOES NOT EXIST'
  return
endif
nfiles = n_elements(files)

; Multiple files to convert
if nfiles gt 1 then begin
  for i=0,nfiles-1 do ps2gif,files(i),rot=rot,eps=eps
  return
end

;file = first_el(strsplit(filename,'/',/extract),/last)
;lo = strpos(filename,file)
;dir = strmid(filename,0,lo)

; PS or EPS
if keyword_set(eps) then begin
  hi = strpos(filename,'.eps')
  gfile = strmid(filename,0,hi)+'.gif'
endif else begin
  hi = strpos(filename,'.ps')
  gfile = strmid(filename,0,hi)+'.gif'
endelse

if (hi(0) eq -1) then begin
  print,'FILE '+filename+' IS NOT A PS or EPS FILE'
  return
end

;f = strmid(file,0,hi)

; Converting to GIF
;spawn,'convert '+filename+' '+gfile
spawn,['convert','-flatten',filename,gfile],out,outerr,/noshell

; Rotating
;if (n_elements(rot) ne 0) then spawn,'mogrify -rotate '+strtrim(rot,2)+' '+gfile
if (n_elements(rot) ne 0) then spawn,['mogrify','-rotate',strtrim(rot,2),gfile],out,errout,/noshell

end
