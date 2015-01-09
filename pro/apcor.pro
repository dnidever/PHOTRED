pro apcor,input,outfile,stp=stp

;+
;
; APCOR
;
; This takes a list of .del files and makes a
; final file of aperture corrections
; This is an IDL version of apcor.cl
;
; INPUTS:
;  input    A list of the .del files.
;  outfile  The file to output the final aperture
;           correction to.  By default outfile='apcor.lst'
;
; OUTPUTS:
;  A final of final aperture corrections.
;
; USAGE:
;  IDL>apcor,'*.del'
;
; By D.Nidever   December 2006 (copy of apcor.cl)
;-

; No input
ninput = n_elements(input)
if ninput eq 0 then begin
  print,'Syntax - apcor,input,outfile'
  return
endif

noutfile = n_elements(outfile)
if noutfile eq 0 then outfile='apcor.lst'

; Loading the input
LOADINPUT,input,names
nnames = n_elements(names)

; How long are the names
len = strlen(names)
maxlen = max(len)

openw,unit,/get_lun,outfile

for i=0,nnames-1 do begin
  READCOL,names[i],id,del,format='I,F',/silent

  med = median(del)
  RESISTANT_MEAN,del,2.0,mean,sigma_mean,num_rejected

  printf,unit,format='(A-'+strtrim(maxlen+3,2)+',F12.8)',names[i],mean
end

close,unit
free_lun,unit

if keyword_set(stp) then stop

end
