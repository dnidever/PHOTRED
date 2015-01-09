pro writeals,outfile,str,head,stp=stp

;+
;
; WRITEALS
;
; This writes a photometry structure in the ALLSTAR output
; format.
;
; INPUTS:
;  outfile  The name of the ALS output file
;  str      The input photometry structure.  This should have
;            the following data (in the same order):
;            ID, X, Y, MAG, ERR, SKY, NITER, CHI, SHARP
;  head     The ALS header
;
; OUTPUTS:
;  The ALS file to "outfile"
;
; USAGE:
;  IDL>writeals,'temp.als',str,head
;
; By D.Nidever  September 2007
;-

; Not enough inputs
if n_elements(str) eq 0 or n_elements(outfile) eq 0 or n_elements(head) eq 0 then begin
  print,'Syntax - writeals,outfile,str,head,stp=stp'
  return
endif

; Checking inputs
zparcheck,'WRITEALS',outfile,1,7,0,'OUTFILE'
zparcheck,'WRITEALS',str,2,8,1,'STR'
zparcheck,'WRITEALS',head,3,7,1,'HEAD'


; Head must have at least 2 elements
if n_elements(head) lt 2 then begin
  print,'HEAD must have at least 2 elements'
  return
end

tags = tag_names(str)
ntags = n_elements(tags)
nstr = n_elements(str)

; Getting the tag types
type = lonarr(ntags)
for i=0,ntags-1 do type[i] = size(str[0].(i),/type)

; Checking each tag type
zparcheck,'WRITEALS',str[0].(0),2,[2,3,4,5],0,'ID not correct type'
zparcheck,'WRITEALS',str[0].(1),2,[2,3,4,5],0,'X not correct type'
zparcheck,'WRITEALS',str[0].(2),2,[2,3,4,5],0,'Y not correct type'
zparcheck,'WRITEALS',str[0].(3),2,[2,3,4,5],0,'MAG not correct type'
zparcheck,'WRITEALS',str[0].(4),2,[2,3,4,5],0,'ERR not correct type'
zparcheck,'WRITEALS',str[0].(5),2,[2,3,4,5],0,'SKY not correct type'
zparcheck,'WRITEALS',str[0].(6),2,[2,3,4,5],0,'NITER not correct type'
zparcheck,'WRITEALS',str[0].(7),2,[2,3,4,5],0,'CHI not correct type'
zparcheck,'WRITEALS',str[0].(8),2,[2,3,4,5],0,'SHARP not correct type'


; Opening the file
openw,unit,/get_lun,outfile

; Print header
printf,unit,head[0]
printf,unit,head[1]
printf,unit,' '

; Setting up the command
; From allstar.f
; 321   FORMAT (I7, 2A9, F9.3, F9.4, A9, F9.0, 2F9.3)
; ID, X, Y, MAG, ERR, SKY, NITER, CHI, SHARP
;format = '(I7,2A9,F9.3,F9.4,A9,F9.0,2F9.3)'
;format = '(I7,2F9.3,F9.3,F9.4,F9.3,F9.0,2F9.3)'
format = '(I7,2F9.3,F9.3,F9.4,F10.3,F9.0,2F9.3)'
; Sky needs at least 10 characters.

for i=0.,nstr-1 do $
  printf,unit,format=format,str[i].(0),str[i].(1),str[i].(2),str[i].(3),str[i].(4),str[i].(5),str[i].(6),str[i].(7),str[i].(8)

; Closing the file
close,unit
free_lun,unit

if keyword_set(stp) then stop

end
