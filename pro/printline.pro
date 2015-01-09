pro printline,lines,number=number,boundary=boundary,stp=stp

;+
; This prints a string array to the screen. Similar to WRITELINE.PRO
;
; INPUTS:
;  lines     String array
;  /number   Print the lines number at the beginning of the line
;  /boundary Adds boundary characters "|" at the beg/end to show you
;            where exactly the string begins/ends.
;  /stp      Stop at the end of the program
;
; OUTPUTS:
;  The output written to the screen
;
; By D.Nidever  Feb.2007
;-

if n_elements(lines) eq 0 then begin
  print,'Syntax - printline,lines,stp=stp'
  return
endif

nlines = n_elements(lines)

if keyword_set(boundary) then bnd='|' else bnd=''

; Regular
if not keyword_set(number) then begin
  for i=0.,nlines-1 do print,bnd+lines[i]+bnd

; With line numbers
endif else begin
  for i=0.,nlines-1 do print,strtrim(long(i),2)+' - ',bnd+lines[i]+bnd
endelse

if keyword_set(stp) then stop

end
