;+
;
; PDFCOMBINE
;
; Combines multiple single-page PDF files into a multi-page PDF.
;
; INPUTS:
;  input    An array of filenames, a wild card expression (e.g., *.pdf)
;             or the name of a file with a list of names (e.g., @list.txt).
;             The normal mode is to use PDF file but if eps file names
;             are input and the respective PDFs don't exist then they
;             will be created.
;  outfile  The output file name.
;  /reconv  Redo the conversion to PDF if EPS filenames input.
;  /clobber Overwrite the output file if it already exists.
;
; OUTPUTS:
;  The combined multi-page PDF is written to OUTFILE.
;
; USAGE:
;  IDL>pdfcombine,['file1.pdf','file2.pdf'],'outfile.pdf'
;
; By D. Nidever,  August 2018
;-

pro pdfcombine,input,outfile,clobber=clobber,reconv=reconv

;; Not enough inputs
if n_elements(input) eq 0 or n_elements(outfile) eq 0 then begin
  print,'Syntax - pdfcombine,input,outfile,clobber=clobber,reconv=reconv'
  return
endif

;; Output file exists already and /clobber not set
if file_test(outfile) eq 1 and not keyword_set(clobber) then begin
  print,outfile,' already exist and /clobber NOT set'
  return
endif

;; Load the input
LOADINPUT,input,files,count=nfiles

;; Create PDFs if necessary
pdffiles = files
for i=0,nfiles-1 do begin
  file = files[i]
  ext = first_el(strsplit(file_basename(file),'.',/extract),/last)
  if ext eq 'eps' then begin
    pdffile = repstr(file,'.eps','.pdf')
    if file_test(pdffile) eq 0 or keyword_set(reconv) then begin
      print,'Converting ',file,' to PDF'
      spawn,['epstopdf',file],/noshell 
      pdffiles[i] = pdffile
    endif
  endif
endfor

;; Combine the PDFs
cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+outfile+' '
cmd += strjoin(pdffiles,' ')
spawn,cmd,out,errout
if errout[0] ne '' then begin
  print,'There was an error combining the PDFs:  ',errout
endif

;; Check that the output exists
if file_test(outfile) eq 0 then begin
  print,outfile,' NOT FOUND'
endif else begin
  print,'Combined PDF written to ',outfile
endelse

;stop

end
