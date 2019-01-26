;+
;
; FILE_ISFITs
;
; Check if this file is a FITS file or not.
;
; INPUTS:
;  filename   The name of the file to check.
;
; OUTPUTS:
;  return     1 if the file is a FITS file and 0 otherwise.
;
; USAGE:
;  IDL>test = file_isfits(filename)
;
; By D. Nidever, Jan 2019
;   Based partially from is_fits.pro by Dave Bazell
;-

function file_isfits,filename

  ;; Does the file exist
  if file_test(filename) eq 0 then return,0

  ;; Four possible possibilities:
  ;; 1) Regular FITS file (this includes fpacked FITS files)
  ;; 2) Gzipped FITS file
  ;; 3) ASCII file
  ;; 4) Some other binary file

  ;; Open the file normally
  OPENR,unit,filename,/get_lun
  ;; Check if the file is empty 
  if EOF(unit) then return,0
  ;; Read the beginning of the file
  hdr = bytarr(80,/nozero)
  ON_IOERROR,badread
  READU,unit,hdr
  head = string(hdr)
  if strmid(head,0,9) eq 'SIMPLE  =' then return,1

  ;; Check if its a gzipped FITS file
  OPENR,unit,filename,/get_lun,/compress
  ;; Read the beginning of the file
  hdr = bytarr(80,/nozero)
  ON_IOERROR,badread  
  READU,unit,hdr
  head = string(hdr)
  if strmid(head,0,9) eq 'SIMPLE  =' then return,1 

  ;; Not a FITS file
  badread:
  return,0

  end
