;-----------------------------------------------------------------------
;+
; NAME:
;   LOCKFILE
;
; PURPOSE:
;   Test if a file is already "locked", and lock it if not.
;
; CALLING SEQUENCE:
;   res = lockfile( filename )
;
; INPUT:
;   filename:   File name
;
; OUTPUTS:
;   res:        Return 0 if file already locked, or 1 if not in which case
;               we would have just locked it.
;
; COMMENTS:
;   A lock file is created with a single byte written to it.
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel, APO
;   26-Nov-2017  Modified by D. Nidever to take the lock filename as input
;                 and only create the file no link.
;-
;-----------------------------------------------------------------------
function lockfile, filename

   if (n_elements(filename) NE 1 OR size(filename,/tname) NE 'STRING') $
    then begin
      print, 'FILENAME must be specified'
      return, 0
   endif

   openw, olun, filename, /append, /get_lun, error=err
   if ((fstat(olun)).size EQ 0) or (err ne 0) then begin
      writeu, olun, 1B
      flush, olun ; Flush output immediately
      res = 1
   endif else begin
      res = 0
   endelse
   close, olun ; This will flush output to this file
   free_lun, olun

   return, res
end
;-----------------------------------------------------------------------
