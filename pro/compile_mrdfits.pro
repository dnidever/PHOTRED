pro compile_mrdfits

;; Make sure MRDFITS is compiled
;; NOTE, you should use RESOLVE_ROUTINE,'mrdfits',/compile_full_file,/either  INSTEAD

;; Create a small temporary file and load it
tmpfile = MKTEMP('tmp',/nodot) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
tmpfile += '.fits'
str = replicate({a:'',b:0.0},2)
MWRFITS,str,tmpfile,/create,/silent
dum = MRDFITS(tmpfile,1,/silent)
file_delete,tmpfile

end
