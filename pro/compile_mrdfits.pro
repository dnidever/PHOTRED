pro compile_mrdfits

;; Make sure MRDFITS is compiled

;; Create a small temporary file and load it
tmpfile = MKTEMP('tmp',/nodot) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
tmpfile += '.fits'
str = replicate({a:'',b:0.0},2)
MWRFITS,str,tmpfile,/create,/silent
dum = MRDFITS(tmpfile,1,/silent)

end
