;+
;
; PHOTRED_READFILE
;
; Generic file reading program for PHOTRED
;
; INPUTS:
;  filename  The name of the file to read from.  The extension
;              and type of the file is used to figure out how
;              to read it.
;  =exten    The extension to read
;
; OUTPUTS:
;  results   The primary data in the file
;  =meta     Meta-data in the 1st extension if it exists.
;  =count    The number of elements.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>cat = photred_readfile('F4-20440011_09.als')
;
; By D. Nidever, Jan 2019
;-

function photred_readfile,filename,meta=meta,count=count,error=error

  undefine,error
  undefine,meta
  count = -1
  
  ;; Check that file exists
  if file_test(filename[0]) eq 0 then begin
    error = filename+' NOT FOUND'
    print,error
    return,-1
  endif

  ;; Is this a FITS file
  isfits = file_isfits(filename[0])

  ;; Break down the file into its components
  dir = file_dirname(filename[0])
  base = file_basename(filename[0])
  ext = first_el(strsplit(base,'.',/extract),/last)
  if isfits eq 1 and ext eq 'gz' then ext='fits'  ; gzipped FITS file
  if isfits eq 1 and ext eq 'fz' then ext='fits'  ; fpacked FITS file  
  
  ;; Go through the various options
  CASE ext of
     'opt': begin
         READCOL,filename,opttags,dum,optvals,format='A,A,F',/silent
         count = n_elements(opttags)
         str = replicate({name:'',value:0.0},count)
         str.name = opttags
         str.value = optvals
         return,str
       end
    'coo': begin
         LOADCOO,filename,str,meta
         count = n_elements(str)
         return,str
       end
    'ap': begin
         LOADAPER,filename,phot,meta
         count = n_elements(phot)
         return,phot
       end
    'mch': begin
         LOADMCH,filename,files,trans,magoff,count=count
         ntrans = n_elements(trans[0,*])
         str = replicate({file:'',trans:fltarr(ntrans),magoff:fltarr(2)},count)
         str.file = files
         str.trans = transpose(trans)
         str.magoff = magoff
         return,str
       end
    'raw': begin
         LOADRAW,filename,phot,meta
         count = n_elements(phot)
         return,raw
       end
    'tfr': begin
         LOADTFR,filename,files,str
         count = n_elements(files)
         meta = files
         return,str
      end
    'als': begin
         LOADALS,filename,phot,meta,count=count
         return,phot
       end
    'alf': begin
         LOADALS,filename,phot,meta,count=count
         return,phot
       end
    'mag': begin
         if isfits eq 1 then begin
           phot = MRDFITS(filename,1,/silent)
           dum = HEADFITS(filename,exten=2,errmsg=errmsg,/silent) 
           if errmsg ne '' then meta=MRDFITS(filename,2,/silent)
           return,phot
         endif else begin
           LOADMAG,filename,phot
           count = n_elements(phot)
           return,phot
         endelse
       end
    'ast': begin
         if isfits eq 1 then begin
           phot = MRDFITS(filename,1,/silent)
           dum = HEADFITS(filename,exten=2,errmsg=errmsg,/silent) 
           if errmsg ne '' then meta=MRDFITS(filename,2,/silent)
           return,phot
         endif else begin
           phot = IMPORTASCII(filename,/header,/silent)
           count = n_elements(phot)
           return,phot
         endelse
       end
    'trans': begin
         READ_TRANS,filename,trans,error=error
         count = n_elements(trans)
         return,trans
       end
    'phot': begin
         if isfits eq 1 then begin
           phot = MRDFITS(filename,1,/silent)
           dum = HEADFITS(filename,exten=2,errmsg=errmsg,/silent) 
           if errmsg ne '' then meta=MRDFITS(filename,2,/silent)
           return,phot
         endif else begin
           ; Get the fieldnames and fieldtypes
           ; We need ID to be STRING not LONG
           tempfile = MKTEMP('temp')
           READLINE,filename,lines,nlineread=40 
           WRITELINE,tempfile,lines
           temp = IMPORTASCII(tempfile,/header,/silent)
           FILE_DELETE,tempfile,/allow,/quiet
           fieldnames = TAG_NAMES(temp)
           nfieldnames = n_elements(fieldnames)
           fieldtypes = lonarr(nfieldnames)
           for k=0,nfieldnames-1 do fieldtypes[k] = SIZE(temp[0].(k),/type)
           idind = where(fieldnames eq 'ID',nidind)
           fieldtypes[idind[0]] = 7
           phot = IMPORTASCII(filename,fieldnames=fieldnames,fieldtypes=fieldtypes,skip=1,/noprint)
           count = n_elements(phot) 
         endelse
       end
    'cmb': begin
         if isfits eq 1 then begin
           phot = MRDFITS(filename,1,/silent)
           dum = HEADFITS(filename,exten=2,errmsg=errmsg,/silent) 
           if errmsg ne '' then meta=MRDFITS(filename,2,/silent)
           return,phot
         endif else begin
           phot = IMPORTASCII(filename,/header,/silent)
           count = n_elements(phot)           
           return,phot          
         endelse
       end
    'dered': begin
         if isfits eq 1 then begin
           phot = MRDFITS(filename,1,/silent)
           dum = HEADFITS(filename,exten=2,errmsg=errmsg,/silent) 
           if errmsg ne '' then meta=MRDFITS(filename,2,/silent)
           return,phot
        endif else begin
           phot = IMPORTASCII(filename,/header,/silent)
           count = n_elements(phot)
           return,phot           
         endelse
       end
    'final': begin
         if isfits eq 1 then begin
           phot = MRDFITS(filename,1,/silent)
           dum = HEADFITS(filename,exten=2,errmsg=errmsg,/silent) 
           if errmsg ne '' then meta=MRDFITS(filename,2,/silent)
           return,phot
         endif else begin
           phot = IMPORTASCII(filename,/header,/silent)
           count = n_elements(phot)
           return,phot
         endelse
       end
    'fits': begin
         if n_elements(exten) gt 0 then begin
            result = MRDFITS(filename,exten,meta,/silent)
            count = n_elements(result)
            return,result
         endif else begin
            ;; Figure out whether to read HDU0 or HDU1
            hd0 = HEADFITS(filename,exten=0,errmsg=errmsg0,/silent)
            hd1 = HEADFITS(filename,exten=1,errmsg=errmsg1,/silent)            
            ;; Only HDU0 exists
            if errmsg1 ne '' then begin
               result = MRDFITS(filename,0,meta,/silent)
               count = n_elements(result)
               return,result
            endif
            ;; Both exist, HDU0 has not data but HDU1 does
            if sxpar(hd0,'NAXIS') eq 0 and sxpar(hd1,'NAXIS') gt 0 then begin
               result = MRDFITS(filename,1,meta,/silent)
               count = n_elements(result)
               return,result
            endif    
            ;; Both exist, both HDU0 and HDU1 have data (meta-data in HDU1)
            if sxpar(hd0,'NAXIS') gt 0 and sxpar(hd1,'NAXIS') gt 0 then begin
               result = MRDFITS(filename,0,/silent)
               meta = MRDFITS(filename,1,/silent)
               count = n_elements(result)
               return,result
            endif  
         endelse
       end
     else:  begin
       error = 'Extension '+ext+' not supported'
       print,error
       return,-1
    end
  ENDCASE

  return,-1

  end
