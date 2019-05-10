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
;  /header   Only return the header/metadata.
;  /nowrite  Don't write the file if using the resource file.
;
; OUTPUTS:
;  results   The primary data in the file
;  meta     Meta-data in the 1st extension if it exists.
;  =count    The number of elements.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>cat = photred_readfile('F4-20440011_09.als',meta)
;
; By D. Nidever, Jan 2019
;-

function photred_readfile,filename,meta,exten=exten,count=count,header=header,nowrite=nowrite,error=error

  undefine,error
  undefine,meta
  count = -1
  
  ;; Not enough inputs
  if n_elements(filename) eq 0 then begin
    error = 'Not enough inputs'
    print,"Syntax - cat = photred_readfile('F4-20440011_09.als',meta)"
    return,-1
  endif

  ;; Check that file exists
  if file_test(filename[0]) eq 0 then begin
    error = filename+' NOT FOUND'
    print,error
    return,-1
  endif

  ;; Header-only
  if keyword_set(header) then count=0

  ;; Is this a FITS file
  isfits = file_isfits(filename[0])

  ;; Break down the file into its components
  info = file_info(filename)
  dir = file_dirname(filename[0])
  base = file_basename(filename[0])
  ext = first_el(strsplit(base,'.',/extract),/last)
  if isfits eq 1 and ext eq 'gz' then ext='fits'  ; gzipped FITS file
  if isfits eq 1 and ext eq 'fz' then ext='fits'  ; fpacked FITS file  
  if isfits eq 1 then ext='fits'      ; fits file with different extension
  
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
         str = replicate({file:'',trans:dblarr(ntrans),magoff:fltarr(2)},count)
         str.file = files
         str.trans = transpose(trans)
         str.magoff = magoff
         return,str
       end
    'raw': begin
         if keyword_set(header) then begin   ;; only read header
           READLINE,filename,meta,nreadline=2
           return,meta
         endif
         LOADRAW,filename,phot,meta
         count = n_elements(phot)
         return,phot
       end
    'tfr': begin
         LOADTFR,filename,files,str
         count = n_elements(files)
         meta = files
         return,str
      end
    'als': begin
         if keyword_set(header) then begin   ;; only read header
           READLINE,filename,meta,nreadline=2
           return,meta
         endif
         LOADALS,filename,phot,meta,count=count
         return,phot
       end
    'alf': begin
         if keyword_set(header) then begin   ;; only read header
           READLINE,filename,meta,nreadline=2
           return,meta
         endif
         LOADALS,filename,phot,meta,count=count
         return,phot
       end
    'mag': begin
         if isfits eq 1 then begin
           if keyword_set(header) then begin   ;; only read header
             meta = HEADFITS(filename)
             return,meta
           endif
           phot = MRDFITS(filename,1,/silent)
           count = n_elements(phot)
           fits_open,filename,fcb & fits_close,fcb
           if fcb.nextend ge 2 then meta=MRDFITS(filename,2,/silent)
           return,phot
         endif else begin
           if keyword_set(header) then begin   ;; only read header
             READLINE,filename,meta,nreadline=1
             return,meta
           endif
           LOADMAG,filename,phot
           count = n_elements(phot)
           return,phot
         endelse
       end
    'ast': begin
         if isfits eq 1 then begin
           fits_open,filename,fcb & fits_close,fcb
           if fcb.nextend ge 2 then meta=MRDFITS(filename,2,/silent)
           if keyword_set(header) then begin   ;; only read header
             if n_elements(meta) eq 0 then meta = HEADFITS(filename)
             return,meta
           endif
           phot = MRDFITS(filename,1,/silent)
           count = n_elements(phot)
           return,phot
         endif else begin
           if keyword_set(header) then begin   ;; only read header
             READLINE,filename,meta,nreadline=1
             return,meta
           endif
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
           fits_open,filename,fcb & fits_close,fcb
           if fcb.nextend ge 2 then meta=MRDFITS(filename,2,/silent)
           if keyword_set(header) then begin   ;; only read header
             if n_elements(meta) eq 0 then meta = HEADFITS(filename)
             return,meta
           endif
           phot = MRDFITS(filename,1,/silent)
           count = n_elements(phot)
           return,phot
         endif else begin
           if keyword_set(header) then begin   ;; only read header
             READLINE,filename,meta,nreadline=1
             return,meta
           endif
           ; Get the fieldnames and fieldtypes
           ; We need ID to be STRING not LONG
           tempfile = MKTEMP('temp',outdir='/tmp')
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
           return,phot
         endelse
       end
    'cmb': begin
         if isfits eq 1 then begin
           fits_open,filename,fcb & fits_close,fcb
           if fcb.nextend ge 2 then meta=MRDFITS(filename,2,/silent)
           if keyword_set(header) then begin   ;; only read header
             if n_elements(meta) eq 0 then meta = HEADFITS(filename)
             return,meta
           endif
           phot= MRDFITS(filename,1,/silent)
           count = n_elements(phot)
           return,phot
         endif else begin
           if file_test(filename+'.meta') eq 1 then meta=IMPORTASCII(filename+'.meta',/header,/silent)
           if keyword_set(header) then begin   ;; only read header
             if n_elements(meta) eq 0 then READLINE,filename,meta,nreadline=1
             return,meta
           endif
           phot = IMPORTASCII(filename,/header,/silent)
           count = n_elements(phot)           
           return,phot          
         endelse
       end
    'dered': begin
         if isfits eq 1 then begin
           fits_open,filename,fcb & fits_close,fcb
           if fcb.nextend ge 2 then meta=MRDFITS(filename,2,/silent)
           if keyword_set(header) then begin   ;; only read header
             if n_elements(meta) eq 0 then meta = HEADFITS(filename)
             return,meta
           endif
           phot = MRDFITS(filename,1,/silent)
           count = n_elements(phot)
           return,phot
        endif else begin
           if file_test(filename+'.meta') eq 1 then meta=IMPORTASCII(filename+'.meta',/header,/silent)
           if keyword_set(header) then begin   ;; only read header
             if n_elements(meta) eq 0 then READLINE,filename,meta,nreadline=1
             return,meta
           endif
           phot = IMPORTASCII(filename,/header,/silent)
           count = n_elements(phot)
           return,phot           
         endelse
       end
    'final': begin
         if isfits eq 1 then begin
           fits_open,filename,fcb & fits_close,fcb
           if fcb.nextend ge 2 then meta=MRDFITS(filename,2,/silent)
           if keyword_set(header) then begin   ;; only read header
             if n_elements(meta) eq 0 then meta = HEADFITS(filename)
             return,meta
           endif
           phot = MRDFITS(filename,1,/silent)
           count = n_elements(phot)
           return,phot
         endif else begin
           if file_test(filename+'.meta') eq 1 then meta=IMPORTASCII(filename+'.meta',/header,/silent)
           if keyword_set(header) then begin   ;; only read header
             if n_elements(meta) eq 0 then READLINE,filename,meta,nreadline=1
             return,meta
           endif
           phot = IMPORTASCII(filename,/header,/silent)
           count = n_elements(phot)
           return,phot
         endelse
       end
    'fits': begin
         ;; Load using resource file
         rfilename = dir+'/.'+base
         rinfo = file_info(rfilename)
         ;;  only use if FITS file does not exist
         if info.size le 1 and rinfo.exists eq 1 then begin
           result = FITS_READ_RESOURCE(filename,meta,header=header,nowrite=nowrite)
           count = 1
           return,result
         endif
         ;; Regular FITS file
         if n_elements(exten) gt 0 then begin
            if keyword_set(header) then begin
              meta = HEADFITS(filename,exten=exten)
              return,meta
            endif
            result = MRDFITS(filename,exten,meta,/silent)
            count = n_elements(result)
            return,result
         endif else begin
            ;; Figure out whether to read HDU0 or HDU1
            fits_open,filename,fcb & fits_close,fcb
            hd0 = HEADFITS(filename,exten=0,/silent)
            if fcb.nextend ge 1 then hd1=HEADFITS(filename,exten=1,/silent)
            ;; Only HDU0 exists
            if fcb.nextend eq 0 then begin
               if keyword_set(header) then begin
                 meta = HEADFITS(filename)
                 return,meta
               endif
               result = MRDFITS(filename,0,meta,/silent)
               count = n_elements(result)
               return,result
            endif
            ;; Both exist, HDU0 has no data but HDU1 does
            if sxpar(hd0,'NAXIS') eq 0 and sxpar(hd1,'NAXIS') gt 0 then begin
               if keyword_set(header) then begin
                 meta = HEADFITS(filename,exten=1)
                 return,meta
               endif
               result = MRDFITS(filename,1,meta,/silent)
               count = n_elements(result)
               return,result
            endif    
            ;; Both exist, both HDU0 and HDU1 have data (meta-data in HDU1)
            if sxpar(hd0,'NAXIS') gt 0 and sxpar(hd1,'NAXIS') gt 0 then begin
               meta = MRDFITS(filename,1,/silent)
               if keyword_set(header) then return,meta
               result = MRDFITS(filename,0,/silent)
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
