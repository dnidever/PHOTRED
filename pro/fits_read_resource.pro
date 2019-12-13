;+
;
; FITS_READ_RESOURCE
;
; Read a FITS file by using the resource file.
;
; INPUTS:
;  file      The FITS filename to read.  This is the actual file, not
;               the resource file.
;  /header   Return the header only.  The header will be returned
;               as the output of the function.
;  /nowrite  Only return the data and don't write the file to the disk.
;
; OUTPUTS:
;  im        The 2D FITS image.
;  meta     The header array for IM.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>im = fits_read_resource(file,meta=meta,header=header,nowrite=nowrite)
;
; By D. Nidever  Feb 2019
;-

function fits_read_resource,file,meta,header=header,nowrite=nowrite,error=error

undefine,im,meta,error

;; Not enough inputs
if n_elements(file) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - im=fits_read_resource(file,meta,header=header,nowrite=nowrite)'
  return,-1
endif

;; It takes about 1 sec

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   error = !ERROR_STATE.MSG
   PHOTRED_ERRORMSG,logfile=logf
   CATCH, /CANCEL 
   return,-1
endif

t0 = systime(1)

dir = file_dirname(file)
base = file_basename(file)
rfile = dir+'/.'+base

;; Check if there's a lock file
lockinfo = file_info(file+'.lock')
if lockinfo.exists eq 1 and t0-lockinfo.ctime lt 100. then begin
  print,'Lock file.  Waiting'
  while (file_test(file+'.lock') eq 1 and systime(1)-t0 lt 100.) do wait,5
endif

;; If the file exists and is not empty then just return it
;;=======================================================
info = file_info(file)
if info.exists eq 1 and info.size gt 1 then begin
  if keyword_set(header) then begin
    meta = HEADFITS(file)
    return,meta
  endif
  FITS_READ,file,im,meta
  return,im
endif

;; If the file is empty then use the resource information
;;=======================================================

;; Create lock file
if not keyword_set(header) then TOUCHZERO,file+'.lock'

;; Load the resource file
READLINE,rfile,rlines,count=nlines,/noblank,comment='#'
arr = strtrim(strsplitter(rlines,'=',/extract),2)
names = reform(arr[0,*])
vals = reform(arr[1,*])
rstr = create_struct(names[0],vals[0])
if nlines gt 1 then for k=1,nlines-1 do rstr=create_struct(rstr,names[k],vals[k])


;; Only the header, local
;; get it from the resource file or a stand-alone file
if keyword_set(header) and tag_exist(rstr,'HEADER') then begin
  if strmid(rstr.header,0,1) eq '/' then hfile=rstr.header else hfile = dir+'/'+rstr.header
  if file_test(hfile) eq 0 then begin
    error = 'Local header file '+hfile+' NOT FOUND'
    print,error
  endif else begin
    READLINE,hfile,meta
    return,meta
  endelse
endif


;; Create a temporary directory to work in
if not keyword_set(header) then begin
  tmpdir = MKTEMP('rsrc',/directory,/nodot)
  FILE_CHMOD,tmpdir,/a_execute
endif

;; Get the flux file
;;   /mss1/archive/pipe/20170328/ct4m/2014B-0404/c4d_170329_043921_ooi_g_v1.fits.fz[39]
lo = strpos(rstr.fluxfile,'[')
hi = strpos(rstr.fluxfile,']')
fluxfile = strmid(rstr.fluxfile,0,lo)
fext = strmid(rstr.fluxfile,lo+1,hi-lo-1)
if not keyword_set(header) then begin
  tfluxfile = tmpdir+'/flux.fits'
  SPAWN,['funpack','-E',fext,'-O',tfluxfile,fluxfile],/noshell
endif

;; Construct the header from the extension and the main headers:
;;                       <required keywords for the extension:
;;                       XTENSION, BITPIX,
;;                               NAXIS, ...>
;;                       BEGIN MAIN HEADER
;;                       --------------------------------
;;                       <PDU header keyword and history less required
;;                       keywords:
;;                               SIMPLE, BITPIX, NAXIS, ...>
;;                       BEGIN EXTENSION HEADER
;;                       ---------------------------
;;                       <extension header less required keywords that
;;                       were
;;                               placed at the beginning of the header.
;;                       END
;; Need PDU header with exposure information
mhead0 = HEADFITS(fluxfile,exten=0,errmsg=errmsg0)
if keyword_set(header) then ehead0 = HEADFITS(fluxfile,exten=fext,errmsg=errmsg1) else $
  ehead0 = HEADFITS(tfluxfile,exten=0,errmsg=errmsg1)
;; Required keywords
;XTENSION= 'IMAGE   '           /extension type                                  
;   SIMPE = T
;BITPIX  =                  -32 /bits per data value                             
;NAXIS   =                    2 /number of axes                                  
;NAXIS1  =                 2046 /                                                
;NAXIS2  =                 4094 /                                                
;PCOUNT  =                    0 /Number of group parameters                      
;GCOUNT  =                    1 /Number of groups  
;; Start the final header
undefine,head
head = ['SIMPLE  =                    T / file does conform to FITS standard             ']
sxaddpar,head,'BITPIX',sxpar(ehead0,'BITPIX'),'bits per data value'
sxaddpar,head,'NAXIS',sxpar(ehead0,'NAXIS'),'number of data axes'
sxaddpar,head,'NAXIS1',sxpar(ehead0,'NAXIS1'),'length of data axis 1'
sxaddpar,head,'NAXIS2',sxpar(ehead0,'NAXIS2'),'length of data axis 2'
sxaddpar,head,'PCOUNT',0,'Number of group parameters'
sxaddpar,head,'GCOUNT',1,'Number of groups'
sxdelpar,head,'END'  ;; sxaddpar adds an "END" line automatically
;; Remove required keywords from the main header
mhead = mhead0
todel = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','EXTEND','DATASUM','CHECKSUM','END']
for j=0,n_elements(todel)-1 do sxdelpar,mhead,todel[j]
;; Remove required keywords from the extension header
ehead = ehead0
todel = ['SIMPLE','XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','EXTEND','DATASUM','CHECKSUM','END']
for j=0,n_elements(todel)-1 do sxdelpar,ehead,todel[j]
;; Combine them all
PUSH,head,'COMMENT BEGIN MAIN HEADER ---------------------------------                     '
PUSH,head,mhead
PUSH,head,'COMMENT BEGIN EXTENSION HEADER ----------------------------                     '
PUSH,head,ehead
PUSH,head,'END                                                                             '
;; Remove any blank lines
bd = where(strtrim(head,2) eq '',nbd)
if nbd gt 0 then REMOVE,bd,head

;; Fix/remove FPACK keywords
if keyword_set(header) then begin
  sxaddpar,head,'BITPIX',sxpar(head,'ZBITPIX')
  sxdelpar,head,'ZBITPIX'
  sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
  sxdelpar,head,'ZNAXIS'
  sxdelpar,head,'ZNAXIS1'
  sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
  todel = ['ZNAXIS','ZNAXIS1','ZNAXIS2','ZIMAGE','ZTILE1','ZTILE2','ZCMPTYPE','ZNAME1','ZVAL1',$
           'ZNAME2','ZVAL2','ZPCOUNt','ZGCOUNT','ZQUANTIZ','ZDITHER0','ZTENSION','TFIELDS','TTYPE1','TFORM1',$
           'TTYPE2','TFORM2','TTYPE3','TFORM3']
  for j=0,n_elements(todel)-1 do sxdelpar,head,todel[j]

  ;; DAOPHOT_IMPREP adds GAIN and RDNOISE
  if keyword_set(header) then begin
    gainA = sxpar(head,'GAINA')
    gainB = sxpar(head,'GAINB')
    gain = (gainA+gainB)*0.5
    rdnoiseA = sxpar(head,'RDNOISEA')
    rdnoiseB = sxpar(head,'RDNOISEB')
    rdnoise = (rdnoiseA+rdnoiseB)*0.5
    sxaddpar,head,'GAIN',gain
    sxaddpar,head,'RDNOISE',rdnoise
  endif
endif

;; This is done in daophot_imprep.pro but if we are only returning
;; the header then do it here
if keyword_set(header) then begin
  saturate = sxpar(head,'saturate',count=nsaturate)
  if nsaturate gt 0 then saturate<=64500.0 else saturate=64500.0  ; set it slightly lower than 65000 for DAOPHOT
  sxaddpar,head,'saturate',saturate
endif

;; Returning only the header
if keyword_set(header) then return,head


;; Now update the fluxfile with this new header
MODFITS,tfluxfile,0,head

;; Get the mask file
lo = strpos(rstr.maskfile,'[')
hi = strpos(rstr.maskfile,']')
maskfile = strmid(rstr.maskfile,0,lo)
mext = strmid(rstr.maskfile,lo+1,hi-lo-1)
tmaskfile = tmpdir+'/mask.fits'
SPAWN,['funpack','-E',mext,'-O',tmaskfile,maskfile],/noshell

;; Create the DAOPHOT file
;;   taken from smashred_imprep_single.pro
DAOPHOT_IMPREP,tfluxfile,tmaskfile,im,meta,header=header,error=error
if n_elements(error) gt 0 then return,-1


;; Use the local header
;;   to write to the file and return
if tag_exist(rstr,'HEADER') then begin
  if strmid(rstr.header,0,1) eq '/' then hfile=rstr.header else hfile = dir+'/'+rstr.header
  if file_test(hfile) eq 0 then begin
    error = 'Local header file '+hfile+' NOT FOUND'
    print,error
  endif else begin
    READLINE,hfile,meta
    meta0 = meta
    READLINE,rstr.header,meta
  endelse
endif

;;; Put in FPACK parameters
;if keyword_set(fpack) then begin
;  ; Remove all past FZ parameters
;  bd = where(strmid(newhead,0,2) eq 'FZ',nbd)
;  if nbd gt 0 then REMOVE,bd,newhead
;  sxaddpar,newhead,'FZALGOR','RICE_1'
;  sxaddpar,newhead,'FZQMETHD','SUBTRACTIVE_DITHER_1'
;  sxaddpar,newhead,'FZQVALUE',8
;  sxaddpar,newhead,'FZDTHRSD','CHECKSUM'
;endif

;; Write new image
if not keyword_set(nowrite) and not keyword_set(header) then $
  MWRFITS,im,file,meta,/create

;; Delete the lock file
FILE_DELETE,file+'.lock',/allow

;; Fpack compress
;if keyword_set(fpack) then begin
;  SPAWN,['fpack',newfile],/noshell
;  if file_test(newfile+'.fz') eq 0 then begin
;    print,'ERROR: Failure creating '+newfile+'.fz'
;    retall
;  endif
;  ; Delete normal FITS file
;  FILE_DELETE,newfile
;endif


;; Clean up
FILE_DELETE,tfluxfile,tmaskfile,/allow
FILE_DELETE,tmpdir

dt = systime(1)-t0
;print,dt

;; Return header only
if keyword_set(header) then return,meta

return,im

end
