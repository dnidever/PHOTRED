;+
;
; DAOMATCH
;
; This matches stars and finds the transformations using
; MATCHSTARS.PRO (originally DAOMATCH was used) and then
; combining them with DAOMASTER.  INputs need to be ALS files.
;
; INPUTS:
;  files     Array of ALS files,  The first file will be used
;            as the reference.  If /fake set then the first ALS file
;            in "files" should already have an associated MCH file.
;  =maxshift Constraints on the initial X/Y shifts.
;  /usewcs   Use the WCS for initial matching.  This is the default.
;  /fake     Run for artificial stars.  The MCH file should be input
;              and daomaster run to create raw/tfr files.
;  /verbose  Verbose output
;  /stp      Stop at the end of the program
;  =hi       Not used anymore.
;
; OUTPUTS:
;  An MCH, TFR and RAW file with basename of the first file.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>daomatch,['obj1034_1.als','obj1035_1.als','obj1036_1.als']
;
; Add options to freeze the scale and/or rotation.
; 
; By D. Nidever   December 2006
;-

pro daomatch_dummy
FORWARD_FUNCTION test_trans, trans_coo
end

;---------------------------------------------------------------

function test_trans,trans

; This function tests if a transformation equation from
; daomatch is good or not.  The scale should be nearly 1.0
; and the rotation should be near 0.
;
; Return value:
;  1  Good
;  0  Bad
; -1  Fail
;

test = -1

if n_elements(trans) eq 0 then return,-1

; The test mainly looks at the rotation/scale values
; and not the xoff/yoff values.

sz = size(trans)
if sz[0] ne 2 or sz[1] ne 2 or sz[2] ne 6 then return,-1


xoff = (reform(trans[1,0]))(0)
yoff = (reform(trans[1,1]))(0)
c = (reform(trans[1,2]))(0)
e = (reform(trans[1,3]))(0)
d = (reform(trans[1,4]))(0) 
f = (reform(trans[1,5]))(0)
; from ccdpck.txt
;              x(1) = A + C*x(n) + E*y(n)
;              y(1) = B + D*x(n) + F*y(n)

; C=F~1 and D=E~0
test = 1
if abs(c-f) gt 0.1 then test=0
if abs(d-e) gt 0.1 then test=0
if abs(c-1.0) gt 0.1 then test=0
if abs(e) gt 0.1 then test=0

return,test

end

;---------------------------------------------------------------

pro daomatch,files,usewcs=usewcs,stp=stp,verbose=verbose,hi=hi,logfile=logfile,error=error,$
             maxshift=maxshift,fake=fake

t0 = systime(1)

undefine,error

nfiles = n_elements(files)
if nfiles eq 0 then begin
  print,'Syntax - daomatch,files,usewcs=usewcs,fake=fake,stp=stp,verbose=verbose'
  error = 'Not enough inputs'
  return
end

; Default parameters
if n_elements(usewcs) eq 0 then usewcs=1  ; use WCS by default

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1

; Only one file, can't match
if nfiles eq 1 then begin
  printlog,logf,'ONLY ONE FILE INPUT.  No matching, simply creating .mch and .raw file'
  ;printlog,logf,'ONLY ONE FILE INPUT.  NEED AT *LEAST* TWO'
  ;error = 'ONLY ONE FILE INPUT.  NEED AT *LEAST* TWO'
  ;return
endif

; Compile MATCHSTARS.PRO
RESOLVE_ROUTINE,'matchstars',/compile_full_file

; Current directory
CD,current=curdir

dir = FILE_DIRNAME(files[0])
CD,dir

files2 = FILE_BASENAME(files,'.als')

; FAKE, running for artificial star tests
if keyword_set(fake) then begin

  ; Check that MCH file exists
  if file_test(files2[0]+'.mch') eq 0 then begin
    error = '/fake set but '+files2[0]+'.mch NOT FOUND'
    printlog,logf,error
    return
  endif
   
  ; Keep backup of original mch file
  FILE_DELETE,files2[0]+'.mch.orig',/allow_nonexistent
  FILE_COPY,files2[0]+'.mch',files2[0]+'.mch.orig' 
   
  ; Remove the output files
  FILE_DELETE,files2[0]+'.raw',/allow_nonexistent
  FILE_DELETE,files2[0]+'.tfr',/allow_nonexistent
   
  goto,rundaomaster
endif

; Remove the output files
FILE_DELETE,files2[0]+'.mch',/allow_nonexistent
FILE_DELETE,files2[0]+'.raw',/allow_nonexistent
FILE_DELETE,files2[0]+'.tfr',/allow_nonexistent

undefine,mchfinal

; Check that the reference file exists
test = FILE_TEST(files[0])
if (test eq 0) then begin
  error = 'REFERENCE FILE '+files[0]+' NOT FOUND'
  printlog,logf,error
  return
endif

; Load the reference data
LOADALS,files[0],refals,count=count
if (count lt 1) then begin
  error = 'PROBLEM LOADING '+files[0]
  printlog,logf,error
  return
endif


;; Get initial guess for X/Y shifts from WCS
;if keyword_set(initwcs) then begin
;  fitsfiles = file_basename(files,'.als')+'.fits'
;  if total(file_test(fitsfiles)) eq nfiles then begin
;    raarr = dblarr(nfiles) & decarr=dblarr(nfiles)
;    getpixscale,fitsfiles[0],pixscale
;    for i=0,nfiles-1 do begin
;      head = headfits(fitsfiles[i])
;      head_xyad,head,0,0,a,d,/deg
;      raarr[i]=a & decarr[i]=d
;    endfor
;    initwcs_xoff = (raarr-raarr[0])*3600*cos(decarr[0]/!radeg)/pixscale
;    initwcs_yoff = (decarr-decarr[0])*3600/pixscale
;    print,'Initial offsets from WCS'
;    for i=0,nfiles-1 do print,files[i],initwcs_xoff[i],initwcs_yoff[i]
;stop
;
;  endif else print,'Not all FITS files found'
;endif


; Use WCS
if keyword_set(usewcs) then begin
  ; Checking WCS of first file
  fitsfile1 = file_basename(files[0],'.als')+'.fits'
  if file_test(fitsfile1) eq 0 then fitsfile1=file_basename(files[0],'.als')+'.fits.fz'
  if file_test(fitsfile1) eq 0 then begin
    print,fitsfile1,' NOT FOUND. Cannot use WCS for matching'
  endif else begin
    if strmid(fitsfile1,6,7,/reverse_offset) eq 'fits.fz' then begin
      head1 = headfits(fitsfile1,exten=1)
      ; Fix the NAXIS1/2 in the header
      sxaddpar,head1,'NAXIS1',sxpar(head1,'ZNAXIS1')
      sxaddpar,head1,'NAXIS2',sxpar(head1,'ZNAXIS2')
    endif else begin
      head1 = headfits(fitsfile1)
    endelse
    EXTAST,head1,astr1,noparams1
    if noparams1 lt 1 then print,fitsfile1,' has NO WCS.  Cannot use WCS for matching'
  endelse
endif


format = '(A2,A-30,A1,2F10.2,4F10.5,2F10.3)'
newline = STRING("'",files[0],"'",0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, format=format)
PUSH,mchfinal,newline



if keyword_set(verbose) then $
  printlog,logf,'Initial Transformations:'

; Printing the first line
if keyword_set(verbose) then $
  printlog,logf,format='(A-20,2F10.4,4F12.8)',files[0], 0.0, 0.0, 1.0, 0.0, 0.0, 1.0



; Run DAOMATCH for each pair (N-1 times)
for i=1,nfiles-1 do begin

  undefine,als,alshead,trans,ind1,ind2,count

  ; Check that the file exists
  test = FILE_TEST(files[i])
  if (test eq 0) then begin
    error = 'FILE '+files[i]+' NOT FOUND'
    printlog,logf,error
    return
  endif

  ; Load the current data
  LOADALS,files[i],als,alshead,count=count
  if (count lt 1) then begin
    error = 'PROBLEM LOADING '+files[i]
    printlog,logf,error
    return
  endif


  ; Getting FRAD
  headarr = strsplit(alshead[1],' ',/extract)
  frad = float(first_el(headarr,/last))

  ; Make CHI, SHARP and ERR cuts here

  ; CUTS for REFALS
  gdref = where(abs(refals.sharp) lt 1.0 and refals.chi lt 2.0 and refals.mag lt 50.0 and $
                refals.err lt 0.2,ngdref)
  if (ngdref lt 100) then begin
    gdref = where(abs(refals.sharp) lt 1.5 and refals.chi lt 3.0 and refals.mag lt 50.0 and $
                  refals.err lt 0.5,ngdref)
  endif
  if (ngdref lt 100) then begin
    gdref = where(abs(refals.sharp) lt 1.5 and refals.chi lt 3.0 and refals.mag lt 50.0 and $
                  refals.err lt 1.0,ngdref)
  endif
  if (ngdref eq 0) then begin
    error = 'NO good reference stars '+files[0]
    printlog,logf,error
    return
  endif
  ; Cuts for ALS
  gdals = where(abs(als.sharp) lt 1.0 and als.chi lt 2.0 and als.mag lt 50.0 and $
                als.err lt 0.2,ngdals)
  if (ngdals lt 100) then begin
    gdals = where(abs(als.sharp) lt 1.5 and als.chi lt 3.0 and als.mag lt 50.0 and $
                  als.err lt 0.5,ngdals)
  endif
  if (ngdals lt 100) then begin
    gdals = where(abs(als.sharp) lt 1.5 and als.chi lt 3.0 and als.mag lt 50.0 and $
                  als.err lt 1.0,ngdals)
  endif
  if (ngdals eq 0) then begin
    error = 'NO good stars for '+files[i]
    printlog,logf,error
    return
  endif


  ; --- Use WCS ---
  if keyword_set(usewcs) and noparams1 ge 1 then begin
    ; Checking WCS of second file
    fitsfile2 = file_basename(files[i],'.als')+'.fits'
    if file_test(fitsfile2) eq 0 then fitsfile2 = file_basename(files[i],'.als')+'.fits.fz'
    if file_test(fitsfile2) eq 0 then begin
      printlog,logf,fitsfile2+' NOT FOUND. Cannot use WCS for matching'
      goto,BOMB1
    endif
    if strmid(fitsfile2,6,7,/reverse_offset) eq 'fits.fz' then begin
      head2 = headfits(fitsfile2,exten=1)
      ; Fix the NAXIS1/2 in the header
      sxaddpar,head2,'NAXIS1',sxpar(head2,'ZNAXIS1')
      sxaddpar,head2,'NAXIS2',sxpar(head2,'ZNAXIS2')      
    endif else begin
      head2 = headfits(fitsfile2)
    endelse
    EXTAST,head2,astr2,noparams2
    if noparams2 lt 1 then begin
      printlog,logf,fitsfile2+' has NO WCS.  Cannot use WCS for matching'
      goto,BOMB1
    endif

    ; Get coordinates for the stars
    head_xyad,head1,refals.x-1,refals.y-1,a1,d1,/deg
    head_xyad,head2,als.x-1,als.y-1,a2,d2,/deg

    SRCMATCH,a1[gdref],d1[gdref],a2[gdals],d2[gdals],1.0,ind1,ind2,count=count,/sph
  
    ; If no matches, try with looser cuts
    if count lt 3 then begin
      printlog,logf,'No matches, trying looser cuts'
      gdref = where(refals.mag lt 50.0 and refals.err lt 1.0,ngdref)
      gdals = where(als.mag lt 50.0 and als.err lt 1.0,ngdals)
      SRCMATCH,a1[gdref],d1[gdref],a2[gdals],d2[gdals],1.0,ind1,ind2,count=count,/sph
      if count lt 3 then $
        SRCMATCH,a1[gdref],d1[gdref],a2[gdals],d2[gdals],3.0,ind1,ind2,count=count,/sph
    endif

    if count gt 0 then begin 
      xdiff = refals[gdref[ind1]].x-als[gdals[ind2]].x
      ydiff = refals[gdref[ind1]].y-als[gdals[ind2]].y
      xmed = median([xdiff],/even)
      ymed = median([ydiff],/even)
      ; Fit rotation with linear fits if enough points
      if count gt 1 then begin
        coef1 = robust_poly_fitq(als[gdals[ind2]].y,xdiff,1)  ; fit rotation term
        coef1b = dln_poly_fit(als[gdals[ind2]].y,xdiff,1,measure_errors=xdiff*0+0.1,sigma=coef1err,/bootstrap)
        coef2 = robust_poly_fitq(als[gdals[ind2]].x,ydiff,1)  ; fit rotation term
        coef2b = dln_poly_fit(als[gdals[ind2]].x,ydiff,1,measure_errors=ydiff*0+0.1,sigma=coef2err,/bootstrap)
        ;theta = mean([-coef1[1],coef2[1]])
        WMEANERR,[-coef1[1],coef2[1]],[coef1err[1],coef2err[1]],theta,thetaerr

        ; [xoff, yoff, cos(th), sin(th), -sin(th), cos(th)]
        ;trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0]
        trans = [xmed, ymed, 1.0-theta^2, theta, -theta, 1.0-theta^2]
        ; Adjust Xoff, Yoff with this transformation
        xyout = trans_coo(als[gdals[ind2]].x,als[gdals[ind2]].y,trans)
        trans[0] += median([refals[gdref[ind1]].x - xyout[0,*]],/even) 
        trans[1] += median([refals[gdref[ind1]].y - xyout[1,*]],/even)
      endif else trans=[xmed, ymed, 1.0, 0.0, 0.0, 1.0]

      ; Fit full six parameters if there are enough stars
      if count gt 10 then begin
        fa = {x1:refals[gdref[ind1]].x,y1:refals[gdref[ind1]].y,x2:als[gdals[ind2]].x,y2:als[gdals[ind2]].y}
        ;initpar = fltarr(6)
        ;initpar = [xmed, ymed, 1.0, 0.0, 0.0, 1.0]
        ;initpar = [xmed, ymed, 1.0-theta^2, theta, -theta, 1.0-theta^2]
        initpar = trans
        fpar = MPFIT('trans_coo_dev',initpar,functargs=fa, perror=perror, niter=iter, status=status,$
                      bestnorm=chisq, dof=dof, autoderivative=1, /quiet) 
        trans = fpar
      endif
    endif
;print,files[i],' ',count
;stop
    BOMB1:
  endif ; use WCS

  ; Match stars with X/Y coordinates
  if (count lt 1) then begin
    ;MATCHSTARS,refals.x,refals.y,als.x,als.y,ind1,ind2,trans,count=count,/silent
    MATCHSTARS,refals[gdref].x,refals[gdref].y,als[gdals].x,als[gdals].y,ind1,ind2,trans,count=count,/silent
  endif

  ; No good matches, try srcmatch with "small" shifts
  if (count lt 1) then begin
    SRCMATCH,refals[gdref].x,refals[gdref].y,als[gdals].x,als[gdals].y,100,ind1a,ind2a,count=count1
    if count1 gt 0 then begin
      xdiff1 = refals[gdref[ind1a]].x-als[gdals[ind2a]].x
      ydiff1 = refals[gdref[ind1a]].y-als[gdals[ind2a]].y
      xmed1 = median([xdiff1],/even)
      ymed1 = median([ydiff1],/even)
      ; redo the search
      SRCMATCH,refals[gdref].x,refals[gdref].y,als[gdals].x+xmed1,als[gdals].y+ymed1,20,ind1,ind2,count=count
      if count eq 0 then begin
        SRCMATCH,refals[gdref].x,refals[gdref].y,als[gdals].x+xmed1,als[gdals].y+ymed1,100,ind1,ind2,count=count
      endif
      xdiff = refals[gdref[ind1]].x-als[gdals[ind2]].x
      ydiff = refals[gdref[ind1]].y-als[gdals[ind2]].y
      xmed = median([xdiff],/even)
      ymed = median([ydiff],/even)
      trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0]
    endif
  endif

  ; No good match
  if (count lt 1) then begin
    printlog,logf,'NO MATCHES.  Using XSHIFT=YSHIFT=ROTATION=0'
    trans = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
  endif

  ; Shift too large
  if keyword_set(maxshift) then begin
    if max(abs(trans[0:1])) gt maxshift then begin
      printlog,logf,'SHIFTS TOO LARGE. ',strtrim(trans[0:1],2),' > ',strtrim(maxshift,2),$
                    ' Using XSHIFT=YSHIFT=ROTATION=0'
      trans = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    endif
  endif

  ; The output is:
  ; filename, xshift, yshift, 4 trans, FRAD (from als file), 0.0
  format = '(A2,A-30,A1,2F10.2,4F10.5,2F10.3)'
  newline = STRING("'",files[i],"'",trans, frad, 0.0, format=format)
  PUSH,mchfinal,newline

  ; Printing the transformation
  if keyword_set(verbose) then $
    printlog,logf,format='(A-20,2F10.4,4F12.8)',files[i],trans

endfor  ; ALS file loop


; Writing the final mchfile
WRITELINE,files2[0]+'.mch',mchfinal

;stop



;#####################
; Running DAOMASTER
;#####################
rundaomaster:

; DAOMASTER has problems with files that have extra dots in them 
; (i.e. F1.obj1123_1.mch).
; Do everything with a temporary file, then rename the output files
; at the end.
;tempbase = MAKETEMP('temp','')
tempbase = FILE_BASENAME(MKTEMP('temp'))
FILE_DELETE,tempbase,/allow       ; remove empty file
tempbase = REPSTR(tempbase,'.','')   ; remove the dot
tempmch = tempbase+'.mch'
FILE_COPY,files2[0]+'.mch',tempmch,/overwrite,/allow

; Make the DAOMASTER script
;--------------------------
undefine,cmdlines
PUSH,cmdlines,'#!/bin/csh'
PUSH,cmdlines,'set input=${1}'
PUSH,cmdlines,'daomaster << DONE'
PUSH,cmdlines,'${input}.mch'
PUSH,cmdlines,'1,1,1'
PUSH,cmdlines,'99.'
PUSH,cmdlines,'6'
PUSH,cmdlines,'10'
PUSH,cmdlines,'5'
PUSH,cmdlines,'4'
PUSH,cmdlines,'3'
PUSH,cmdlines,'2'
PUSH,cmdlines,'1'
PUSH,cmdlines,'1'
PUSH,cmdlines,'1'
PUSH,cmdlines,'1'
PUSH,cmdlines,'0'
PUSH,cmdlines,'y'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'y'
PUSH,cmdlines,''
PUSH,cmdlines,'y'
PUSH,cmdlines,''
PUSH,cmdlines,''
PUSH,cmdlines,'y'
PUSH,cmdlines,''
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'DONE'
;tempscript = MAKETEMP('daomaster','.sh')
tempscript = MKTEMP('daomaster')   ; absolute filename
WRITELINE,tempscript,cmdlines
FILE_CHMOD,tempscript,'755'o

; Run DAOMASTER
;---------------
;cmd2 = '/net/home/dln5q/bin/daomaster.sh '+files2[0]
;cmd2 = './daomaster.sh '+tempbase
cmd2 = tempscript+' '+tempbase
SPAWN,cmd2,out2,errout2

;stop

; Remove temporary DAOMASTER script
;-----------------------------------
FILE_DELETE,tempscript,/allow_non


; Rename the outputs
;-------------------
; MCH file
mchfile = FILE_SEARCH(tempbase+'.mch',count=nmchfile)
if (nmchfile gt 0) then begin
  FILE_COPY,mchfile[0],files2[0]+'.mch',/overwrite,/allow
  FILE_DELETE,mchfile,/allow
endif else begin
  error = 'NO FINAL MCH FILE'
  printlog,logf,error
  return
endelse
; TFR file
tfrfile = FILE_SEARCH(tempbase+'.tfr',count=ntfrfile)
if (ntfrfile gt 0) then begin
  FILE_COPY,tfrfile[0],files2[0]+'.tfr',/overwrite,/allow
  FILE_DELETE,tfrfile,/allow
endif else begin
  error = 'NO FINAL TFR FILE'
  printlog,logf,error
  return
endelse
; RAW file
rawfile = FILE_SEARCH(tempbase+'.raw',count=nrawfile)
if (nrawfile gt 0) then begin
  FILE_COPY,rawfile[0],files2[0]+'.raw',/overwrite,/allow
  FILE_DELETE,rawfile,/allow
endif else begin
  error = 'NO FINAL RAW FILE'
  printlog,logf,error
  return
endelse

; FAKE, copy back the original MCH file
if keyword_set(fake) then begin
  FILE_COPY,files2[0]+'.mch',files2[0]+'.mch.daomaster',/overwrite,/allow
  FILE_MOVE,files2[0]+'.mch.orig',files2[0]+'.mch',/overwrite,/allow
endif

; Were there any errors
if (errout2[0] ne '') then begin
  printlog,logf,'DAOMASTER.SH ERROR'
  printlog,logf,errout2
  error = errout2
  return
endif


; Print out the final transformations
if keyword_set(verbose) then begin
  LOADMCH,files2[0]+'.mch',files,trans
  nfiles = n_elements(files)
  printlog,logf,''
  printlog,logf,'Final DAOMASTER Transformations:'
  for i=0,nfiles-1 do begin
    printlog,logf,format='(A-20,2F10.4,4F12.8)',files[i],transpose(trans[i,*])
  end
endif


; Back to the original directory
CD,curdir

if keyword_set(stp) then stop

end
