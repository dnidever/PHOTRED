pro stdred_matchcat,redo=redo,stp=stp

;+
;
; STDRED_MATCHCAT
;
; This gets the calibrated standard star photometry
; for standard stars in the each image.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  CAT files for each frame that contains some standard
;  stars.  The CAT files has the same columns as the
;  AST files plus the calibrated data.
;
; By D.Nidever  May 2008
;-

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING STDRED_MATCHCAT'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'MATCHCAT'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting STDRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','add_tag','printstr','photred_getinput',$
         'photred_updatelists','photred_loadsetup','push','undefine','printlog','srcmatch',$
         'first_el','minloc','range','strsplitter','touchzero','writeline','mktemp','stress',$
         'strep']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count,/std
  if count lt 1 then return
endif


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')
matchdist = READPAR(setup,'MATCHDIST')
if matchdist eq '0' or matchdist eq '-1' or matchdist eq '' then $
  undefine,matchdist else matchdist=float(matchdist)


;###################
; GETTING INPUTLIST
;###################
; INLIST         PHOT files
; OUTLIST        AST files
; SUCCESSLIST    PHOT files

; Get input
;-----------
precursor = 'ASTROM'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='ast')
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines




; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif


; Get the list of standard star fields
; Check for a local file first
printlog,logfile,'Loading the >>standards.lst<< file'
stdtest = FILE_TEST('standards.lst')

; NO Local "standards.lst" file
if (stdtest eq 0) then begin

  printlog,logfile,'NO local >>standards.lst<< file.  Checking '+scriptsdir
  stdtest = FILE_TEST(scriptsdir+'/standards.lst')
  if (stdtest eq 0) then begin
    printlog,logfile,'NO >>standards.lst<< in ',scriptsdir
    return
  endif
 
  ; Copy to current directory
  printlog,logfile,'Found '+scriptsdir+'/standards.lst.  Copying to current directory.'
  FILE_COPY,scriptsdir+'/standards.lst',curdir,/allow,/over

endif

; Now read in the standards.lst file
READLINE,'standards.lst',stdfiles1,count=nstdfiles1,comment='#',/noblank
if (nstdfiles1 eq 0) then begin
  printlog,logfile,'NO standard star catalog file names in >>standards.lst<<'
  return
endif
ui = uniq(stdfiles1,sort(stdfiles1))  ; getting unique ones
stdfiles1 = stdfiles1[ui]

; Copy ALL the standard star catalogs to the current directory
printlog,logfile,''
stdtest2 = FILE_TEST(scriptsdir+'/'+stdfiles1)
bd = where(stdtest2 eq 0,nbd)
if (nbd gt 0) then begin
  printlog,logfile,stdfiles1[bd],' NOT FOUND'
endif
gd = where(stdtest2 eq 1,ngd)
; Some good ones
if (ngd gt 0) then begin
  printlog,logfile,strtrim(ngd,2)+' standard star catalogs found'
  printlog,logfile,'Copying '+stdfiles1[gd]
  FILE_COPY,scriptsdir+'/'+stdfiles1[gd],'.',/overwrite,/allow_same
  stdfiles = stdfiles1[gd]

; No standard star catalog files
endif else begin
  printlog,logfile,'NO standard star catalog files'
  return
endelse



; Restoring the standard star catalogs
; Naming the structures 'std1', 'std2', etc.
nstdfiles = n_elements(stdfiles)
;stdnames = strarr(nstdfiles)        ; names of the fields, i.e. SA98, etc.
;stdvarnames = strarr(nstdfiles)     ; names of structures, i.e. std1, std2, ...
;stdra = fltarr(nstdfiles)           ; RA center of fields
;stddec = fltarr(nstdfiles)          ; DEC center of fields
stdstr = replicate({name:'',varname:'',ra:0.0d0,dec:0.0,size:0.0},nstdfiles)

for i=0,nstdfiles-1 do begin
  dum = strsplit(stdfiles[i],'.',/extract)
  ext = first_el(dum,/last)
  if ext eq 'gz' then begin
    pos = strpos(stdfiles[i],'.')
    ext = strmid(stdfiles[i],pos+1)
  endif
  case ext of
  'cat': std = IMPORTASCII(stdfiles[i],/header,/noprint)
  'fits': std = MRDFITS(stdfiles[i],1,/silent)
  'fits.gz': std = MRDFITS(stdfiles[i],1,/silent)
  else: stop,'Extension ',ext,' NOT SUPPORTED'
  endcase

  fname = FILE_BASENAME(stdfiles[i],'.cat')
  ;(SCOPE_VARFETCH(fname+'cat',/enter)) = std
  ;;stdnames[i] = fname+'cat'
  ;stdnames[i] = fname
  ;stdvarnames[i] = 'std'+strtrim(i+1,2)
  stdstr[i].name = fname
  stdstr[i].varname = 'std'+strtrim(i+1,2)
  (SCOPE_VARFETCH(stdstr[i].varname,/enter)) = std

  ; Get center of standard fields
  ;  check if the fields wraps at RA=0
  hi = where(std.ra gt 350,nhi)
  lo = where(std.ra lt 10,nlo)
  if nlo gt 0 and nhi gt 0 then begin
    ra = std.ra
    ra[hi] -= 360
    medra = median(ra)
    if medra lt 0 then medra+=360
    ;stdra[i] = medra
    stdstr[i].ra = medra
  endif else begin
    ;stdra[i] = median(std.ra)
    stdstr[i].ra = median(std.ra)
  endelse
  ;stddec[i] = median(std.dec)
  stdstr[i].dec = median(std.dec)
  dist = sphdist(stdstr[i].ra,stdstr[i].dec,std.ra,std.dec,/deg)
  stdstr[i].size = max(dist)

  print,stdstr[i].name,stdstr[i].ra,stdstr[i].dec

endfor


;##################################################
;#  PROCESSING THE FILES
;##################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'

undefine,outlist,successlist,failurelist

; Loop through the input PHOT files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  base = FILE_BASENAME(file,'.ast')
  filedir = FILE_DIRNAME(longfile)
  fitsfile = base+'.fits'

  printlog,logfile,''
  printlog,logfile,'================================================='
  printlog,logfile,'Getting calibrated photometry for ',file
  printlog,logfile,'================================================='


  CD,filedir

  ; Test that the AST file exists
  test = FILE_TEST(file)
  if test eq 1 then nlines = FILE_LINES(file)
  if (test eq 0) or (nlines eq 0) then begin
    PUSH,failurelist,longfile
    if test eq 0 then printlog,logfile,file,' NOT FOUND'
    if test eq 1 and nlines eq 0 then printlog,logfile,file,' HAS 0 LINES'
    goto,BOMB
  endif

  ; Loading the AST file
  ;---------------------
  ; MAGFAP is the aperture magnitude for the final aperture for this star
  ; APCORR is the aperture correction for MAGFAP
  ; FINALAP is the number for the final aperture for this star
  ;            i.e. 5 for the 5th aperture
  ; MAG = MAGFAP + APCORR
  phot = IMPORTASCII(file,/header,/noprint)
  nphot = n_elements(phot)
  printlog,logfile,'Nstars = ',strtrim(nphot,2)

  medra = median(phot.ra)
  meddec = median(phot.dec)

  ; Figure out if this is within 2 degree of any of the standard star fields
  dist = sphdist(medra,meddec,stdstr.ra,stdstr.dec,/deg)
  mindist = min(dist)
  minind = first_el(minloc(dist))
  dist_thresh = 2.0 > stdstr[minind].size

  ; This frame is of a standard field
  if (mindist lt dist_thresh) then begin

    printlog,logfile,'This frame is of standard star field: ',stdstr[minind].name

    ; Get the standard star structure
    std = (SCOPE_VARFETCH(stdstr[minind].varname))
    nstd = n_elements(std)

    ; Now match the stars
    dcr = 0.8
    if keyword_set(matchdist) then dcr = matchdist
    dcr = (dcr > 0.2) < 2.0    ; keep it reasonable
    SRCMATCH,phot.ra,phot.dec,std.ra,std.dec,dcr,ind1,ind2,/sph,count=count

    ;print,sphdist(phot[ind1].ra,phot[ind1].dec,std[ind2].ra,std[ind2].dec,/deg)*3600.

    printlog,logfile,'There are ',strtrim(count,2),'/',strtrim(nstd,2),' matches for the standard stars'

    ; Add the calibrated photometry to the structure
    if (count gt 0) then begin

      ; Matched catalogs
      photmatch = phot[ind1]
      stdmatch = std[ind2]

      std_tags = TAG_NAMES(std)
      nstd_tags = N_TAGS(std)
      types = intarr(nstd_tags)
      for k=0,nstd_tags-1 do types[k]=size(std[0].(k),/type)
      new_tags = 'C_'+std_tags

      ; Make new structure
      dum = photmatch[0]
      for k=0,nstd_tags-1 do dum=CREATE_STRUCT(dum,new_tags[k],fix(0,type=types[k]))
      cat = REPLICATE(dum,count)
      ncat_tags = N_TAGS(cat)

      ; Copy in the instrumental data
      phot_tags = TAG_NAMES(photmatch)
      nphot_tags = N_TAGS(photmatch)
      for k=0,nphot_tags-1 do cat.(k)=photmatch.(k)

      ; Copy in the calibrated data
      for k=0,nstd_tags-1 do cat.(k+nphot_tags)=stdmatch.(k)

      ; Now save the matched catalog
      catfile = filedir+'/'+base+'.cat'
      printlog,logfile,'Saving matched catalog to ',catfile
      PRINTSTR,cat,catfile,/silent

      ; Success
      testcat = FILE_TEST(catfile)
      if (testcat eq 1) then begin
        PUSH,successlist,longfile
        PUSH,outlist,catfile
      endif else begin
        printlog,logfile,catfile,' NOT FOUND'
        PUSH,failurelist,longfile
      endelse

    ; No standard star matches for this frame
    endif else begin

      printlog,logfile,'NO MATCHES for ',file
      PUSH,failurelist,longfile

    endelse


  ; This frame is NOT of a standard field
  endif else begin

    printlog,logfile,'This frame is NOT of a standard star field'
    PUSH,failurelist,longfile

  endelse


  BOMB:

  CD,curdir


  ;##########################################
  ;#  UPDATING LIST FILES
  ;##########################################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent

END


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'STDRED_MATCHCAT Finished  ',systime(0)

if keyword_set(stp) then stop




end
