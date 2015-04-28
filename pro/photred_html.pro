pro photred_html,redo=redo,stp=stp

;+
;
; PHOTRED_HTML
;
; This looks at all of the final FIELD files and makes a separate
; HTML page for each as well as an INDEX page.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  HTML Summary pages
;
; By D.Nidever  June 2008
;-

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_HTML'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'HTML'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','photred_getinput','photred_updatelists','photred_loadsetup',$
         'photred_getairmass','photred_getdate','photred_getexptime','photred_getfilter','photred_getuttime',$
         'airmass','maxloc','remove','printlog','push','range','ten2sexig','sexig2ten',$
         'loadals','importascii','minloc','mad','scale_vector','loadmch','mktemp','loadinput',$
         'undefine','first_el','strsplitter','touchzero','writeline','scale','stringize','signs','stress',$
         'strep','strmult','strtrim0','sign','display','imgscl','displayc','colorbar','tv_24',$
         'psym8','ps_open','ps_close','ps2gif','closest','odd','zscale','badpar','goodpoly']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif

; LOAD THE SETUP FILE
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
PHOTRED_LOADSETUP,setup,count=count
if count lt 1 then begin
  print,'NO >>photred.setup<< file'
  return
endif


; LOAD information from the "photred.setup" file
;-----------------------------------------------
; REDO
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1
; TELESCOPE
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(strtrim(telescope,2))
if (telescope eq '0' or telescope eq '' or telescope eq '-1') then begin
  printlog,logfile,'NO TELESCOPE FOUND.  Please add to >>photred.setup<< file'
  return
endif
; INSTRUMENT
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(strtrim(instrument,2))
if (instrument eq '0' or instrument eq '' or instrument eq '-1') then begin
  printlog,logfile,'NO INSTRUMENT FOUND.  Please add to >>photred.setup<< file'
  return
endif
; CMD2CDAXES
cmd2cdaxes = READPAR(setup,'CMD2CDAXES')
cmd2cdaxes = strupcase(strtrim(cmd2cdaxes,2))
if (cmd2cdaxes eq '0' or cmd2cdaxes eq '' or cmd2cdaxes eq '-1') then begin
  printlog,logfile,'CMD2CDAXES NOT INPUT.  Please add to >>photred.setup<< file'
  printlog,logfile,"For example, CMD2CDAXES  'M,M-T,M-D' will produce M vs. M-T CMD and M-D vs. M-T 2CD"
  return
endif

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif


; LOAD THE "imagers" FILE
;----------------------------
printlog,logfile,'Loading imager information'
imagerstest = FILE_TEST(scriptsdir+'/imagers')
if (imagerstest eq 0) then begin
  printlog,logfile,'NO >>imagers<< file in '+scriptsdir+'  PLEASE CREATE ONE!'
  return
endif
; The columns need to be: Telescope, Instrument, Naps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(scriptsdir+'/imagers',fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#')
imagers.telescope = strupcase(strtrim(imagers.telescope,2))
imagers.instrument = strupcase(strtrim(imagers.instrument,2))
imagers.observatory = strupcase(strtrim(imagers.observatory,2))
singleind = where(imagers.namps eq 1,nsingle)
if nsingle gt 0 then imagers[singleind].separator = ''
if (n_tags(imagers) eq 0) then begin
  printlog,logfile,'NO imagers in '+scriptsdir+'/imagers'
  return
endif



; What IMAGER are we using??
;---------------------------
ind_imager = where(imagers.telescope eq telescope and imagers.instrument eq instrument,nind_imager)
if nind_imager eq 0 then begin
  printlog,logfile,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
  return
endif
thisimager = imagers[ind_imager[0]]
; print out imager info
printlog,logfile,''
printlog,logfile,'USING IMAGER:'
printlog,logfile,'Telescope = '+thisimager.telescope
printlog,logfile,'Instrument = '+thisimager.instrument
printlog,logfile,'Namps = '+strtrim(thisimager.namps,2)
printlog,logfile,"Separator = '"+thisimager.separator+"'"
printlog,logfile,''

; How many amps are there?
namps = thisimager.namps




;###################
; GETTING INPUTLIST
;###################
; INLIST         DAT files
; OUTLIST        HTML files
; SUCCESSLIST    DAT files


; Get input
;-----------
precursor = 'SAVE'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='dat',/noempty)
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines





; Make the HTML directory
if FILE_TEST('html',/directory) eq 0 then $
  FILE_MKDIR,'html'

printlog,logfile,'Making SUMMARY HTML pages for ',strtrim(ninputlines,2),' Pointings/Fields'
printlog,logfile,''

; Load the "fields" file
READCOL,'fields',shfields,fields,format='A,A',/silent

; Load the "apcor.lst" file
; Maybe the CALIB input files instead
apcor = IMPORTASCII('apcor.lst',fieldnames=['name','value'],/noprint)
; Remove the 'a.del' endings for the names
apcor_orig = apcor
apcor.name = repstr(apcor.name,'a.del','')  ; base names


; Magnitudes and colors to use for CMD and 2CD
;----------------------------------------------
; CMD2CDAXES = 'M,M-T,M-D'
; gives CMD of M vs. M-T
;   and 2CD of M-D vs. M-T
magname = ''
col1name = ['','']
col2name = ['','']
if n_elements(cmd2cdaxes) gt 0 then begin
  axes = strsplit(cmd2cdaxes,',',/extract)
  naxes = n_elements(axes)

  ; Need at least 2 axes
  if naxes lt 2 then begin
    printlog,logfile,'ONLY '+strtrim(naxes,2)+' Color/Magnitude for CMD/2CD.  Need at LEAST 2'
    goto,AXESBOMB
  endif

  ; First axis must NOT be a color
  if (strpos(axes[0],'-') ne -1) then begin
    printlog,logfile,'ERROR '+axes[0]+' FIRST CMD/2CD axis CANNOT be a COLOR and MUST be MAGNITUDE'
    goto,AXESBOMB
  endif

  ; Second axis MUST be a color
  if (strpos(axes[1],'-') eq -1) then begin
    printlog,logfile,'ERROR '+axes[1]+' SECOND CMD/2CD axis MUST a COLOR, i.e. M-T'
    goto,AXESBOMB
  endif

  ; NO 3rd color, CAN'T make 2CD
  if (naxes lt 3) then begin
    printlog,logfile,'ONLY 2 CMD/2CD AXES. CANNOT make 2CD'
  endif

  ; Magnitude name
  magname = strtrim(axes[0],2)
    
  ; First color
  col1axis = axes[1]
  col1axisarr = strsplit(col1axis,'-',/extract)
  col1name = [col1axisarr[0],col1axisarr[1]]
  col1name = strtrim(col1name,2)

  ; Incomplete color
  if (col1name[0] eq '' or col1name[1] eq '') then begin
    printlog,logfile,col1axis+' is an INCOMPLETE COLOR.  CANNOT make CMD/2CD'
    col1name=['','']
    goto,AXESBOMB
  endif


  ; Second color
  if (naxes ge 3) then begin
    col2axis = axes[2]

    ; This is a color
    if (strpos(col2axis,'-') ne -1) then begin
      col2axisarr = strsplit(col2axis,'-',/extract)
      col2name = [col2axisarr[0],col2axisarr[1]]
      col2name = strtrim(col2name,2)

      ; Incomplete color
      if (col2name[0] eq '' or col2name[1] eq '') then begin
        printlog,logfile,col2axis+' is an INCOMPLETE COLOR.  CANNOT make CMD/2CD'
        col2name=['','']
        goto,AXESBOMB
      endif

    ; NOT a color
    endif else begin
      printlog,logfile,'The 3rd CMD/2CD axis is NOT a COLOR.  CANNOT make 2CD'
    endelse

  endif
endif
AXESBOMB:


;; Load the black-white color table
;loadcol,3
loadct,0,/silent
psym8

fieldarr = strarr(ninputlines)
htmlfilearr = strarr(ninputlines)
undefine,allinfo

; Start the field info structure
infostruct = {base:'',field:'',uttime:'',date:'',airmass:0.0,nsources:0L,nstars:0L,fwhm:0.0,type:'',$
              ra:0.0d0,dec:0.0d0,cmdfile:'',imagefile:'',skyfile:'',histfile:'',errorfile:'',$
              chisharpfile:'',chimagfile:'',cmdcolampfile:'',cmdgiantsfile:'',skygiantsfile:'',$
              htmlfile:''}
fieldinfo = REPLICATE(infostruct,ninputlines)

undefine,outlist,successlist,failurelist,allinfo

FOR i=0,ninputlines-1 do begin
;FOR i=0,1 do begin

  idatfile = inputlines[i]
  ifield = FILE_BASENAME(idatfile,'.dat')
  fieldinfo[i].field = ifield

  ; Get the short field name
  gdsh = where(fields eq ifield,ngdsh)
  if (ngdsh eq 0) then begin
    printlog,logfile,ifield+' is not in >>fields<<'
    PUSH,failurelist,inputlines[i]
    goto,FIELD_BOMB
  endif
  ishortfield = shfields[gdsh[0]]

  ; Get the DERED file
  deredfile = FILE_SEARCH(ishortfield+'-'+'*.dered',count=nderedfile)
  if (nderedfile eq 0) then begin
    printlog,logfile,'NO DERED file for '+ifield
    PUSH,failurelist,inputlines[i]
    goto,FIELD_BOMB
  endif

  ; Main BASE name, i.e. F10-obj1108
  base = FILE_BASENAME(deredfile,'.dered')
  base = base[0]
  fieldinfo[i].base = base

  printlog,logfile,'--- MAKING HTML SUMMARY PAGE for '+base+' - '+ifield+' ---'


  ; Load the DAT file
  dattest = FILE_TEST(idatfile)
  if (dattest eq 0) then begin
    printlog,logfile,idatfile+' NOT FOUND'
    PUSH,failurelist,inputlines[i]
    goto,FIELD_BOMB
  endif
  RESTORE,idatfile               ; the structure name is "final"
  nfinal = n_elements(final)


  ; Getting the Magnitude/Error field indices
  ;------------------------------------------
  ; NO magindarr input
  tags = tag_names(final)

  ; Assuming that the file is in the format:
  ; ID, X, Y, MAG1, ERR1, MAG2, ERR2, ..., CHI, SHARP, other tags
  gdy = where(tags eq 'Y',ngdy)
  if (ngdy eq 0) then begin
    printlog,logfile,'No "Y" TAG FOUND'
    error = 'No "Y" TAG FOUND'
    return
  endif
  gdchi = where(tags eq 'CHI',ngdchi)
  if (ngdchi eq 0) then begin
    printlog,logfile,'NO "CHI" TAG FOUND'
    error = 'NO "CHI" TAG FOUND'
    return
  endif
  lo = gdy[0] + 1         ; first magnitude field is AFTER Y
  hi = gdchi[0] - 1       ; last error field is BEFORE CHI

  ; Need an even number of magnitude/error fields
  if odd(hi-lo+1) eq 1 then begin
    printlog,logfile,'NOT THE RIGHT NUMBER OF MAGNITUDE/ERROR TAGS BETWEEN "Y" AND "CHI"'
    error = 'NOT THE RIGHT NUMBER OF MAGNITUDE/ERROR TAGS BETWEEN "Y" AND "CHI"'
    return
  endif
  nfilters = (hi-lo+1)/2
  magindarr = lindgen(nfilters)*2+lo      ; indices for the magnitudes fields
  nmagindarr = n_elements(magindarr)
  errindarr = magindarr+1                  ; indices for the error fields
  filters = tags[magindarr]

  ; Bands to use for the CMD and 2CD
  ;---------------------------------
  ; Do we have AXES information
  nocmd=0 & no2cd=0
  if magname eq '' or col1name[0] eq '' or col1name[1] eq '' then begin
    nocmd=1 & no2cd=1
  endif
  if col2name[0] eq '' or col2name[1] eq '' then no2cd=1
  ; Get Magnitude index in the structure
  if not keyword_set(nocmd) then begin
    magfilterind = where(filters eq magname,nmagnameind)
    if (nmagnameind eq 0) then begin
      printlog,logfile,magname+' NOT FOUND. CANNOT MAKE CMD/2CD'
      nocmd=1 & no2cd=1
      magind = -1
    endif else magind = magindarr[magfilterind[0]]
  endif
  ; Get Color1 index in the structure
  if not keyword_set(nocmd) then begin
    col1mag1filterind = where(filters eq col1name[0],ncol1mag1filterind)
    if (ncol1mag1filterind eq 0) then begin
      printlog,logfile,col1name[0]+' NOT FOUND. CANNOT MAKE CMD/2CD'
      nocmd=1 & no2cd=1
      col1mag1ind = -1
    endif else col1mag1ind = magindarr[col1mag1filterind[0]]
    col1mag2filterind = where(filters eq col1name[1],ncol1mag2filterind)
    if (ncol1mag2filterind eq 0) then begin
      printlog,logfile,col1name[1]+' NOT FOUND. CANNOT MAKE CMD/2CD'
      nocmd=1 & no2cd=1
      col1mag2ind = -1
    endif else col1mag2ind = magindarr[col1mag2filterind[0]]
  endif
  ; Get Color2 index in the structure
  if not keyword_set(nocmd) and not keyword_set(no2cd) then begin
    col2mag1filterind = where(filters eq col2name[0],ncol2mag1filterind)
    if (ncol2mag1filterind eq 0) then begin
      printlog,logfile,col2name[0]+' NOT FOUND. CANNOT MAKE 2CD'
      no2cd=1   
      col2mag1ind = -1
    endif else col2mag1ind = magindarr[col2mag1filterind[0]]
    col2mag2filterind = where(filters eq col2name[1],ncol2mag2filterind)
    if (ncol2mag2filterind eq 0) then begin
      printlog,logfile,col2name[1]+' NOT FOUND. CANNOT MAKE 2CD'
      no2cd=1
      col2mag2ind = -1
    endif else col2mag2ind = magindarr[col2mag2filterind[0]]
  endif
  ; Get ERROR index in structure
  if not keyword_set(nocmd) then begin
    magerrind = where(tags eq strupcase(magname)+'ERR',nmagerrind)
  endif

  ; NO Magnitude index
  ; Need these for plots that need a magnitude
  if (magind eq -1) then begin
    printlog,logfile,'NO input magnitude to work with.  Using FIRST magnitude - '+filters[0]
    magind = magindarr[0]
    magerrind = magind+1
  endif
  ; NO color1mag1 index
  ; Need this for selecting "good" stars
  if (col1mag1ind eq -1) then begin
    printlog,logfile,'NO input color to work with.  Using FIRST magnitude COL1MAG1- '+filters[0]
    magind = magindarr[0]
  endif
  ; NO color1mag2 index
  ; Need this for selecting "good" stars
  ; THIS WON'T WORK IF THERE IS ONLY *ONE* FILTER INDEX!!!
  if (col1mag2ind eq -1) then begin
    printlog,logfile,'NO input color to work with.  Using SECOND magnitude COL1MAG2 - '+filters[1]
    magind = magindarr[1]
  endif

  ; Only stars with decent photometry
  gdfinal = where(final.(magind) lt 50. and final.(col1mag1ind) lt 50. and $
                  final.(col1mag2ind) lt 50. and abs(final.sharp) lt 1.0,ngdfinal)
  final2 = final[gdfinal]

  ; Can we select GIANTS?
  Mfilterind = where(filters eq 'M',nMfilterind)
  Tfilterind = where(filters eq 'T',nTfilterind)
  Dfilterind = where(filters eq 'D',nDfilterind)
  candogiants=0
  if nMfilterind gt 0 and nTfilterind gt 0 and nDfilterind gt 0 then candogiants=1
  ; WE HAVE Washington+DDO51 photometry, CAN select giants
  if (candogiants eq 1) then begin
    printlog,logfile,ifield+' has WASHINGTON+DDO51 photometry. CAN select giants'
    ; Use dereddened magnitudes if possible
    M0filterind = where(filters eq 'M0',nM0filterind)
    T0filterind2 = where(filters eq 'T0',nT0filterind)
    D0filterind2 = where(filters eq 'D0',nD0filterind)
    ; Using dereddened magnitudes
    if (nM0filterind gt 0 and nT0filterind gt 0 and nD0filterind gt 0) then begin
      Mmagind = magindarr[M0filterind[0]]
      Tmagind = magindarr[T0filterind[0]]
      Dmagind = magindarr[D0filterind[0]]
    ; Using undereddened magnitudes
    endif else begin
      Mmagind = magindarr[Mfilterind[0]]
      Tmagind = magindarr[Tfilterind[0]]
      Dmagind = magindarr[Dfilterind[0]]
    endelse
  ; NO Washington+DDO51 photometry, can't select giants
  endif else begin
    if nMfilterind eq 0 then PUSH,missingfilters,'M'
    if nTfilterind eq 0 then PUSH,missingfilters,'T'
    if nDfilterind eq 0 then PUSH,missingfilters,'D'
    outmissingfilters = strjoin(missingfilters,',')
    printlog,logfile,'CANNOT SELECT GIANTS. MISSING '+outmissingfilters+' BANDS'
  endelse


  ; Some field info
  fieldinfo[i].ra = mean(final.ra)
  fieldinfo[i].dec = mean(final.dec)
  fieldinfo[i].nsources = nfinal
  fieldinfo[i].nstars = ngdfinal


  ; Make fieldampinfo structure
  ampinfostruct = {amp:0,base:'',field:'',uttime:'',date:'',airmass:0.0,nsources:0L,nstars:0L,fwhm:0.0,type:'',$
                   ra:0.0d0,dec:0.0d0,cmdfile:'',imagefile:'',skyfile:'',histfile:'',errorfile:'',$
                   chisharpfile:'',chimagfile:'',cmdcolampfile:'',cmdgiantsfile:'',skygiantsfile:'',$
                   htmlfile:''}
  fieldampinfo = REPLICATE(ampinfostruct,namps)

  ;--------------------------------
  ; Looping through the Amplifiers
  ;--------------------------------
  For j=0,namps-1 do begin

    iamp = j+1

    fieldampinfo[j].amp = iamp

    ; Get the MCH file
    if (namps gt 1) then begin
      printlog,logfile,'AMP='+strtrim(iamp,2)
      ampbase = base+thisimager.separator+strtrim(iamp,2)
      mchfile = ampbase+'.mch'
    endif else begin
      ampbase = base
      mchfile = base+'.mch'
    endelse

    mchtest = FILE_TEST(mchfile)
    if mchtest eq 0 then begin
      printlog,logfile,mchfile+' NOT FOUND'
      PUSH,failurelist,inputlines[i]
      goto,AMP_BOMB
    endif

    ; Load the MCH file
    LOADMCH,mchfile,alsfiles,trans,count=nalsfiles
    alsbases = FILE_BASENAME(alsfiles,'.als')
    printlog,logfile,strtrim(nalsfiles,2)+' ALS files'

    ; Start the frame info structure
    dum = {filename:'',base:'',field:'',object:'',nx:0L,ny:0L,filter:'',exptime:0.0,uttime:'',$
           date:'',airmass:0.0,skymed:0.0,skysig:0.0,wcs:0,wcstype:'',wcsrms:0.0,ra:0.0d0,$
           dec:0.0d0,trans:'',nsources:0L,nstars:0L,psfchi:'',$
           fwhm:0.0,saturation:0.0,apcor:0.0,imagefile:'',skyfile:'',histfile:'',errmagfile:'',$
           chisharpfile:'',chimagfile:''}
    info = REPLICATE(dum,nalsfiles)
    
    ;---------------------------------------
    ; Loop through the ALS separate files
    ;---------------------------------------
    For k=0,nalsfiles-1 do begin

      ialsbase = alsbases[k]
      printlog,logfile,ialsbase

      ; What is the base?
      if (namps gt 1) then begin
        ; remove the amp ending 
        dum = strsplit(ialsbase,thisimager.separator,/extract)
        thisbase = dum[0]
        info[k].base = thisbase
      endif else begin
        info[k].base = ialsbase
      endelse

      ; Load FITS file
      fitsfile = ialsbase+'.fits'
      FITS_READ,fitsfile,im,head
      sz = size(im)
      nx = sz[1] & ny = sz[2]

      ; Save the HEADER
      WRITELINE,'html/'+ialsbase+'.header',head

      ; Load ALS file
      LOADALS,ialsbase+'.als',als
      nals = n_elements(als)
      info[k].nsources = nals
      gdstars = where(als.chi lt 2.0 and abs(als.sharp) lt 1.0 and als.mag lt 50.,nstars)
      info[k].nstars = nstars

      ; Getting the FILE information
      filter = PHOTRED_GETFILTER(fitsfile)
      exptime = PHOTRED_GETEXPTIME(fitsfile)
      uttime = PHOTRED_GETUTTIME(fitsfile)
      date = PHOTRED_GETDATE(fitsfile)
      airmass = PHOTRED_GETAIRMASS(fitsfile)
      object = SXPAR(head,'OBJECT',/silent)
      info[k].filename = ialsbase
      info[k].field = ifield
      info[k].object = object
      info[k].nx = nx
      info[k].ny = ny
      info[k].filter = filter
      info[k].exptime = exptime
      info[k].uttime = uttime
      info[k].date = date
      info[k].airmass = airmass

      ; Image statistics
      info[k].skymed = median(im,/even)
      info[k].skysig = mad(im)

      ; WCS stuff
      EXTAST,head,astro
      info[k].wcstype = '--'
      info[k].ra = 99.9999
      info[k].dec = 99.9999
      if n_elements(astro) gt 0 then begin
        info[k].wcs=1
        ctype = astro[0].ctype[0]
        wcstype = first_el(strsplit(ctype,'-',/extract),/last)
        info[k].wcstype = wcstype
        info[k].ra = astro.crval[0]
        info[k].dec = astro.crval[1]
      endif
      gwcs = where(stregex(head,'WCS',/boolean) eq 1 and $
                   stregex(head,'RMS',/boolean) eq 1,ngwcs)
      gwcs = first_el(gwcs,/last)
      info[k].wcsrms = 99.99
      if ngwcs gt 0 then begin
        wcsarr = strsplit(head[gwcs[0]],' ',/extract)
        strwcs = wcsarr[2]
        wcsrms = first_el(strsplit(strwcs,'=',/extract),/last)
        wcsrms = float(wcsrms)
        info[k].wcsrms = wcsrms
      endif

      ; Load the OPT file
      info[k].fwhm = 99.99
      info[k].saturation = 999999.
      if FILE_TEST(ialsbase+'.opt') eq 1 then begin
        READLINE,ialsbase+'.opt',optlines
        gfwhm = where(stregex(optlines,'^FW = ',/boolean) eq 1,ngfwhm)
        if ngfwhm gt 0 then begin
          fwhmarr = strsplit(optlines[gfwhm[0]],' ',/extract)
          fwhm = float(fwhmarr[2])
          info[k].fwhm = fwhm
        endif
        gsat = where(stregex(optlines,'^HI = ',/boolean) eq 1,ngsat)
        if ngsat gt 0 then begin
          satarr = strsplit(optlines[gsat[0]],' ',/extract)
          sat = float(satarr[2])
          info[k].saturation = sat
        endif
      endif

      ; Get transformation equations
      itrans = trans[k,*]
      info[k].trans = string(itrans,format='(2F8.2,4F8.4)')

      ; Aperture Correction
      apind = where(apcor.name eq ialsbase,napind)
      info[k].apcor = 99.99
      if napind gt 0 then $
        info[k].apcor = apcor[apind[0]].value


      ; Make an IMAGE of the FITS file
      ;-------------------------------
      ; USE MSCDISPLAY.PRO FOR MOSAIC, IMACS OR LBT!!!
      ; The larger size should be ~10 inches
      imagefile1 = ialsbase+'_image.gif'
      test = FILE_TEST('html/'+imagefile1)
      if test eq 0 or keyword_set(redo) then begin
        fac = 10./float(nx) < 10./float(ny)
        xsize = nx*fac
        ysize = ny*fac + 1.5     ; plus 1 for the colorbar
        ps_open,'html/'+ialsbase+'_image',thick=3,/portrait
        device,xsize=xsize,ysize=ysize,/inches
        displayc,im,/z,xtit='X',ytit='Y',tit=fitsfile,charsize=1.1
        ps_close
        ps2gif,'html/'+ialsbase+'_image.ps'
        FILE_DELETE,'html/'+ialsbase+'_image.ps',/allow,/quiet
      endif
      info[k].imagefile = imagefile1

      ; Spatial distribution
      ;----------------------------
      skyfile1 = ialsbase+'_sky.gif'
      test = FILE_TEST('html/'+skyfile1)
      if test eq 0 or keyword_set(redo) then begin
        ps_open,'html/'+ialsbase+'_sky',thick=5
        plot,als.x,als.y,ps=8,sym=0.5,xr=minmax(als.x),yr=minmax(als.y),xs=1,ys=1,$
                 xtit='X (pixels)',ytit='Y (pixels)',tit=ialsbase+'.als',charsize=1.2
        ps_close
        ps2gif,'html/'+ialsbase+'_sky.ps',rot=-90
        FILE_DELETE,'html/'+ialsbase+'_sky.ps',/allow,/quiet
      endif
      info[k].skyfile = skyfile1

      ; Make histogram plot
      ;--------------------
      histfile1 = ialsbase+'_hist.gif'
      test = FILE_TEST('html/'+histfile1)
      if test eq 0 or keyword_set(redo) then begin
        loadct,39,/silent
        ps_open,'html/'+ialsbase+'_hist',thick=5,/color
        binsize = 3.5*stddev(als.mag)/n_elements(als.mag)^(0.33)  ; Scott's choice binsize formula
        binsize = binsize/3.0   ; make a little smaller
        binsize = round(long(binsize*20.))/20.
        binsize = 0.05 > binsize < 1.0 
        plothist,als.mag,bin=binsize,xr=minmax(als.mag),xtit='Mag',ytit='N',tit=ialsbase+'.als',$
                 charsize=1.2
        plothist,als[gdstars].mag,bin=binsize,co=250,/over
        al_legend,['ALL','STARS'],textcolor=[0,250]
        ps_close
        ps2gif,'html/'+ialsbase+'_hist.ps',rot=-90
        FILE_DELETE,'html/'+ialsbase+'_hist.ps',/allow,/quiet
        loadct,0,/silent
      endif
      info[k].histfile = histfile1

      ; Make error vs. mag plot
      ;------------------------
      errmagfile1 = ialsbase+'_errmag.gif'
      test = FILE_TEST('html/'+errmagfile1)
      if test eq 0 or keyword_set(redo) then begin
        ps_open,'html/'+ialsbase+'_errmag',thick=5
        plot,als.mag,als.err,ps=8,sym=0.5,xr=minmax(als.mag),yr=minmax(als.err),xtit='Mag',$
                 ytit='Err',tit=ialsbase+'.als',charsize=1.2
        ps_close
        ps2gif,'html/'+ialsbase+'_errmag.ps',rot=-90
        FILE_DELETE,'html/'+ialsbase+'_errmag.ps',/allow,/quiet
      endif
      info[k].errmagfile = errmagfile1

      ; Make chi vs. sharp plot
      ;------------------------
      chisharpfile1 = ialsbase+'_chisharp.gif'
      test = FILE_TEST('html/'+chisharpfile1)
      if test eq 0 or keyword_set(redo) then begin
        ps_open,'html/'+ialsbase+'_chisharp',thick=5
        plot,als.chi,als.sharp,ps=8,sym=0.5,xr=minmax(als.chi),yr=minmax(als.sharp),xtit='Chi',$
                 ytit='Sharp',tit=ialsbase+'.als',charsize=1.2
        ps_close
        ps2gif,'html/'+ialsbase+'_chisharp.ps',rot=-90
        FILE_DELETE,'html/'+ialsbase+'_chisharp.ps',/allow,/quiet
      endif
      info[k].chisharpfile = chisharpfile1

      ; Make chi vs. mag plot
      ;----------------------
      chimagfile1 = ialsbase+'_chimag.gif'
      test = FILE_TEST('html/'+chimagfile1)
      if test eq 0 or keyword_set(redo) then begin
        ps_open,'html/'+ialsbase+'_chimag',thick=5
        plot,als.mag,als.chi,ps=8,sym=0.5,xr=minmax(als.mag),yr=[0,10],xtit='Mag',$
                 ytit='Chi',tit=ialsbase+'.als',charsize=1.2
        ps_close
        ps2gif,'html/'+ialsbase+'_chimag.ps',rot=-90
        FILE_DELETE,'html/'+ialsbase+'_chimag.ps',/allow,/quiet
      endif
      info[k].chimagfile = chimagfile1

      ; PSF CHI values
      ;---------------
      psflogfile = ialsbase+'.psf.log'
      info[k].psfchi = '???'
      if FILE_TEST(psflogfile) eq 1 then begin
        READLINE,psflogfile,psflines
        gpsf = where(stregex(psflines,'^>> ',/boolean) eq 1,ngpsf)
        if ngpsf gt 0 then begin
          chilines = psflines[gpsf]
          arr = strsplitter(chilines,' ',/extract)
          chi = float(arr[1,*])
          bestchi = min(chi)
          info[k].psfchi = strtrim(string(bestchi,format='(F10.4)'),2)
        endif
      endif

      ;stop

    Endfor ; loop through the separate ALS files

    ; Add to giant ALLINFO structure
    PUSH,allinfo,info



    ;----------------------------------
    ; FINAL RESULTS for this AMPLIFIER
    ;----------------------------------
    ; If namps=1 then this is the field final results

    ; Get field average RA, DEC, UT-Time, Date, Airmass and FWHM
    fieldampinfo[j].ra = mean(info.ra)
    fieldampinfo[j].dec = mean(info.dec)
    fieldampinfo[j].uttime = info[0].uttime
    fieldampinfo[j].date = info[0].date
    fieldampinfo[j].airmass = mean(info.airmass)
    fieldampinfo[j].fwhm = mean(info.fwhm)

    ; Get FINAL structure
    ; OR use the amplifier .AST file, must use undereddened photometry
    if (namps gt 1) then begin
      gdampstars = where(stregex(final.id,'^'+ifield+'_'+strtrim(iamp,2)+'\.',/boolean) eq 1,ngdampstars)
      str = final[gdampstars]
    endif else begin
      str = final
    endelse
    nstr = n_elements(str)
    fieldampinfo[j].nsources = nstr

    ; TYPE
    fieldampinfo[j].type = 'ALLSTAR'
    if FILE_TEST(ampbase+'.makemag') eq 1 and TAG_EXIST(str,'PROB') eq 1 then $
       fieldampinfo[j].type='ALLFRAME/SExtractor'
    if FILE_TEST(ampbase+'.makemag') eq 1 and TAG_EXIST(str,'PROB') eq 0 then $
       fieldampinfo[j].type='ALLFRAME/DAOPHOT'
    ;if TAG_EXIST(str,'PROB') eq 1 then fieldampinfo[j].type='ALLFRAME' else fieldampinfo[j].type='ALLSTAR'

   
    ; Only getting "good" sources
    gd = where(str.(magind) lt 50. and str.(col1mag1ind) lt 50. and $
               str.(col1mag2ind) lt 50. and abs(str.sharp) lt 1.0,ngd)
    str2 = str[gd]
    fieldampinfo[j].nstars = ngd

    ; Output filename
    if (namps gt 1) then namebase = ifield+'_'+strtrim(iamp,2) else $
      namebase = ifield

    ; Making CMD and 2CD
    ;-------------------
    cmdfile = namebase+'_cmd.gif'
    test = FILE_TEST('html/'+cmdfile)
    if (test eq 0 or keyword_set(redo)) and not keyword_set(nocmd) then begin
      ps_open,'html/'+namebase+'_cmd',thick=5
      device,xsize=9,ysize=10,/inches
      charsize=1.2
      mag = str2.(magind)
      err = str2.(magerrind)
      col1 = str2.(col1mag1ind) - str2.(col1mag2ind)
      col1range = median(col1)+[-3.0,3.0]*stddev(col1)
      plot,col1,mag,ps=8,sym=0.5,xr=col1range,yr=reverse(minmax(mag)),xs=1,ys=1,$
           xtit=' ',xtickformat='(A1)',ytit=magname,tit=namebase,$
           position=[0.10,0.35,0.98,0.97],charsize=charsize
      ; Getting faint limit for 2CD
      FITEXY,10.^(mag/2.5),err, A, B, X_SIG=mag*0.0,Y_SIG=sqrt(err)
      coef = [a,b]
      x = scale_vector(findgen(1000),min(mag),max(mag))
      modelerr = poly(10.^(x/2.5),coef)
      hi = first_el(where(modelerr gt 0.07))
      lim = x[hi]
      lim = round(lim*2.0)/2.0
      bright = where(mag le lim,nbright)
      if not keyword_set(no2cd) then begin
        col2 = str2.(col2mag1ind) - str2.(col2mag2ind)
        gdcol2 = where(str.(magind) le lim and str.(col2mag1ind) lt 50. and str.(col2mag2ind) lt 50.,ngdcol2)
        col2gd = str[gdcol2].(col2mag1ind)-str[gdcol2].(col2mag2ind)
        col2range = median(col2gd)+[-3.0,3.0]*stddev(col2gd)
        plot,[0],[0],/nodata,xr=col1range,yr=col2range,xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=col2name[0]+'-'+col2name[1],tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        oplot,col1[bright],col2[bright],ps=8,sym=0.5
        twocd_x = col1range[0]+range(col1range)*0.07
        twocd_y = col2range[0]+range(col2range)*0.85
        xyouts,twocd_x,twocd_y,magname+'<'+stringize(lim,ndec=1),charthick=4,charsize=1.2
      endif else begin
        plot,[0],[0],/nodata,xr=col1range,yr=[-1,1],xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=' ',tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        xyouts,mean(col1range),0.0,'CANNOT MAKE 2CD',charthick=4,charsize=2.0
      endelse
      ps_close
      ps2gif,'html/'+namebase+'_cmd.ps',rot=-90
      FILE_DELETE,'html/'+namebase+'_cmd.ps',/allow,/quiet
    endif
    fieldampinfo[j].cmdfile = cmdfile

    ; Image of combined FITS file
    ;----------------------------
    combfits = ampbase+'_comb.fits'
    if FILE_TEST(combfits) eq 1 then begin
      FITS_READ,combfits,combim,combhead
      sz2 = size(combim)
      nx2 = sz2[1]
      ny2 = sz2[2]

      combimagefile = namebase+'_comb_image.gif'
      test = FILE_TEST('html/'+combimagefile)
      if test eq 0 or keyword_set(redo) then begin
        fac = 10./float(nx2) < 10./float(ny2)
        xsize2 = nx2*fac
        ysize2 = ny2*fac + 1.5     ; plus 1 for the colorbar
        ps_open,'html/'+namebase+'_comb_image',thick=3,/portrait
        device,xsize=xsize2,ysize=ysize2,/inches
        displayc,combim,/z,xtit='X',ytit='Y',tit=combfits,charsize=1.1
        ps_close
        ps2gif,'html/'+namebase+'_comb_image.ps'
        FILE_DELETE,'html/'+namebase+'_comb_image.ps',/allow,/quiet
      endif
      fieldampinfo[j].imagefile = combimagefile
    endif

    ; CMD and 2CD with giant cut
    ;-------------------------------
    cmdgiantsfile = namebase+'_cmdgiants.gif'
    test = FILE_TEST('html/'+cmdgiantsfile)
    if (test eq 0 or keyword_set(redo)) and candogiants eq 1 and not keyword_set(nocmd) then begin
      loadct,39,/silent
      ps_open,'html/'+namebase+'_cmdgiants',thick=5,/color
      device,xsize=9,ysize=10,/inches
      charsize=1.2
      ; Select giants
      MTcol = str2.(Mmagind)-str2.(Tmagind)
      MDcol = str2.(Mmagind)-str2.(Dmagind)
      gdgiants = where( MTcol gt 0.8 and MTcol lt 2.5 and MDcol gt -0.18 and MDcol lt 0.30 and $
                        MDcol gt (-0.354750*MTcol+0.49) and MDcol gt (0.280605*MTcol-0.77),ngdgiants)
      mag = str2.(magind)
      err = str2.(magerrind)
      col1 = str2.(col1mag1ind) - str2.(col1mag2ind)
      col1range = median(col1)+[-3.0,3.0]*stddev(col1)
      magrange = reverse(minmax(mag))
      plot,[0],[0],/nodata,xr=col1range,yr=magrange,xs=1,ys=1,$
           xtit=' ',xtickformat='(A1)',ytit=magname,tit=namebase+' stars (giants in red)',$
           position=[0.10,0.35,0.98,0.97],charsize=charsize
      ps=8 & sym=0.5
      if n_elements(mag) gt 3000. then sym=0.2
      oplot,col1,mag,ps=ps,sym=sym
      if ngdgiants gt 0 then $
        oplot,col1[gdgiants],mag[gdgiants],ps=ps,sym=0.5,co=250
      cmd_x = col1range[0]+range(col1range)*0.07
      cmd_y = magrange[0]-range(magrange)*0.85
      xyouts,cmd_x,cmd_y,'Giants',charthick=4,charsize=1.2,color=250
      ; Getting faint limit for 2CD
      FITEXY,10.^(mag/2.5),err, A, B, X_SIG=mag*0.0,Y_SIG=sqrt(err)
      coef = [a,b]
      x = scale_vector(findgen(1000),min(mag),max(mag))
      modelerr = poly(10.^(x/2.5),coef)
      hi = first_el(where(modelerr gt 0.07))
      lim = x[hi]
      lim = round(lim*2.0)/2.0
      bright = where(mag le lim,nbright)
      gdbrightgiants = where( MTcol gt 0.8 and MTcol lt 2.5 and MDcol gt -0.18 and MDcol lt 0.30 and $
                              MDcol gt (-0.354750*MTcol+0.49) and MDcol gt (0.280605*MTcol-0.77) and $
                              mag le lim,ngdbrightgiants)
      if not keyword_set(no2cd) then begin
        col2 = str2.(col2mag1ind) - str2.(col2mag2ind)
        gdcol2 = where(str.(magind) le lim and str.(col2mag1ind) lt 50. and str.(col2mag2ind) lt 50.,ngdcol2)
        col2gd = str[gdcol2].(col2mag1ind)-str[gdcol2].(col2mag2ind)
        col2range = median(col2gd)+[-3.0,3.0]*stddev(col2gd)
        plot,[0],[0],/nodata,xr=col1range,yr=col2range,xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=col2name[0]+'-'+col2name[1],tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        oplot,col1[bright],col2[bright],ps=ps,sym=sym
        if ngdbrightgiants gt 0 then $
          oplot,col1[gdbrightgiants],col2[gdbrightgiants],ps=ps,sym=0.5,co=250
        twocd_x = col1range[0]+range(col1range)*0.07
        twocd_y = col2range[0]+range(col2range)*0.85
        xyouts,twocd_x,twocd_y,magname+'<'+stringize(lim,ndec=1),charthick=4,charsize=1.2
      endif else begin
        plot,[0],[0],/nodata,xr=col1range,yr=[-1,1],xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=' ',tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        xyouts,mean(col1range),0.0,'CANNOT MAKE 2CD',charthick=4,charsize=2.0
      endelse
      ps_close
      ps2gif,'html/'+namebase+'_cmdgiants.ps',rot=-90
      FILE_DELETE,'html/'+namebase+'_cmdgiants.ps',/allow,/quiet
    endif
    fieldampinfo[j].cmdgiantsfile = cmdgiantsfile


    ; Sky distribution with giant cut
    ;----------------------------------
    skygiantsfile = namebase+'_skygiants.gif'
    test = FILE_TEST('html/'+skygiantsfile)
    if (test eq 0 or keyword_set(redo)) and candogiants eq 1 then begin
      loadct,39,/silent
      ps_open,'html/'+namebase+'_skygiants',thick=5,/color
      ps=8 & sym=0.5
      if n_elements(str2) gt 3000 then sym=0.2
      ; Select giants
      MTcol = str2.(Mmagind)-str2.(Tmagind)
      MDcol = str2.(Mmagind)-str2.(Dmagind)
      gdgiants = where( MTcol gt 0.8 and MTcol lt 2.5 and MDcol gt -0.18 and MDcol lt 0.30 and $
                        MDcol gt (-0.354750*MTcol+0.49)and MDcol gt (0.280605*MTcol-0.77),ngdgiants)
      rarange = reverse(minmax(str.ra))
      decrange = minmax(str.dec)
      plot,str2.ra,str2.dec,ps=ps,sym=sym,xr=rarange,yr=decrange,xs=1,ys=1,$
               xtit='RA (deg)',ytit='DEC (deg)',tit=namebase+' stars (giants in red)',charsize=1.2
      if ngdgiants gt 0 then $
        oplot,str2[gdgiants].ra,str2[gdgiants].dec,ps=ps,sym=0.5,co=250
      sky_x = rarange[0]-range(rarange)*0.07
      sky_y = decrange[0]+range(decrange)*0.85
      xyouts,sky_x,sky_y,'Giants',charthick=4,charsize=1.2,color=250
      ps_close
      ps2gif,'html/'+namebase+'_skygiants.ps',rot=-90
      FILE_DELETE,'html/'+namebase+'_skygiants.ps',/allow,/quiet
    endif
    fieldampinfo[j].skygiantsfile = skygiantsfile

    ; Spatial distribution
    ;----------------------------
    skyfile = namebase+'_sky.gif'
    test = FILE_TEST('html/'+skyfile)
    if test eq 0 or keyword_set(redo) then begin
      ps_open,'html/'+namebase+'_sky',thick=5
      plot,str2.ra,str2.dec,ps=8,sym=0.5,xr=reverse(minmax(str.ra)),yr=minmax(str.dec),xs=1,ys=1,$
               xtit='RA (deg)',ytit='DEC (deg)',tit=namebase+' stars',charsize=1.2
      ps_close
      ps2gif,'html/'+namebase+'_sky.ps',rot=-90
      FILE_DELETE,'html/'+namebase+'_sky.ps',/allow,/quiet
    endif
    fieldampinfo[j].skyfile = skyfile

    ; Histogram of ALL and GOOD sources
    ;----------------------------------
    histfile = namebase+'_hist.gif'
    test = FILE_TEST('html/'+histfile)
    if test eq 0 or keyword_set(redo) then begin
      loadct,39,/silent
      ps_open,'html/'+namebase+'_hist',thick=5,/color
      mag = str.(magind)
      mag2 = str2.(magind)
      binsize = 3.5*stddev(mag2)/n_elements(mag2)^(0.33)  ; Scott's choice binsize formula
      binsize = binsize/3.0   ; make a little smaller
      binsize = round(long(binsize*20.))/20.
      binsize = 0.05 > binsize < 1.0 
      nobad = where(mag lt 50.,nbad)
      plothist,mag[nobad],bin=binsize,xr=minmax(mag[nobad]),xtit=colname,ytit='N',tit=namebase,$
               charsize=1.2
      plothist,mag2,bin=binsize,co=250,/over
      al_legend,['All Sources','Stars'],textcolor=[0,250]
      ps_close
      ps2gif,'html/'+namebase+'_hist.ps',rot=-90
      FILE_DELETE,'html/'+namebase+'_hist.ps',/allow,/quiet
    endif
    fieldampinfo[j].histfile = histfile

    ; Make error vs. mag plot
    ;------------------------
    errorfile = namebase+'_error.gif'
    test = FILE_TEST('html/'+errorfile)
    if test eq 0 or keyword_set(redo) then begin
      ps_open,'html/'+namebase+'_error',thick=5
      xsize3 = 10
      ysize3 = nfilters*1.7
      device,xsize=xsize3,ysize=ysize3,/inches
      charsize = 1.5
      xr = minmax(str2.(magind))
      yr = [0.0,0.48]
      y0 = 0.07
      dy = (0.98-y0)/nfilters
      for f=0,nfilters-1 do begin
        if f eq 0 then begin
          xtit=magname
          xform=''
          noerase=0
        endif else begin
          xtit=' '
          xform='(A1)'
          noerase=1
        endelse
        ymin = y0+f*dy
        ymax = ymin+dy
        plot,mag,str.(errindarr[f]),ps=8,sym=0.5,xr=xr,yr=yr,xs=1,ys=1,charsize=charsize,$
             xtit=xtit,ytit=tags[errindarr[f]],xtickformat=xform,position=[0.08,ymin,0.98,ymax],$
             noerase=noerase,yminor=2,xticklen=0.05
        xyouts,xr[0]+0.5,0.3,filters[f],charsize=charsize*1.3
      end
      ps_close
      ps2gif,'html/'+namebase+'_error.ps',rot=-90
      FILE_DELETE,'html/'+namebase+'_error.ps',/allow,/quiet
    endif
    fieldampinfo[j].errorfile = errorfile

    ; Make chi vs. sharp plot
    ;------------------------
    chisharpfile = namebase+'_chisharp.gif'
    test = FILE_TEST('html/'+chisharpfile)
    if test eq 0 or keyword_set(redo) then begin
      ps_open,'html/'+namebase+'_chisharp',thick=5
      plot,str.chi,str.sharp,ps=8,sym=0.5,xr=minmax(str.chi),yr=minmax(str.sharp),xtit='Chi',$
               ytit='Sharp',tit=namebase,charsize=1.2
      ps_close
      ps2gif,'html/'+namebase+'_chisharp.ps',rot=-90
      FILE_DELETE,'html/'+namebase+'_chisharp.ps',/allow,/quiet
    endif
    fieldampinfo[j].chisharpfile = chisharpfile

    ; Make chi vs. mag plot
    ;----------------------
    chimagfile = namebase+'_chimag.gif'
    test = FILE_TEST('html/'+chimagfile)
    if test eq 0 or keyword_set(redo) then begin
      ps_open,'html/'+namebase+'_chimag',thick=5
      mag = str.(magind)
      nobad = where(mag lt 50.,nbad)
      plot,mag[nobad],str[nobad].chi,ps=8,sym=0.5,xr=minmax(mag[nobad]),yr=[0,10],xtit=magname,$
               ytit='Chi',tit=namebase,charsize=1.2
      ps_close
      ps2gif,'html/'+namebase+'_chimag.ps',rot=-90
      FILE_DELETE,'html/'+namebase+'_chimag.ps',/allow,/quiet
    endif
    fieldampinfo[j].chimagfile = chimagfile



    ;-----------------
    ; Make HTML file
    ;-----------------
    undefine,hlines
    PUSH,hlines,'<html>'
    PUSH,hlines,'<head>'
    PUSH,hlines,'<title>'+namebase
    PUSH,hlines,'</title'
    PUSH,hlines,'</head>'
    PUSH,hlines,'<body>'
    PUSH,hlines,'<center><H1>PHOTRED Summary page for '+namebase+' ('+ampbase+')</H1></center>'
    PUSH,hlines,'<hr>'
    PUSH,hlines,'<p>'
    PUSH,hlines,'<p>'
    PUSH,hlines,'<center><H2>Final Results</H2>'
    PUSH,hlines,'<p>'
    PUSH,hlines,'<table border=1>'
    PUSH,hlines,'<tr><td><b>Field Name</b></td><td><center>'+ifield+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>Base Name</b></td><td><center>'+ampbase+'</center></td></tr>'
    if TAG_EXIST(str,'PROB') then type='ALLFRAME' else type='ALLSTAR'
    PUSH,hlines,'<tr><td><b>Type</b></td><td><center>'+type+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>RA</b></td><td><center>'+ten2sexig(fieldampinfo[j].ra/15.0)+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>DEC</b></td><td><center>'+ten2sexig(fieldampinfo[j].dec)+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>Total Sources</b></td><td><center>'+strtrim(nstr,2)+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>Total Stars</b></td><td><center>'+strtrim(ngd,2)+'</center></td></tr>'
    ; --CMD--
    PUSH,hlines,'<tr><td><b>CMD and 2CD</b></td>'
    test = FILE_TEST('html/'+fieldampinfo[j].cmdfile)
    if test eq 1 then begin
      PUSH,hlines,'<td><center><a href="'+fieldampinfo[j].cmdfile+'"><img src="'+$
                  fieldampinfo[j].cmdfile+'" height=400></a></center></td>'
    endif else begin
      PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
    endelse
    PUSH,hlines,'</tr>'
    ; --Combined Image--
    if FILE_TEST('html/'+combimagefile) eq 1 then begin
      result = query_gif('html/'+combimagefile,gifstr)
      dims = gifstr.dimensions
      height = 400 < dims[1]
      sheight = strtrim(long(height),2)
      PUSH,hlines,'<tr><td><b>Combined Image</b></td>'
      PUSH,hlines,'<td><center><a href="'+combimagefile+'"><img src="'+combimagefile+$
                  '" height='+sheight+'></a></center></td>'
      PUSH,hlines,'</tr>'
    endif
    ; --CMD and 2CD with giants--
    PUSH,hlines,'<tr><td><b>CMD and 2CD with giants</b></td>'
    test = FILE_TEST('html/'+fieldampinfo[j].cmdgiantsfile)
    if test eq 1 then begin
      PUSH,hlines,'<td><center><a href="'+fieldampinfo[j].cmdgiantsfile+'"><img src="'+$
                  fieldampinfo[j].cmdgiantsfile+'" height=300></a></center></td>'
    endif else begin
      PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
    endelse
    PUSH,hlines,'</tr>'
    ; --Sky plot with giants--
    PUSH,hlines,'<tr><td><b>Sky plot with giants</b></td>'
    test = FILE_TEST('html/'+fieldampinfo[j].skygiantsfile)
    if test eq 1 then begin
      PUSH,hlines,'<td><center><a href="'+fieldampinfo[j].skygiantsfile+'"><img src="'+$
                  fieldampinfo[j].skygiantsfile+'" height=300></a></center></td>'
    endif else begin
      PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
    endelse
    PUSH,hlines,'</tr>'
    ; --Sky plot--
    PUSH,hlines,'<tr><td><b>Sky Distribution</b></td>'
    PUSH,hlines,'<td><center><a href="'+skyfile+'"><img src="'+skyfile+'" height=300></a></center></td>'
    PUSH,hlines,'</tr>'
    ; --Histogram plot--
    PUSH,hlines,'<tr><td><b>Histogram</b></td>'
    PUSH,hlines,'<td><center><a href="'+histfile+'"><img src="'+histfile+'" height=300></a></center></td>'
    PUSH,hlines,'</tr>'
    ; --Error plots--
    PUSH,hlines,'<tr><td><b>Errors</b></td>'
    PUSH,hlines,'<td><center><a href="'+errorfile+'"><img src="'+errorfile+'" height=350></a></center></td>'
    PUSH,hlines,'</tr>'
    ; --Chi vs. Sharp plot--
    PUSH,hlines,'<tr><td><b>Chi vs. Sharp</b></td>'
    PUSH,hlines,'<td><center><a href="'+chisharpfile+'"><img src="'+chisharpfile+'" height=300></a></center></td>'
    PUSH,hlines,'</tr>'
    ; --Chi vs. Mag plot--
    PUSH,hlines,'<tr><td><b>Chi vs. Mag</b></td>'
    PUSH,hlines,'<td><center><a href="'+chimagfile+'"><img src="'+chimagfile+'" height=300></a></center></td>'
    PUSH,hlines,'</tr>'
    PUSH,hlines,'</table>'
    PUSH,hlines,''
    PUSH,hlines,''
    PUSH,hlines,''
    PUSH,hlines,'<center><H2>Individual File Information:</H2></center>'
    PUSH,hlines,'<table border=1>'
    PUSH,hlines,'<tr><td><b>Filename</b></td>'
    for l=0,nalsfiles-1 do begin
      PUSH,hlines,'<td><center><b><a href="'+info[l].filename+'.header">'+info[l].filename+'</a></b></center></td>'
    end
    PUSH,hlines,'</tr>'
    ; Image Sise
    PUSH,hlines,'<tr><td><b>Image Size</b></td>'
    for l=0,nalsfiles-1 do begin
      PUSH,hlines,'<td><center>'+strtrim(info[l].nx,2)+'x'+strtrim(info[l].ny,2)+'</center></td>'
    end
    PUSH,hlines,'</tr>'
    ; Filter
    PUSH,hlines,'<tr><td><b>Filter</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+info[l].filter+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Exptime
    PUSH,hlines,'<tr><td><b>Exptime</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+strtrim(info[l].exptime,2)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; UT-TIME
    PUSH,hlines,'<tr><td><b>UT-Time</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+info[l].uttime+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Date
    PUSH,hlines,'<tr><td><b>Date</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+info[l].date+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Airmass
    PUSH,hlines,'<tr><td><b>Airmass</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+strtrim(info[l].airmass,2)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Trans
    PUSH,hlines,'<tr><td><b>XY-Trans</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+info[l].trans+'</center></td>'
    PUSH,hlines,'</tr>'
    ; RA
    PUSH,hlines,'<tr><td><b>RA</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+ten2sexig(info[l].ra/15.0)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; DEC
    PUSH,hlines,'<tr><td><b>DEC</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+ten2sexig(info[l].dec)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; PSF Chi
    PUSH,hlines,'<tr><td><b>PSF Chi</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+info[l].psfchi+'</center></td>'
    PUSH,hlines,'</tr>'
    ; FWHM
    PUSH,hlines,'<tr><td><b>FWHM</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+strtrim(info[l].fwhm,2)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Saturation Level
    PUSH,hlines,'<tr><td><b>Saturation</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+strtrim(info[l].saturation,2)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Aperture Correction
    PUSH,hlines,'<tr><td><b>Aperture Correction</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+stringize(info[l].apcor,ndec=4)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Nsources
    PUSH,hlines,'<tr><td><b>Nsources</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+strtrim(info[l].nsources,2)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Nstars
    PUSH,hlines,'<tr><td><b>Nstars</b></td>'
    for l=0,nalsfiles-1 do PUSH,hlines,'<td><center>'+strtrim(info[l].nstars,2)+'</center></td>'
    PUSH,hlines,'</tr>'
    ; Image
    PUSH,hlines,'<tr><td><b>Image</b></td>'
    for l=0,nalsfiles-1 do begin
      PUSH,hlines,'<td><center><a href="'+info[l].imagefile+'"><img src="'+info[l].imagefile+$
                  '" height=250></a></center></td>'
    end
    PUSH,hlines,'</tr>'
    ; Histogram
    PUSH,hlines,'<tr><td><b>ALS Histogram</b></td>'
    for l=0,nalsfiles-1 do begin
      PUSH,hlines,'<td><center><a href="'+info[l].histfile+'"><img src="'+info[l].histfile+$
                  '" height=250></a></center></td>'
    end
    PUSH,hlines,'</tr>'
    ; Err vs. Mag plot
    PUSH,hlines,'<tr><td><b>ALS Err vs. Mag</b></td>'
    for l=0,nalsfiles-1 do begin
      PUSH,hlines,'<td><center><a href="'+info[l].errmagfile+'"><img src="'+info[l].errmagfile+$
                  '" height=250></a></center></td>'
    end
    PUSH,hlines,'</tr>'
    ; Chi vs. Sharp plot
    PUSH,hlines,'<tr><td><b>ALS Chi vs. Sharp</b></td>'
    for l=0,nalsfiles-1 do begin
      PUSH,hlines,'<td><center><a href="'+info[l].chisharpfile+'"><img src="'+info[l].chisharpfile+$
                  '" height=250></a></center></td>'
    end
    PUSH,hlines,'</tr>'
    ; Chi vs. Mag plot
    PUSH,hlines,'<tr><td><b>ALS Chi vs. Mag</b></td>'
    for l=0,nalsfiles-1 do begin
      PUSH,hlines,'<td><center><a href="'+info[l].chimagfile+'"><img src="'+info[l].chimagfile+$
                  '" height=250></a></center></td>'
    end
    PUSH,hlines,'</tr>'
    PUSH,hlines,'</table>'
    PUSH,hlines,''
    PUSH,hlines,'</center>'
    PUSH,hlines,''
    PUSH,hlines,'</body>'
    PUSH,hlines,'</html>'

    htmlfile = 'html/'+namebase+'.html'
    WRITELINE,htmlfile,hlines
    FILE_CHMOD,htmlfile,'755'o
    printlog,logfile,'Writing ',htmlfile
    htmlfilearr[i] = namebase+'.html'
    ;htmlfilearr[j] = namebase+'.html'
    fieldampinfo[j].htmlfile = namebase+'.html'

    AMP_BOMB:


  Endfor ; amps loop



  ;----------------------------------
  ; FINAL FIELD RESULTS
  ;----------------------------------
  ; only for multi-amp systems
  if (namps gt 1) then begin

    ; Get field average RA, DEC, UT-Time, Date, Airmass and FWHM
    fieldinfo[i].uttime = fieldampinfo[0].uttime
    fieldinfo[i].date = fieldampinfo[0].date
    fieldinfo[i].airmass = mean(fieldampinfo.airmass)
    fieldinfo[i].fwhm = mean(fieldampinfo.fwhm)
    fieldinfo[i].type = fieldampinfo[0].type

    ; Making CMD and 2CD
    ;-------------------
    cmdfile = ifield+'_cmd.gif'
    test = FILE_TEST('html/'+cmdfile)
    if (test eq 0 or keyword_set(redo)) and not keyword_set(nocmd) then begin
      ps_open,'html/'+ifield+'_cmd',thick=5
      device,xsize=9,ysize=10,/inches
      charsize=1.2
      mag = final2.(magind)
      err = final2.(magerrind)
      col1 = final2.(col1mag1ind) - final2.(col1mag2ind)
      col1range = median(col1)+[-3.0,3.0]*stddev(col1)
      plot,[0],[0],/nodata,xr=col1range,yr=reverse(minmax(mag)),xs=1,ys=1,$
           xtit=' ',xtickformat='(A1)',ytit=magname,tit=ifield,$
           position=[0.10,0.35,0.98,0.97],charsize=charsize
      ps=8 & sym=0.5
      if n_elements(mag) gt 3000. then sym=0.2
      oplot,col1,mag,ps=ps,sym=sym
      ; Getting faint limit for 2CD
      FITEXY,10.^(mag/2.5),err, A, B, X_SIG=mag*0.0,Y_SIG=sqrt(err)
      coef = [a,b]
      x = scale_vector(findgen(1000),min(mag),max(mag))
      modelerr = poly(10.^(x/2.5),coef)
      hi = first_el(where(modelerr gt 0.07))
      lim = x[hi]
      lim = round(lim*2.0)/2.0
      bright = where(mag le lim,nbright)
      if not keyword_set(no2cd) then begin
        col2 = final2.(col2mag1ind) - final2.(col2mag2ind)
        gdcol2 = where(final.(magind) le lim and final.(col2mag1ind) lt 50. and final.(col2mag2ind) lt 50.,ngdcol2)
        col2gd = final[gdcol2].(col2mag1ind)-final[gdcol2].(col2mag2ind)
        col2range = median(col2gd)+[-3.0,3.0]*stddev(col2gd)
        plot,[0],[0],/nodata,xr=col1range,yr=col2range,xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=col2name[0]+'-'+col2name[1],tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        oplot,col1[bright],col2[bright],ps=ps,sym=sym
        twocd_x = col1range[0]+range(col1range)*0.07
        twocd_y = col2range[0]+range(col2range)*0.85
        xyouts,twocd_x,twocd_y,magname+'<'+stringize(lim,ndec=1),charthick=4,charsize=1.2
      endif else begin
        plot,[0],[0],/nodata,xr=col1range,yr=[-1,1],xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=' ',tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        xyouts,mean(col1range),0.0,'CANNOT MAKE 2CD',charthick=4,charsize=2.0
      endelse
      ps_close
      ps2gif,'html/'+ifield+'_cmd.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_cmd.ps',/allow,/quiet
    endif
    fieldinfo[i].cmdfile = cmdfile

    ; CMD and 2CD color-coded by amplifier
    ;--------------------------------------
    cmdcolampfile = ifield+'_cmdcolamp.gif'
    test = FILE_TEST('html/'+cmdcolampfile)
    if (test eq 0 or keyword_set(redo)) and not keyword_set(nocmd) then begin
      loadct,39,/silent
      ps_open,'html/'+ifield+'_cmdcolamp',thick=5,/color
      device,xsize=9,ysize=10,/inches
      charsize=1.2
      mag = final2.(magind)
      err = final2.(magerrind)
      col1 = final2.(col1mag1ind) - final2.(col1mag2ind)
      col1range = median(col1)+[-3.0,3.0]*stddev(col1)
      plot,[0],[0],/nodata,xr=col1range,yr=reverse(minmax(mag)),xs=1,ys=1,$
           xtit=' ',xtickformat='(A1)',ytit=magname,tit=ifield,$
           position=[0.10,0.35,0.98,0.91],charsize=charsize
      ps=8 & sym=0.5
      if n_elements(mag) gt 3000. then sym=0.2
      ampcol = SCALE(final2.ext,[1,namps],[60,240])
      plots,col1,mag,ps=ps,sym=sym,color=ampcol,noclip=0
      ; Getting faint limit for 2CD
      FITEXY,10.^(mag/2.5),err, A, B, X_SIG=mag*0.0,Y_SIG=sqrt(err)
      coef = [a,b]
      x = scale_vector(findgen(1000),min(mag),max(mag))
      modelerr = poly(10.^(x/2.5),coef)
      hi = first_el(where(modelerr gt 0.07))
      lim = x[hi]
      lim = round(lim*2.0)/2.0
      bright = where(mag le lim,nbright)
      if not keyword_set(no2cd) then begin
        col2 = final2.(col2mag1ind) - final2.(col2mag2ind)
        gdcol2 = where(final.(magind) le lim and final.(col2mag1ind) lt 50. and final.(col2mag2ind) lt 50.,ngdcol2)
        col2gd = final[gdcol2].(col2mag1ind)-final[gdcol2].(col2mag2ind)
        col2range = median(col2gd)+[-3.0,3.0]*stddev(col2gd)
        plot,[0],[0],/nodata,xr=col1range,yr=col2range,xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=col2name[0]+'-'+col2name[1],tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        ampcol = SCALE(final2[bright].ext,[1,namps],[60,240])
        plots,col1[bright],col2[bright],ps=ps,sym=sym,color=ampcol,noclip=0
        twocd_x = col1range[0]+range(col1range)*0.07
        twocd_y = col2range[0]+range(col2range)*0.85
        xyouts,twocd_x,twocd_y,magname+'<'+stringize(lim,ndec=1),charthick=4,charsize=1.2
      endif else begin
        plot,[0],[0],/nodata,xr=col1range,yr=[-1,1],xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=' ',tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        xyouts,mean(col1range),0.0,'CANNOT MAKE 2CD',charthick=4,charsize=2.0
      endelse
      colpos = [0.10,0.97,0.98,0.99]
      COLORBAR,minrange=1,maxrange=namps,position=colpos,bottom=60,ncolors=240-60+1,format='(I3)'
      ps_close
      ps2gif,'html/'+ifield+'_cmdcolamp.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_cmdcolamp.ps',/allow,/quiet
    endif
    fieldinfo[i].cmdcolampfile = cmdcolampfile

    ; CMD and 2CD with giant cut
    ;-------------------------------
    cmdgiantsfile = ifield+'_cmdgiants.gif'
    test = FILE_TEST('html/'+cmdgiantsfile)
    if (test eq 0 or keyword_set(redo)) and candogiants eq 1 and not keyword_set(nocmd) then begin
      loadct,39,/silent
      ps_open,'html/'+ifield+'_cmdgiants',thick=5,/color
      device,xsize=9,ysize=10,/inches
      charsize=1.2
      ; Select giants
      MTcol = final2.(Mmagind)-final2.(Tmagind)
      MDcol = final2.(Mmagind)-final2.(Dmagind)
      gdgiants = where( MTcol gt 0.8 and MTcol lt 2.5 and MDcol gt -0.18 and MDcol lt 0.30 and $
                        MDcol gt (-0.354750*MTcol+0.49) and MDcol gt (0.280605*MTcol-0.77),ngdgiants)
      mag = final2.(magind)
      err = final2.(magerrind)
      col1 = final2.(col1mag1ind) - final2.(col1mag2ind)
      col1range = median(col1)+[-3.0,3.0]*stddev(col1)
      magrange = reverse(minmax(mag))
      plot,[0],[0],/nodata,xr=col1range,yr=magrange,xs=1,ys=1,$
           xtit=' ',xtickformat='(A1)',ytit=magname,tit=ifield+' stars (giants in red)',$
           position=[0.10,0.35,0.98,0.97],charsize=charsize
      ps=8 & sym=0.5
      if n_elements(mag) gt 3000. then sym=0.2
      oplot,col1,mag,ps=ps,sym=sym
      if ngdgiants gt 0 then $
        oplot,col1[gdgiants],mag[gdgiants],ps=ps,sym=0.5,co=250
      cmd_x = col1range[0]+range(col1range)*0.07
      cmd_y = magrange[0]-range(magrange)*0.85
      xyouts,cmd_x,cmd_y,'Giants',charthick=4,charsize=1.2,color=250
      ; Getting faint limit for 2CD
      FITEXY,10.^(mag/2.5),err, A, B, X_SIG=mag*0.0,Y_SIG=sqrt(err)
      coef = [a,b]
      x = scale_vector(findgen(1000),min(mag),max(mag))
      modelerr = poly(10.^(x/2.5),coef)
      hi = first_el(where(modelerr gt 0.07))
      lim = x[hi]
      lim = round(lim*2.0)/2.0
      bright = where(mag le lim,nbright)
      gdbrightgiants = where( MTcol gt 0.8 and MTcol lt 2.5 and MDcol gt -0.18 and MDcol lt 0.30 and $
                              MDcol gt (-0.354750*MTcol+0.49) and MDcol gt (0.280605*MTcol-0.77) and $
                              mag le lim,ngdbrightgiants)
      if not keyword_set(no2cd) then begin
        col2 = final2.(col2mag1ind) - final2.(col2mag2ind)
        gdcol2 = where(final.(magind) le lim and final.(col2mag1ind) lt 50. and final.(col2mag2ind) lt 50.,ngdcol2)
        col2gd = final[gdcol2].(col2mag1ind)-final[gdcol2].(col2mag2ind)
        col2range = median(col2gd)+[-3.0,3.0]*stddev(col2gd)
        plot,[0],[0],/nodata,xr=col1range,yr=col2range,xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=col2name[0]+'-'+col2name[1],tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        oplot,col1[bright],col2[bright],ps=ps,sym=sym
        if ngdbrightgiants gt 0 then $
          oplot,col1[gdbrightgiants],col2[gdbrightgiants],ps=ps,sym=0.5,co=250
        twocd_x = col1range[0]+range(col1range)*0.07
        twocd_y = col2range[0]+range(col2range)*0.85
        xyouts,twocd_x,twocd_y,magname+'<'+stringize(lim,ndec=1),charthick=4,charsize=1.2
      endif else begin
        plot,[0],[0],/nodata,xr=col1range,yr=[-1,1],xs=1,ys=1,$
             xtit=col1name[0]+'-'+col1name[1],ytit=' ',tit=' ',/noerase,$
             position=[0.10,0.06,0.98,0.35],charsize=charsize
        xyouts,mean(col1range),0.0,'CANNOT MAKE 2CD',charthick=4,charsize=2.0
      endelse
      ps_close
      ps2gif,'html/'+ifield+'_cmdgiants.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_cmdgiants.ps',/allow,/quiet
    endif
    fieldinfo[i].cmdgiantsfile = cmdgiantsfile


    ; Sky distribution with giant cut
    ;----------------------------------
    skygiantsfile = ifield+'_skygiants.gif'
    test = FILE_TEST('html/'+skygiantsfile)
    if (test eq 0 or keyword_set(redo)) and candogiants eq 1 then begin
      loadct,39,/silent
      ps_open,'html/'+ifield+'_skygiants',thick=5,/color
      ps=8 & sym=0.5
      if n_elements(final2) gt 3000 then sym=0.2
      ; Select giants
      MTcol = final2.(Mmagind)-final2.(Tmagind)
      MDcol = final2.(Mmagind)-final2.(Dmagind)
      gdgiants = where( MTcol gt 0.8 and MTcol lt 2.5 and MDcol gt -0.18 and MDcol lt 0.30 and $
                        MDcol gt (-0.354750*MTcol+0.49)and MDcol gt (0.280605*MTcol-0.77),ngdgiants)
      rarange = reverse(minmax(final.ra))
      decrange = minmax(final.dec)
      plot,final2.ra,final2.dec,ps=ps,sym=sym,xr=rarange,yr=decrange,xs=1,ys=1,$
               xtit='RA (deg)',ytit='DEC (deg)',tit=ifield+' stars (giants in red)',charsize=1.2
      if ngdgiants gt 0 then $
        oplot,final2[gdgiants].ra,final2[gdgiants].dec,ps=ps,sym=0.5,co=250
      sky_x = rarange[0]-range(rarange)*0.07
      sky_y = decrange[0]+range(decrange)*0.85
      xyouts,sky_x,sky_y,'Giants',charthick=4,charsize=1.2,color=250
      ps_close
      ps2gif,'html/'+ifield+'_skygiants.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_skygiants.ps',/allow,/quiet
    endif
    fieldinfo[i].skygiantsfile = skygiantsfile

    ; Spatial distribution
    ;----------------------------
    skyfile = ifield+'_sky.gif'
    test = FILE_TEST('html/'+skyfile)
    if test eq 0 or keyword_set(redo) then begin
      ps_open,'html/'+ifield+'_sky',thick=5
      ps=8 & sym=0.5
      if n_elements(final2) gt 3000 then sym=0.2
      plot,final2.ra,final2.dec,ps=ps,sym=sym,xr=reverse(minmax(final.ra)),yr=minmax(final.dec),xs=1,ys=1,$
               xtit='RA (deg)',ytit='DEC (deg)',tit=ifield+' stars',charsize=1.2
      ps_close
      ps2gif,'html/'+ifield+'_sky.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_sky.ps',/allow,/quiet
    endif
    fieldinfo[i].skyfile = skyfile

    ; Histogram of ALL and GOOD sources
    ;----------------------------------
    histfile = ifield+'_hist.gif'
    test = FILE_TEST('html/'+histfile)
    if test eq 0 or keyword_set(redo) then begin
      loadct,39,/silent
      ps_open,'html/'+ifield+'_hist',thick=5,/color
      mag = final.(magind)
      mag2 = final2.(magind)
      binsize = 3.5*stddev(mag2)/n_elements(mag2)^(0.33)  ; Scott's choice binsize formula
      binsize = binsize/3.0   ; make a little smaller
      binsize = round(long(binsize*20.))/20.
      binsize = 0.05 > binsize < 1.0 
      nobad = where(mag lt 50.,nbad)
      plothist,mag[nobad],bin=0.1,xr=minmax(mag[nobad]),xtit=magname,ytit='N',tit=ifield,$
               charsize=1.2
      plothist,final2.(magind),bin=0.1,co=250,/over
      al_legend,['All Sources','Stars'],textcolor=[0,250]
      ps_close
      ps2gif,'html/'+ifield+'_hist.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_hist.ps',/allow,/quiet
    endif
    fieldinfo[i].histfile = histfile

    ; Make error vs. mag plot
    ;------------------------
    errorfile = ifield+'_error.gif'
    test = FILE_TEST('html/'+errorfile)
    if test eq 0 or keyword_set(redo) then begin
      ps_open,'html/'+ifield+'_error',thick=5
      xsize3 = 10
      ysize3 = nfilters*1.7
      device,xsize=xsize3,ysize=ysize3,/inches
      charsize = 1.5
      xr = minmax(final2.m)
      yr = [0.0,0.48]
      y0 = 0.07
      dy = (0.98-y0)/nfilters
      for f=0,nfilters-1 do begin
        if f eq 0 then begin
          xtit=magname
          xform=''
          noerase=0
        endif else begin
          xtit=' '
          xform='(A1)'
          noerase=1
        endelse
        ymin = y0+f*dy
        ymax = ymin+dy
        mag = final.(magind)
        plot,mag,final.(errindarr[f]),ps=8,sym=0.5,xr=xr,yr=yr,xs=1,ys=1,charsize=charsize,$
             xtit=xtit,ytit=tags[errindarr[f]],xtickformat=xform,position=[0.08,ymin,0.98,ymax],$
             noerase=noerase,yminor=2,xticklen=0.05
        xyouts,xr[0]+0.5,0.3,filters[f],charsize=charsize*1.3
      end
      ps_close
      ps2gif,'html/'+ifield+'_error.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_error.ps',/allow,/quiet
    endif
    fieldinfo[i].errorfile = errorfile

    ; Make chi vs. sharp plot
    ;------------------------
    chisharpfile = ifield+'_chisharp.gif'
    test = FILE_TEST('html/'+chisharpfile)
    if test eq 0 or keyword_set(redo) then begin
      ps_open,'html/'+ifield+'_chisharp',thick=5
      plot,final.chi,final.sharp,ps=8,sym=0.5,xr=minmax(final.chi),yr=minmax(final.sharp),xtit='Chi',$
               ytit='Sharp',tit=ifield,charsize=1.2
      ps_close
      ps2gif,'html/'+ifield+'_chisharp.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_chisharp.ps',/allow,/quiet
    endif
    fieldinfo[i].chisharpfile = chisharpfile

    ; Make chi vs. mag plot
    ;----------------------
    chimagfile = ifield+'_chimag.gif'
    test = FILE_TEST('html/'+chimagfile)
    if test eq 0 or keyword_set(redo) then begin
      ps_open,'html/'+ifield+'_chimag',thick=5
      mag = final.(magind)
      nobad = where(mag lt 50.,nbad)
      plot,mag[nobad],final[nobad].chi,ps=8,sym=0.5,xr=minmax(mag[nobad]),yr=[0,10],xtit=magname,$
               ytit='Chi',tit=ifield,charsize=1.2
      ps_close
      ps2gif,'html/'+ifield+'_chimag.ps',rot=-90
      FILE_DELETE,'html/'+ifield+'_chimag.ps',/allow,/quiet
    endif
    fieldinfo[i].chimagfile = chimagfile

    ;stop


    ;------------------
    ; Make HTML file
    ;------------------
    undefine,hlines
    PUSH,hlines,'<html>'
    PUSH,hlines,'<head>'
    PUSH,hlines,'<title>'+ifield
    PUSH,hlines,'</title'
    PUSH,hlines,'</head>'
    PUSH,hlines,'<body>'
    PUSH,hlines,'<center><H1>PHOTRED Summary page for '+ifield+' ('+base+')</H1></center>'
    PUSH,hlines,'<hr>'
    PUSH,hlines,'<p>'
    PUSH,hlines,'<p>'
    PUSH,hlines,'<center><H2>Final Results</H2>'
    PUSH,hlines,'<p>'
    ; --Amplifier HTML pages--
    PUSH,hlines,'Amplifier HTML Summary Pages<br>'
    ampline = ''
    for l=0,namps-1 do begin
      ampline = ampline+'<a href="'+fieldampinfo[l].htmlfile+'">['+strtrim(fieldampinfo[l].amp,2)+']</a>'
    end
    PUSH,hlines,ampline
    PUSH,hlines,'<p>'
    PUSH,hlines,'<table border=1>'
    PUSH,hlines,'<tr><td><b>Field Name</b></td><td><center>'+ifield+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>Base Name</b></td><td><center>'+base+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>Type</b></td><td><center>'+fieldinfo[i].type+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>RA</b></td><td><center>'+ten2sexig(fieldinfo[i].ra/15.0)+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>DEC</b></td><td><center>'+ten2sexig(fieldinfo[i].dec)+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>Total Sources</b></td><td><center>'+strtrim(fieldinfo[i].nsources,2)+'</center></td></tr>'
    PUSH,hlines,'<tr><td><b>Total Stars</b></td><td><center>'+strtrim(fieldinfo[i].nstars,2)+'</center></td></tr>'
    ; --CMD--
    PUSH,hlines,'<tr><td><b>CMD and 2CD</b></td>'
    test = FILE_TEST('html/'+fieldinfo[i].cmdfile)
    if (test eq 1) then begin
      PUSH,hlines,'<td><center><a href="'+fieldinfo[i].cmdfile+'"><img src="'+$
                  fieldinfo[i].cmdfile+'" height=400></a></center></td>'
    endif else begin
      PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
    endelse
    PUSH,hlines,'</tr>'
    ; --CMD color-coded by amplifier--
    PUSH,hlines,'<tr><td><b>CMD and 2CD<br>color-coded by AMP</b></td>'
    test = FILE_TEST('html/'+fieldinfo[i].cmdcolampfile)
    if (test eq 1) then begin
      PUSH,hlines,'<td><center><a href="'+fieldinfo[i].cmdcolampfile+'"><img src="'+$
                  fieldinfo[i].cmdcolampfile+'" height=400></a></center></td>'
    endif else begin
      PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
    endelse
    PUSH,hlines,'</tr>'
    ; --CMD with giants--
    PUSH,hlines,'<tr><td><b>CMD and 2CD<br>with giants</b></td>'
    test = FILE_TEST('html/'+fieldinfo[i].cmdgiantsfile)
    if test eq 1 then begin
      PUSH,hlines,'<td><center><a href="'+fieldinfo[i].cmdgiantsfile+'"><img src="'+$
                fieldinfo[i].cmdgiantsfile+'" height=400></a></center></td>'
    endif else begin
      PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
    endelse
    PUSH,hlines,'</tr>'
    ; --Sky with giants plot--
    PUSH,hlines,'<tr><td><b>Sky Distribution<br>with giants</b></td>'
    test = FILE_TEST('html/'+fieldinfo[i].skygiantsfile)
    if test eq 1 then begin
      PUSH,hlines,'<td><center><a href="'+fieldinfo[i].skygiantsfile+'"><img src="'+fieldinfo[i].skygiantsfile+$
                  '" height=300></a></center></td>'
    endif else begin
      PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
    endelse
    PUSH,hlines,'</tr>'
    ; --Sky plot--
    PUSH,hlines,'<tr><td><b>Sky Distribution</b></td>'
    PUSH,hlines,'<td><center><a href="'+fieldinfo[i].skyfile+'"><img src="'+fieldinfo[i].skyfile+$
                '" height=300></a></center></td>'
    PUSH,hlines,'</tr>'
    ; --Histogram plot--
    PUSH,hlines,'<tr><td><b>Histogram</b></td>'
    PUSH,hlines,'<td><center><a href="'+fieldinfo[i].histfile+'"><img src="'+fieldinfo[i].histfile+$
                '" height=300></a></center></td>'
    PUSH,hlines,'</tr>'
    ; --Error plots--
    PUSH,hlines,'<tr><td><b>Errors</b></td>'
    PUSH,hlines,'<td><center><a href="'+fieldinfo[i].errorfile+'"><img src="'+fieldinfo[i].errorfile+$
                '" height=350></a></center></td>'
    PUSH,hlines,'</tr>'
    ; --Chi vs. Sharp plot--
    PUSH,hlines,'<tr><td><b>Chi vs. Sharp</b></td>'
    PUSH,hlines,'<td><center><a href="'+fieldinfo[i].chisharpfile+'"><img src="'+$
                fieldinfo[i].chisharpfile+'" height=300></a></center></td>'
    PUSH,hlines,'</tr>'
    ; --Chi vs. Mag plot--
    PUSH,hlines,'<tr><td><b>Chi vs. Mag</b></td>'
    PUSH,hlines,'<td><center><a href="'+fieldinfo[i].chimagfile+'"><img src="'+$
                fieldinfo[i].chimagfile+'" height=300></a></center></td>'
    PUSH,hlines,'</tr>'
    PUSH,hlines,'</table>'
    PUSH,hlines,''
    PUSH,hlines,'</center>'
    PUSH,hlines,''
    PUSH,hlines,'</body>'
    PUSH,hlines,'</html>'

    htmlfile = 'html/'+ifield+'.html'
    WRITELINE,htmlfile,hlines
    FILE_CHMOD,htmlfile,'755'o
    printlog,logfile,'Writing ',htmlfile
    htmlfilearr[i] = ifield+'.html'
    fieldinfo[i].htmlfile = ifield+'.html'

  ; Field Summary, Namp=1
  Endif else begin
    ; Copy fieldampinfo to fieldinfo
    ; the gifs were created in the amp=0 stage
    fieldinfo[i].cmdfile = fieldampinfo[0].cmdfile
    fieldinfo[i].imagefile = fieldampinfo[0].imagefile
    fieldinfo[i].skyfile = fieldampinfo[0].skyfile
    fieldinfo[i].histfile = fieldampinfo[0].histfile
    fieldinfo[i].errorfile = fieldampinfo[0].errorfile
    fieldinfo[i].chisharpfile = fieldampinfo[0].chisharpfile
    fieldinfo[i].chimagfile = fieldampinfo[0].chimagfile
    fieldinfo[i].cmdgiantsfile = fieldampinfo[0].cmdgiantsfile
    fieldinfo[i].skygiantsfile = fieldampinfo[0].skygiantsfile
  Endelse

  ;stop

  FIELD_BOMB:

END



;---------------------------------------
; Make the "digital" logsheet of FRAMES
;---------------------------------------
; only for multi-amp imagers
if (namps gt 1) then begin
  undefine,hlines
  PUSH,hlines,'<html>'
  PUSH,hlines,'<head>'
  PUSH,hlines,'<title>Digital Logsheet of Frames'
  PUSH,hlines,'</title'
  PUSH,hlines,'</head>'
  PUSH,hlines,'<body>'
  CD,current=curdir
  PUSH,hlines,'<center><H1>PHOTRED Digital Logsheet of Frames for<br>'
  PUSH,hlines,curdir+'/</H1></center>'
  PUSH,hlines,'<p>'
  PUSH,hlines,'<center>'
  PUSH,hlines,'<table border=1 cellpadding=6>'
  PUSH,hlines,'<tr><td><center>Filename</center></td><td><center>Field</center></td><td><center>Object</center></td>'+$
              '<td><center>Filter</center></td><td><center>Exptime</center></td><td><center>Date</center></td>'+$
              '<td><center>UT-Time</center></td><td><center>Airmass</center></td>'+$
              '<td><center>AVG Sky Median</center></td><td><center>AVG Sky RMS</center></td>'+$
              '<td><center>AVG WCS</center></td>'+$
              '<td><center>WCS Type</center></td><td><center>AVG WCS RMS</center></td><td><center>RA</center></td>'+$
              '<td><center>DEC</center></td><td><center>AVG PSF Chi</center></td><td><center>AVG FWHM</center></td>'+$
              '<td><center>AVG Saturation</center></td><td><center>AVG Ap.Corr.</center></td>'+$
              '<td><center>Nsources</center></td><td><center>Nstars</center></td></tr>'

  ui = uniq(allinfo.base,sort(allinfo.base))
  bases = allinfo[ui].base
  nbases = n_elements(bases)


  for i=0,nbases-1 do begin

    gbase = where(allinfo.base eq bases[i],ngbase)
    baseinfo = allinfo[gbase]

    ; Filename, object, filter, exptime, date, ut-time, airmass, fwhm, saturation, skymed, skymode, wcs, wcstype
    ; wcsrms, nsources, nstars, plots

    ;dum = {filename:'',field:'',object:'',nx:0L,ny:0L,filter:'',exptime:0.0,uttime:'',date:'',airmass:0.0,skymed:0.0,$
    ;       skysig:0.0,wcs:0,wcstype:'',wcsrms:0.0,trans:'',nsources:0L,nstars:0L,psfchi:'',fwhm:0.0,saturation:0.0,$
    ;       imagefile:'',skyfile:'',histfile:'',errmagfile:'',chisharpfile:'',chimagfile:''}

    thisline = '<tr>'
    ; Filename
    thisline = thisline+'<td><center>'+bases[i]+'</center></td>'
    ; Field
    thisline = thisline+'<td><center>'+baseinfo[0].field+'</center></td>'
    ; Object
    thisline = thisline+'<td><center>'+baseinfo[0].object+'</center></td>'
    ; Filter
    thisline = thisline+'<td><center>'+baseinfo[0].filter+'</center></td>'
    ; Exptime
    thisline = thisline+'<td><center>'+stringize(baseinfo[0].exptime,ndec=2)+'</center></td>'
    ; Date
    thisline = thisline+'<td><center>'+baseinfo[0].date+'</center></td>'
    ; UT-Time
    thisline = thisline+'<td><center>'+baseinfo[0].uttime+'</center></td>'
    ; Airmass
    thisline = thisline+'<td><center>'+stringize(baseinfo[0].airmass,ndec=5)+'</center></td>'
    ;; Image size
    ;thisline = thisline+'<td><center>'+strtrim(baseinfo[0].nx,2)+'x'+strtrim(baseinfo[0].ny,2)+'</center></td>'
    ; Skymed
    skymed = median(baseinfo.skymed)
    thisline = thisline+'<td><center>'+stringize(skymed,ndec=2)+'</center></td>'
    ; Skysig
    skysig = median(baseinfo.skysig)
    thisline = thisline+'<td><center>'+stringize(skysig,ndec=2)+'</center></td>'
    ; WCS
    if baseinfo[0].wcs eq 0 then havewcs='NO' else havewcs='YES'
    thisline = thisline+'<td><center>'+havewcs+'</center></td>'
    ; WCS type
    thisline = thisline+'<td><center>'+baseinfo[0].wcstype+'</center></td>'
    ; WCS RMS
    wcsrms = median(baseinfo.wcsrms)
    swcsrms = stringize(wcsrms,ndec=3)
    if wcsrms eq 99.99 then swcsrms='---'
    thisline = thisline+'<td><center>'+swcsrms+'</center></td>'
    ; RA
    ra = mean(baseinfo.ra)
    thisline = thisline+'<td><center>'+ten2sexig(ra/15.0)+'</center></td>'
    ; DEC
    dec = mean(baseinfo.dec)
    thisline = thisline+'<td><center>'+ten2sexig(dec)+'</center></td>'
    ; PSF Chi
    psfchi = median(float(baseinfo.psfchi))
    thisline = thisline+'<td><center>'+stringize(psfchi,ndec=4)+'</center></td>'
    ; FWHM
    fwhm = median(baseinfo.fwhm)
    thisline = thisline+'<td><center>'+stringize(fwhm,ndec=2)+'</center></td>'
    ; Saturation
    saturation = median(baseinfo.saturation)
    thisline = thisline+'<td><center>'+stringize(saturation,ndec=1)+'</center></td>'
    ; Aperture Correction
    apcor = median(baseinfo.apcor)
    thisline = thisline+'<td><center>'+stringize(apcor,ndec=4)+'</center></td>'
    ; Nsources
    nsources = total(baseinfo.nsources)
    thisline = thisline+'<td><center>'+strtrim(nsources,2)+'</center></td>'
    ; Nstars
    nstars = total(baseinfo.nstars)
    thisline = thisline+'<td><center>'+strtrim(nstars,2)+'</center></td>'
    ;; Image plot
    ;thisline = thisline+'<td><center><a href="'+baseinfo[0].imagefile+'">Image</a></center></td>'
    ;; Sky plot
    ;thisline = thisline+'<td><center><a href="'+baseinfo[0].skyfile+'">Sky</a></center></td>'
    ;; Histogram plot
    ;thisline = thisline+'<td><center><a href="'+baseinfo[0].histfile+'">Hist</a></center></td>'
    ;; Err vs. Mag plot
    ;thisline = thisline+'<td><center><a href="'+baseinfo[0].errmagfile+'">Error</a></center></td>'
    ;; Chi vs. Sharp plot
    ;thisline = thisline+'<td><center><a href="'+baseinfo[0].chisharpfile+'">Chi/Sharp</a></center></td>'
    ;; Chi vs. Mag plot
    ;thisline = thisline+'<td><center><a href="'+baseinfo[0].chimagfile+'">Chi/Mag</a></center></td>'
    thisline = thisline+'</tr>'
    PUSH,hlines,thisline
  end

  PUSH,hlines,'</table>'
  PUSH,hlines,''
  PUSH,hlines,'</center>'
  PUSH,hlines,''
  PUSH,hlines,'</body>'
  PUSH,hlines,'</html>'

  htmlfile = 'html/logsheet_frames.html'
  WRITELINE,htmlfile,hlines
  FILE_CHMOD,htmlfile,'755'o
  printlog,logfile,'Writing ',htmlfile

endif ; namps>1




;-----------------------------
; Make the "digital" logsheet
;-----------------------------
undefine,hlines
PUSH,hlines,'<html>'
PUSH,hlines,'<head>'
PUSH,hlines,'<title>Digital Logsheet'
PUSH,hlines,'</title'
PUSH,hlines,'</head>'
PUSH,hlines,'<body>'
CD,current=curdir
PUSH,hlines,'<center><H1>PHOTRED Digital Logsheet for<br>'
PUSH,hlines,curdir+'/</H1></center>'
PUSH,hlines,'<p>'
PUSH,hlines,'<center>'
PUSH,hlines,'<table border=1 cellpadding=6>'
PUSH,hlines,'<tr><td><center>Filename</center></td><td><center>Field</center></td><td><center>Object</center></td>'+$
            '<td><center>Filter</center></td><td><center>Exptime</center></td><td><center>Date</center></td>'+$
            '<td><center>UT-Time</center></td><td><center>Airmass</center></td><td><center>Size</center></td>'+$
            '<td><center>Sky Median</center></td><td><center>Sky RMS</center></td><td><center>WCS</center></td>'+$
            '<td><center>WCS Type</center></td><td><center>WCS RMS</center></td><td><center>RA</center></td>'+$
            '<td><center>DEC</center></td><td><center>PSF Chi</center></td><td><center>FWHM</center></td>'+$
            '<td><center>Saturation</center></td><td><center>Ap.Corr.</center></td>'+$
            '<td><center>Nsources</center></td><td><center>Nstars</center></td><td><center>Image</center></td>'+$
            '<td><center>Sky</center></td><td><center>Hist</center></td><td><center>Error</center></td>'+$
            '<td><center>Chi/Sharp</center></td><td><center>Chi/Mag</center></td></tr>'

nallinfo = n_elements(allinfo)
bases = strsplitter(allinfo.filename,'-',/extract)
bases = reform(bases[1,*])
si = sort(bases)
allinfo = allinfo[si]
for i=0,nallinfo-1 do begin

  ; Filename, object, filter, exptime, date, ut-time, airmass, fwhm, saturation, skymed, skymode, wcs, wcstype
  ; wcsrms, nsources, nstars, plots

  ;dum = {filename:'',field:'',object:'',nx:0L,ny:0L,filter:'',exptime:0.0,uttime:'',date:'',airmass:0.0,skymed:0.0,$
  ;       skysig:0.0,wcs:0,wcstype:'',wcsrms:0.0,trans:'',nsources:0L,nstars:0L,psfchi:'',fwhm:0.0,saturation:0.0,$
  ;       imagefile:'',skyfile:'',histfile:'',errmagfile:'',chisharpfile:'',chimagfile:''}

  thisline = '<tr>'
  ; Filename
  thisline = thisline+'<td><center><a href="'+allinfo[i].filename+'.header">'+allinfo[i].filename+'</a></center></td>'
  ; Field
  thisline = thisline+'<td><center>'+allinfo[i].field+'</center></td>'
  ; Object
  thisline = thisline+'<td><center>'+allinfo[i].object+'</center></td>'
  ; Filter
  thisline = thisline+'<td><center>'+allinfo[i].filter+'</center></td>'
  ; Exptime
  thisline = thisline+'<td><center>'+stringize(allinfo[i].exptime,ndec=2)+'</center></td>'
  ; Date
  thisline = thisline+'<td><center>'+allinfo[i].date+'</center></td>'
  ; UT-Time
  thisline = thisline+'<td><center>'+allinfo[i].uttime+'</center></td>'
  ; Airmass
  thisline = thisline+'<td><center>'+stringize(allinfo[i].airmass,ndec=5)+'</center></td>'
  ; Image size
  thisline = thisline+'<td><center>'+strtrim(allinfo[i].nx,2)+'x'+strtrim(allinfo[i].ny,2)+'</center></td>'
  ; Skymed
  thisline = thisline+'<td><center>'+stringize(allinfo[i].skymed,ndec=2)+'</center></td>'
  ; Skysig
  thisline = thisline+'<td><center>'+stringize(allinfo[i].skysig,ndec=2)+'</center></td>'
  ; WCS
  if allinfo[i].wcs eq 0 then havewcs='NO' else havewcs='YES'
  thisline = thisline+'<td><center>'+havewcs+'</center></td>'
  ; WCS type
  thisline = thisline+'<td><center>'+allinfo[i].wcstype+'</center></td>'
  ; WCS RMS
  thisline = thisline+'<td><center>'+stringize(allinfo[i].wcsrms,ndec=3)+'</center></td>'
  ; RA
  thisline = thisline+'<td><center>'+ten2sexig(allinfo[i].ra/15.0)+'</center></td>'
  ; DEC
  thisline = thisline+'<td><center>'+ten2sexig(allinfo[i].dec)+'</center></td>'
  ; PSF Chi
  thisline = thisline+'<td><center>'+allinfo[i].psfchi+'</center></td>'
  ; FWHM
  thisline = thisline+'<td><center>'+stringize(allinfo[i].fwhm,ndec=2)+'</center></td>'
  ; Saturation
  thisline = thisline+'<td><center>'+stringize(allinfo[i].saturation,ndec=1)+'</center></td>'
  ; Aperture Correction
  thisline = thisline+'<td><center>'+stringize(allinfo[i].apcor,ndec=4)+'</center></td>'
  ; Nsources
  thisline = thisline+'<td><center>'+strtrim(allinfo[i].nsources,2)+'</center></td>'
  ; Nstars
  thisline = thisline+'<td><center>'+strtrim(allinfo[i].nstars,2)+'</center></td>'
  ; Image plot
  thisline = thisline+'<td><center><a href="'+allinfo[i].imagefile+'">Image</a></center></td>'
  ; Sky plot
  thisline = thisline+'<td><center><a href="'+allinfo[i].skyfile+'">Sky</a></center></td>'
  ; Histogram plot
  thisline = thisline+'<td><center><a href="'+allinfo[i].histfile+'">Hist</a></center></td>'
  ; Err vs. Mag plot
  thisline = thisline+'<td><center><a href="'+allinfo[i].errmagfile+'">Error</a></center></td>'
  ; Chi vs. Sharp plot
  thisline = thisline+'<td><center><a href="'+allinfo[i].chisharpfile+'">Chi/Sharp</a></center></td>'
  ; Chi vs. Mag plot
  thisline = thisline+'<td><center><a href="'+allinfo[i].chimagfile+'">Chi/Mag</a></center></td>'
  thisline = thisline+'</tr>'
  PUSH,hlines,thisline
end

PUSH,hlines,'</table>'
PUSH,hlines,''
PUSH,hlines,'</center>'
PUSH,hlines,''
PUSH,hlines,'</body>'
PUSH,hlines,'</html>'

htmlfile = 'html/logsheet.html'
WRITELINE,htmlfile,hlines
FILE_CHMOD,htmlfile,'755'o
printlog,logfile,'Writing ',htmlfile




;---------------------------------------
; Make the "digital" logsheet PLUS PLOTS
;---------------------------------------
undefine,hlines
PUSH,hlines,'<html>'
PUSH,hlines,'<head>'
PUSH,hlines,'<title>Digital Logsheet Plus PLOTS'
PUSH,hlines,'</title'
PUSH,hlines,'</head>'
PUSH,hlines,'<body>'
CD,current=curdir
PUSH,hlines,'<center><H1>PHOTRED Digital Logsheet Plus PLOTS for<br>'
PUSH,hlines,curdir+'/</H1></center>'
PUSH,hlines,'<p>'
PUSH,hlines,'<center>'
PUSH,hlines,'<table border=1>'

PUSH,hlines,'<tr><td><b>Filename</b></td>'
for j=0,nallinfo-1 do begin
  PUSH,hlines,'<td><center><b><a href="'+allinfo[j].filename+'.header">'+allinfo[j].filename+'</a></b></center></td>'
end
PUSH,hlines,'</tr>'
; Object
PUSH,hlines,'<tr><td><b>Object</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+allinfo[j].object+'</center></td>'
PUSH,hlines,'</tr>'
; Filter
PUSH,hlines,'<tr><td><b>Filter</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+allinfo[j].filter+'</center></td>'
PUSH,hlines,'</tr>'
; Exptime
PUSH,hlines,'<tr><td><b>Exptime</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+strtrim(allinfo[j].exptime,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Date
PUSH,hlines,'<tr><td><b>Date</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+allinfo[j].date+'</center></td>'
PUSH,hlines,'</tr>'
; UT-TIME
PUSH,hlines,'<tr><td><b>UT-Time</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+allinfo[j].uttime+'</center></td>'
PUSH,hlines,'</tr>'
; Airmass
PUSH,hlines,'<tr><td><b>Airmass</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+strtrim(allinfo[j].airmass,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Image size
PUSH,hlines,'<tr><td><b>Size</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+strtrim(allinfo[j].nx,2)+'x'+$
                           strtrim(allinfo[j].ny,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Skymed
PUSH,hlines,'<tr><td><b>Sky-Median</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+stringize(allinfo[j].skymed,ndec=2)+'</center></td>'
PUSH,hlines,'</tr>'
; Skysig
PUSH,hlines,'<tr><td><b>Sky-RMS</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+stringize(allinfo[j].skysig,ndec=2)+'</center></td>'
PUSH,hlines,'</tr>'
; WCS
PUSH,hlines,'<tr><td><b>WCS</b></td>'
for j=0,nallinfo-1 do begin
  if allinfo[j].wcs eq 0 then havewcs='NO' else havewcs='YES'
  PUSH,hlines,'<td><center>'+havewcs+'</center></td>'
end
PUSH,hlines,'</tr>'
; WCS type
PUSH,hlines,'<tr><td><b>WCS-Type</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+allinfo[j].wcstype+'</center></td>'
PUSH,hlines,'</tr>'
; WCS RMS
PUSH,hlines,'<tr><td><b>WCS-RMS</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+stringize(allinfo[j].wcsrms,ndec=3)+'</center></td>'
PUSH,hlines,'</tr>'
; RA
PUSH,hlines,'<tr><td><b>RA</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+ten2sexig(allinfo[j].ra/15.0)+'</center></td>'
PUSH,hlines,'</tr>'
; DEC
PUSH,hlines,'<tr><td><b>DEC</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+ten2sexig(allinfo[j].dec)+'</center></td>'
PUSH,hlines,'</tr>'
; Trans
PUSH,hlines,'<tr><td><b>XY-Trans</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+allinfo[j].trans+'</center></td>'
PUSH,hlines,'</tr>'
; PSF Chi
PUSH,hlines,'<tr><td><b>PSF Chi</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+allinfo[j].psfchi+'</center></td>'
PUSH,hlines,'</tr>'
; FWHM
PUSH,hlines,'<tr><td><b>FWHM</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+strtrim(allinfo[j].fwhm,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Saturation Level
PUSH,hlines,'<tr><td><b>Saturation</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+strtrim(allinfo[j].saturation,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Aperture Correction
PUSH,hlines,'<tr><td><b>Ap.Corr.</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+stringize(allinfo[j].apcor,ndec=4)+'</center></td>'
PUSH,hlines,'</tr>'
; Nsources
PUSH,hlines,'<tr><td><b>Nsources</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+strtrim(allinfo[j].nsources,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Nstars
PUSH,hlines,'<tr><td><b>Nstars</b></td>'
for j=0,nallinfo-1 do PUSH,hlines,'<td><center>'+strtrim(allinfo[j].nstars,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Image
PUSH,hlines,'<tr><td><b>Image</b></td>'
for j=0,nallinfo-1 do begin
  PUSH,hlines,'<td><center><a href="'+allinfo[j].imagefile+'"><img src="'+allinfo[j].imagefile+$
              '" height=250></a></center></td>'
end
PUSH,hlines,'</tr>'
; Histogram
PUSH,hlines,'<tr><td><b>ALS Histogram</b></td>'
for j=0,nallinfo-1 do begin
  PUSH,hlines,'<td><center><a href="'+allinfo[j].histfile+'"><img src="'+allinfo[j].histfile+$
              '" height=250></a></center></td>'
end
PUSH,hlines,'</tr>'
; Err vs. Mag plot
PUSH,hlines,'<tr><td><b>ALS Err vs. Mag</b></td>'
for j=0,nallinfo-1 do begin
  PUSH,hlines,'<td><center><a href="'+allinfo[j].errmagfile+'"><img src="'+allinfo[j].errmagfile+$
              '" height=250></a></center></td>'
end
PUSH,hlines,'</tr>'
; Chi vs. Sharp plot
PUSH,hlines,'<tr><td><b>ALS Chi vs. Sharp</b></td>'
for j=0,nallinfo-1 do begin
  PUSH,hlines,'<td><center><a href="'+allinfo[j].chisharpfile+'"><img src="'+allinfo[j].chisharpfile+$
              '" height=250></a></center></td>'
end
PUSH,hlines,'</tr>'
; Chi vs. Mag plot
PUSH,hlines,'<tr><td><b>ALS Chi vs. Mag</b></td>'
for j=0,nallinfo-1 do begin
  PUSH,hlines,'<td><center><a href="'+allinfo[j].chimagfile+'"><img src="'+allinfo[j].chimagfile+$
              '" height=250></a></center></td>'
end
PUSH,hlines,'</tr>'
PUSH,hlines,'</table>'
PUSH,hlines,''
PUSH,hlines,'</center>'
PUSH,hlines,''
PUSH,hlines,'</body>'
PUSH,hlines,'</html>'

htmlfile = 'html/logsheetplots.html'
WRITELINE,htmlfile,hlines
FILE_CHMOD,htmlfile,'755'o
printlog,logfile,'Writing ',htmlfile






;-----------------------------
; Make the FIELDS Summary page
;-----------------------------
undefine,hlines
PUSH,hlines,'<html>'
PUSH,hlines,'<head>'
PUSH,hlines,'<title>Fields Summary Page'
PUSH,hlines,'</title'
PUSH,hlines,'</head>'
PUSH,hlines,'<body>'
CD,current=curdir
PUSH,hlines,'<center><H1>PHOTRED Fields Summary Page for<br>'
PUSH,hlines,curdir+'/</H1></center>'
PUSH,hlines,'<p>'
PUSH,hlines,'<center>'
PUSH,hlines,'<table border=1>'

PUSH,hlines,'<tr><td><b>Field</b></td>'
for j=0,ninputlines-1 do begin
  PUSH,hlines,'<td><center><b>'+fieldinfo[j].field+'</b></center></td>'
end
PUSH,hlines,'</tr>'
; Amplifiers
if (namps gt 1) then begin
  PUSH,hlines,'<tr><td><b>Amplifier HTML files</b></td>'
  for j=0,ninputlines-1 do begin
    line = '<td><center>'
    for k=0,namps-1 do begin
      htmlfile = fieldinfo[j].field+'_'+strtrim(k+1,2)+'.html'
      test = FILE_TEST('html/'+htmlfile)
      if test eq 1 then line=line+'<a href="'+htmlfile+'">['+strtrim(k+1,2)+']</a>' else $
         line=line+'['+strtrim(k+1,2)+']'
    end
    line = line+'</center></td>'
    PUSH,hlines,line
  end
  PUSH,hlines,'</tr>'
endif
; Base
PUSH,hlines,'<tr><td><b>Base</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+fieldinfo[j].base+'</center></td>'
PUSH,hlines,'</tr>'
; Date
PUSH,hlines,'<tr><td><b>Date</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+fieldinfo[j].date+'</center></td>'
PUSH,hlines,'</tr>'
; UT-TIME
PUSH,hlines,'<tr><td><b>UT-Time</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+fieldinfo[j].uttime+'</center></td>'
PUSH,hlines,'</tr>'
; Airmass
PUSH,hlines,'<tr><td><b>Airmass</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+strtrim(fieldinfo[j].airmass,2)+'</center></td>'
PUSH,hlines,'</tr>'
; FWHM
PUSH,hlines,'<tr><td><b>FWHM</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+strtrim(fieldinfo[j].fwhm,2)+'</center></td>'
PUSH,hlines,'</tr>'
; RA
PUSH,hlines,'<tr><td><b>RA</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+ten2sexig(fieldinfo[j].ra/15.0)+'</center></td>'
PUSH,hlines,'</tr>'
; DEC
PUSH,hlines,'<tr><td><b>DEC</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+ten2sexig(fieldinfo[j].dec)+'</center></td>'
PUSH,hlines,'</tr>'
; Nsources
PUSH,hlines,'<tr><td><b>Nsources</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+strtrim(fieldinfo[j].nsources,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Nstars
PUSH,hlines,'<tr><td><b>Nstars</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+strtrim(fieldinfo[j].nstars,2)+'</center></td>'
PUSH,hlines,'</tr>'
; Type
PUSH,hlines,'<tr><td><b>Type</b></td>'
for j=0,ninputlines-1 do PUSH,hlines,'<td><center>'+fieldinfo[j].type+'</center></td>'
PUSH,hlines,'</tr>'
; CMD
PUSH,hlines,'<tr><td><b>CMD and 2CD</b></td>'
for j=0,ninputlines-1 do begin
  PUSH,hlines,'<td><center><a href="'+fieldinfo[j].cmdfile+'"><img src="'+fieldinfo[j].cmdfile+$
              '" height=400></a></center></td>'
end
; Combined Image
; Only for single amp
if (namps eq 1) then begin
  PUSH,hlines,'<tr><td><b>Combined Image</b></td>'
  for j=0,ninputlines-1 do begin
    if fieldinfo[j].imagefile ne '' then $
      PUSH,hlines,'<td><center><a href="'+fieldinfo[j].imagefile+'"><img src="'+fieldinfo[j].imagefile+$
                  '" height=400></a></center></td>'
  end
  PUSH,hlines,'</tr>'
endif
if (namps gt 1) then begin
  ; CMD color-coded by amplifiers
  PUSH,hlines,'<tr><td><b>CMD and 2CD<br>color-coded by Amplifier</b></td>'
  for j=0,ninputlines-1 do begin
    PUSH,hlines,'<td><center><a href="'+fieldinfo[j].cmdcolampfile+'"><img src="'+fieldinfo[j].cmdcolampfile+$
                '" height=400></a></center></td>'
  end
endif
; CMD and 2CD with giants
PUSH,hlines,'<tr><td><b>CMD and 2CD<br>with giants</b></td>'
for j=0,ninputlines-1 do begin
  test = FILE_TEST('html/'+fieldinfo[j].cmdgiantsfile)
  if test eq 1 then begin
    PUSH,hlines,'<td><center><a href="'+fieldinfo[j].cmdgiantsfile+'"><img src="'+fieldinfo[j].cmdgiantsfile+$
                '" height=400></a></center></td>'
  endif else begin
    PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
  endelse
end
; Sky distribution with giants
PUSH,hlines,'<tr><td><b>Sky distribution<br> with giants</b></td>'
for j=0,ninputlines-1 do begin
  test = FILE_TEST('html/'+fieldinfo[j].skygiantsfile)
  if test eq 1 then begin
    PUSH,hlines,'<td><center><a href="'+fieldinfo[j].skygiantsfile+'"><img src="'+fieldinfo[j].skygiantsfile+$
                '" height=250></a></center></td>'
  endif else begin
    PUSH,hlines,'<td><center>Image NOT FOUND</center></td>'
  endelse
end
; Sky distribution
PUSH,hlines,'<tr><td><b>Sky distribution</b></td>'
for j=0,ninputlines-1 do begin
  PUSH,hlines,'<td><center><a href="'+fieldinfo[j].skyfile+'"><img src="'+fieldinfo[j].skyfile+$
              '" height=250></a></center></td>'
end
; Histogram
PUSH,hlines,'<tr><td><b>Histogram</b></td>'
for j=0,ninputlines-1 do begin
  PUSH,hlines,'<td><center><a href="'+fieldinfo[j].histfile+'"><img src="'+fieldinfo[j].histfile+$
              '" height=250></a></center></td>'
end
PUSH,hlines,'</tr>'
; Err vs. Mag plot
PUSH,hlines,'<tr><td><b>Errors</b></td>'
for j=0,ninputlines-1 do begin
  PUSH,hlines,'<td><center><a href="'+fieldinfo[j].errorfile+'"><img src="'+fieldinfo[j].errorfile+$
              '" height=350></a></center></td>'
end
PUSH,hlines,'</tr>'
; Chi vs. Sharp plot
PUSH,hlines,'<tr><td><b>Chi vs. Sharp</b></td>'
for j=0,ninputlines-1 do begin
  PUSH,hlines,'<td><center><a href="'+fieldinfo[j].chisharpfile+'"><img src="'+fieldinfo[j].chisharpfile+$
              '" height=250></a></center></td>'
end
PUSH,hlines,'</tr>'
; Chi vs. Mag plot
PUSH,hlines,'<tr><td><b>Chi vs. Mag</b></td>'
for j=0,ninputlines-1 do begin
  PUSH,hlines,'<td><center><a href="'+fieldinfo[j].chimagfile+'"><img src="'+fieldinfo[j].chimagfile+$
              '" height=250></a></center></td>'
end
PUSH,hlines,'</tr>'
PUSH,hlines,'</table>'
PUSH,hlines,''
PUSH,hlines,'</center>'
PUSH,hlines,''
PUSH,hlines,'</body>'
PUSH,hlines,'</html>'

htmlfile = 'html/fields.html'
WRITELINE,htmlfile,hlines
FILE_CHMOD,htmlfile,'755'o
printlog,logfile,'Writing ',htmlfile



;--------------------------------------------------
; NIGHT INDEX HTML file for this whole night
;--------------------------------------------------
undefine,hlines
PUSH,hlines,'<html>'
PUSH,hlines,'<head>'
PUSH,hlines,'<title>Index File'
PUSH,hlines,'</title'
PUSH,hlines,'</head>'
PUSH,hlines,'<body>'
PUSH,hlines,'<center><H1>PHOTRED Index Page for<br>'
PUSH,hlines,curdir+'/</H1></center>'
PUSH,hlines,'<hr>'
PUSH,hlines,'<p>'
PUSH,hlines,'<center>'
PUSH,hlines,'<table border=1>'
; Multiple amps
if namps gt 1 then begin
  PUSH,hlines,'<tr<td><b><center>Fields</center></b></td>'
  PUSH,hlines,'<td colspan='+strtrim(namps,2)+'><center>Amplifiers</center></td></tr>'
  for i=0,ninputlines-1 do begin
    ;PUSH,hlines,'<tr><td><center><a href="'+htmlfilearr[i]+'">'+fieldarr[i]+'</a></center></td></tr>'
    ifield = fieldinfo[i].field
    PUSH,hlines,'<tr><td><center><a href="'+fieldinfo[i].htmlfile+'">'+ifield+'</a></center></td>'
    for j=0,namps-1 do begin
      htmlfile = ifield+'_'+strtrim(j+1,2)+'.html'
      test = FILE_TEST('html/'+htmlfile)
      if test eq 1 then $
        PUSH,hlines,'<td><center><a href="'+htmlfile+'">'+strtrim(j+1,2)+'</a></center></td>'
      if test eq 0 then $
        PUSH,hlines,'<td><center>'+strtrim(j+1,2)+'</center></td>'
    end
  end
; Single AmP
endif else begin
  PUSH,hlines,'<tr<td><b><center>Fields</center></b></td></tr>'
  for i=0,ninputlines-1 do begin
    PUSH,hlines,'<tr><td><center><a href="'+fieldinfo[i].htmlfile+'">'+fieldinfo[i].field+'</a></center></td></tr>'
  end
endelse
PUSH,hlines,'</table>'
PUSH,hlines,'<p>'
PUSH,hlines,'<h3><a href="fields.html">Fields Summary</a></h3>'
if (namps gt 1) then begin
  PUSH,hlines,'<p>'
  PUSH,hlines,'<h3><a href="logsheet_frames.html">Digital Logsheet of Frames</a></h3>'
endif
PUSH,hlines,'<p>'
PUSH,hlines,'<h3><a href="logsheet.html">Digital Logsheet</a></h3>'
PUSH,hlines,'<p>'
PUSH,hlines,'<h3><a href="logsheetplots.html">Digital Logsheet plus PLOTS</a></h3>'
PUSH,hlines,'<p>'

; Add the Transformation Equations
PHOTRED_LOADSETUP,setup,count=count
transfile = READPAR(setup,'TRANS')
READLINE,transfile,translines
PUSH,hlines,'Photometric Transformation Equations ('+transfile[0]+'):'
PUSH,hlines,'<pre>'
PUSH,hlines,translines
PUSH,hlines,'</pre>'

; photred.setup
if FILE_TEST('photred.setup') eq 1 then begin
  FILE_COPY,'photred.setup','html/',/over,/allow
  PUSH,hlines,'<p>'
  PUSH,hlines,'<a href="photred.setup">photred.setup</a>'
endif else begin
  PUSH,hlines,'<p>'
  PUSH,hlines,'NO photred.setup file'  
endelse
; fields
if FILE_TEST('fields') eq 1 then begin
  FILE_COPY,'fields','html/',/over,/allow
  PUSH,hlines,'<p>'
  PUSH,hlines,'<a href="fields">fields</a>'
endif else begin
  PUSH,hlines,'<p>'
  PUSH,hlines,'NO fields file'
endelse
; filters
if FILE_TEST('filters') eq 1 then begin
  FILE_COPY,'filters','html/',/over,/allow
  PUSH,hlines,'<p>'
  PUSH,hlines,'<a href="filters">filters</a>'
endif else begin
  PUSH,hlines,'<p>'
  PUSH,hlines,'NO filters file'
endelse
; extinction
if FILE_TEST('extinction') eq 1 then begin
  FILE_COPY,'extinction','html/',/over,/allow
  PUSH,hlines,'<p>'
  PUSH,hlines,'<a href="extinction">extinction</a>'
endif else begin
  PUSH,hlines,'<p>'
  PUSH,hlines,'NO extinction file'
endelse

PUSH,hlines,''
PUSH,hlines,'</center>'
PUSH,hlines,''
PUSH,hlines,'</body>'
PUSH,hlines,'</html>'

htmlfile = 'html/index.html'
WRITELINE,htmlfile,hlines
FILE_CHMOD,htmlfile,'755'o
printlog,logfile,'Writing ',htmlfile


; Maybe put up the PHOTRED log files



;#####################
; SUMMARY of the Lists
;#####################
;PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
;                    failurelist=failurelist


printlog,logfile,'PHOTRED_HTML Finished  ',systime(0)

if keyword_set(stp) then stop

end
