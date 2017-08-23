pro makemag,tfrfile,outfile,stp=stp,error=error

; This combines the ALLFRAME alf photometry output files

undefine,error

ntfrfile = n_elements(tfrfile)
noutfile = n_elements(outfile)
if ntfrfile eq 0 or noutfile eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - makemag,tfrfile,outfile'
  return
endif

test = file_test(tfrfile)
if test eq 0 then begin
  error = tfrfile+' NOT FOUND'
  print,error
  return
endif

; GETTING the number of files and filenames
;------------------------------------------
OPENR,unit,/get_lun,tfrfile
line=''
nfiles=0
undefine,linearr
WHILE (strmid(strtrim(line,2),0,3) ne '===') and ~EOF(unit) do begin
  line=''
  READF,unit,line
  PUSH,linearr,line
  nfiles++
end
CLOSE,unit
FREE_LUN,unit
nfiles = nfiles-1


; No files
if nfiles lt 1 then begin
  error = 'No files in '+tfrfile
  print,error
  return
endif

; Getting filenames
linearr2 = linearr[0:nfiles-1]
linearr2 = strtrim(linearr2,2)
arr = strsplitter(linearr2,' ',/extract)
files = reform(arr[0,*])

; Loading the entire TFR file
;----------------------------
nlines = FILE_LINES(tfrfile)
OPENR,unit,/get_lun,tfrfile
line=''
for i=0,nfiles do READF,unit,line
nrow = nlines-nfiles-1
ncol = 3+nfiles

; Only 9-18 per line, the rest wraps on multiple lines
; but ID/X/Y stay the same
id = lonarr(nrow)
x = fltarr(nrow)
y = fltarr(nrow)
num = lonarr(nrow,nfiles)-1
;arr = fltarr(ncol,nrow)
flag = 0
count = 0LL
istar = 0LL
while ~EOF(unit) do begin
  line = ''
  readf,unit,line
  arr = strsplit(line,' ',/extract)
  id1 = long(arr[0])
  x1 = float(arr[1])
  y1 = float(arr[2])
  num1 = long(arr[3:*])
  nnum1 = n_elements(num1)

  ; First line
  if count eq 0 then begin
    id[istar] = id1
    x[istar] = x1
    y[istar] = y1
    num[istar,0:nnum1-1] = num1
    numcount = nnum1

  ; Second or later lines
  endif else begin

    ; Same star, wrapped line
    if id1 eq id[istar] then begin
      num[istar,numcount:numcount+nnum1-1] = num1
      numcount += nnum1
    ; New star
    endif else begin
      istar++
      id[istar] = id1
      x[istar] = x1
      y[istar] = y1
      num[istar,0:nnum1-1] = num1
      numcount = nnum1
    endelse
  endelse

  count++
endwhile

;READF,unit,arr
CLOSE,unit
FREE_LUN,unit

; Trim excess lines
nstars = istar+1
id = id[0:nstars-1]
x = x[0:nstars-1]
y = y[0:nstars-1]
num = num[0:nstars-1,*]


magarr = fltarr(nstars,nfiles)+99.9999
magerrarr = fltarr(nstars,nfiles)+9.9999
skyarr = fltarr(nstars,nfiles)
iterarr = fltarr(nstars,nfiles)
chiarr = fltarr(nstars,nfiles)+!values.f_nan    ; NANs are ignored by MEDIAN
sharparr = fltarr(nstars,nfiles)+!values.f_nan
countarr = lonarr(nstars)        ; how many good mags for this stars

; Loop through the ALF files
;---------------------------
for i=0,nfiles-1 do begin

  if file_test(files[i]) eq 0 then begin
    error = files[i]+' NOT FOUND'
    print,error
    return
  endif

  ; Load the ALF file
  LOADALS,files[i],alf,alfhead

  ind = reform(num[*,i])
  ;ind = reform(arr[i+3,*])
  gd = where(ind gt 0,ngd)
  bd = where(ind eq 0,nbd)

  alfind = ind[gd]-1   ; idl indices

  if ngd gt 0 then begin
    magarr[gd,i] = alf[alfind].mag
    magerrarr[gd,i] = alf[alfind].err
    skyarr[gd,i] = alf[alfind].sky
    iterarr[gd,i] = alf[alfind].iter
    chiarr[gd,i] = alf[alfind].chi
    sharparr[gd,i] = alf[alfind].sharp
    countarr[gd]++
  endif

endfor

; Calculating the Median CHI and SHARP
if nfiles gt 1 then begin
  chi = MEDIAN([chiarr],dim=2)
  sharp = MEDIAN([sharparr],dim=2)
  ; The fortran makemag code used the mean
  chimean = total(chiarr,2,/nan)/countarr
  sharpmean = total(sharparr,2,/nan)/countarr
; only 1 file
endif else begin
  chi = chiarr
  sharp = sharparr
  chimean = chiarr
  sharpmean = sharparr
endelse

; Creating mag/magerr output array
magoutarr = fltarr(nstars,nfiles*2)
magoutarr[*,indgen(nfiles)*2] = magarr
magoutarr[*,indgen(nfiles)*2+1] = magerrarr


; Print the output
OPENW,unit,/get_lun,outfile

for i=0,nstars-1 do begin
 ;format(1x,I8,2f9.3,100f9.4)
 format = '(A1,I8,2F9.3,'+strtrim(nfiles*2+2,2)+'F9.4)'
 ;printf,unit,format=format,'',id[i],x[i],y[i],reform(magoutarr[i,*]),chi[i],sharp[i]
 printf,unit,format=format,'',id[i],x[i],y[i],reform(magoutarr[i,*]),chimean[i],sharpmean[i]
end


CLOSE,unit
FREE_LUN,unit


if keyword_set(stp) then stop


;chimedarr = chi#(fltarr(nfiles)+1.0)
;chidiff = chiarr - chimedarr
;chistd = mad(chiarr-chidum)
;chidiffstd = chidiff/chistd
;bd = where(abs(chidiffstd) gt 4.0,nbd)
;chiarr2 = chiarr
;if nbd gt 0 then chiarr2[bd] = !values.f_nan 
;chimed2 = MEDIAN(chiarr2,dim=2)

;stop

end

;      program pullout
;
;      implicit none
;      character*40 inpfile(30),tfrfile,outfile,sexinp
;      real jnk,chitot,routot,chi(30,60000),rou(30,60000)
;      real work(2,14),chis(30),rous(30),Xs(30),Ys(30)
;      real x(30,60000),y(30,60000),x1,y1
;      real iter(30,60000),sky(30,60000),M(30,60000),Merr(30,60000)
;      real T(30,60000),Terr(30,60000),chiavg,rouavg
;      integer i,j,k,id(30,60000),count,FMAX,LMAX(30),id1
;      integer num(30),num1(30),ids(30),ass,renumber,offset
;
;      write(6,*)'Enter the .tfr file to be read:'
;      read(5,*) tfrfile
;      open(9,file=tfrfile)
;
;      write(6,*)'Enter the number of files, id offset:'
;      read(5,*) FMAX,offset
;
;      write(6,*)'Enter the filename for output:'
;      read(5,*) outfile
;      open(8,file=outfile)
;
;      write(6,*)'Renumber (1) yes, (2) no?'
;      read(5,*) renumber
;      
;
;C read in the inpfile names  
;      i=11 
;      do j=1,FMAX
;        read(9,*)inpfile(j)
;        open(i,file=inpfile(j))
;        i=i+1
;      enddo
;
;C  REad in all the data
;
;      do i=1,FMAX
;
;C  Get rid of headers
;        do j=1,3
;           read(i+10,*)
;        enddo 
;
;        j=1
;        do while(j.ne.0)
;
;          read(i+10,*,END=300)id(i,j),X(i,j),Y(i,j),m(i,j),merr(i,j),
;     &    sky(i,j),iter(i,j),chi(i,j),rou(i,j)
;
;          j=j+1
;        enddo
;300     LMAX(i)=j-1
;        print*, 'Read ',inpfile(i),' which has ',LMAX(i),' lines.'
;      enddo
;C Get rid of stupid ---- line in .tfr file
;      read(9,*)
;
;      k=1  
;      do while(k.ne.0)
;
;C        read(9,*,END=400) id1,x,y,num1(1),num1(2),num1(3),num1(4)
;C     &  ,num1(5),num1(6),num1(7),num1(8),num1(9)
;
;         read(9,*,END=400) id1,X1,Y1,num(1:FMAX)
;
;C  Compute the mean
;        count=1  
;        chitot=0.
;        routot=0.
;        do i=1,FMAX
;           if(renumber.eq.2) then
;            ids(i)=id1+offset
;           else
;            ids(i)=k+offset
;           endif
;           if(num(i).eq.0) then
;             work(1,i)=99.9999
;             work(2,i)=9.9999
;           endif
;           if(num(i).ne.0) then
;             work(1,i)=m(i,num(i))
;             work(2,i)=merr(i,num(i))
;             Xs(count)=x1
;             Ys(count)=y1
;             chitot=chi(i,num(i))+chitot
;             routot=rou(i,num(i))+routot
;             count=count+1
;           endif
;C           print*, num1(i),i,work(1:2,i),work(1:2,1),count,
;C     &     xs(count),ys(count)
;        enddo
;C        print*, work(1:2,1:FMAX)
;C        pause
;        chiavg=chitot/dble(count-1)
;        rouavg=routot/dble(count-1)
;        write(8,48) ids(1),Xs(1),Ys(1),
;     &  work(1:2,1:FMAX),chiavg,rouavg
;45      format(1x,I8,4f9.4)
;48      format(1x,I8,2f9.3,100f9.4)
;        k=k+1
;       enddo
;400    print*, k-1 
;      stop
;      end
;
