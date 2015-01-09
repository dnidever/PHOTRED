pro goodpsf,lstfile,chifile,outfile

; This does the job of Tony's goodpsf program

minchi = 0.5
nmax = 1000

nlstfile = n_elements(lstfile)
nchifile = n_elements(chifile)
noutfile = n_elements(outfile)

; Not enough inputs
if nlstfile eq 0 or nchifile eq 0 or noutfile eq 0 then begin
  print,'Syntax - goodpsf,lstfile,chifile,outfile'
  return
end

;print,'LSTFILE = ',lstfile
;print,'CHIFILE = ',chifile
;print,'OUTFILE = ',outfile

;;###################################
;;# Reading in the filenames
;lstfile = ''
;chifile = ''
;outfile = ''
;
;print,'File with PSF candidates (.lst)?'
;read,lstfile
;print,'File with chi values of PSF candidates (.lst.chi)?'
;read,chifile
;print,'Output file ?'
;read,outfile

; Initializing the arrays
id = lonarr(nmax)
x = fltarr(nmax)
y = fltarr(nmax)
mag = fltarr(nmax)
err = fltarr(nmax)
sky = fltarr(nmax)

id1 = 0L
x1 = 0.0
y1 = 0.0
mag1 = 0.0
err1 = 0.0
sky1 = 0.0

;###################################
;# Reading in the list file
line11 = ''
line22 = ''
line = ''
count = 0

openr,unit,/get_lun,lstfile
readf,unit,line11             ; There are 3 extra lines at the beginning
readf,unit,line22
readf,unit,line

form = '(I7,F9.3,F9.3,F9.3,F9.3,F9.3)'

; Read in lines until end-of-file
while (~eof(unit)) do begin
  readf,unit,format=form,id1,x1,y1,mag1,err1,sky1

  id[count] = id1
  x[count] = x1
  y[count] = y1
  mag[count] = mag1
  err[count] = err1
  sky[count] = sky1

  count = count+1
end
close,unit
free_lun,unit

nstars = count

; Keeping only the ones read in
id = id[0:nstars-1]
x = x[0:nstars-1]
y = y[0:nstars-1]
mag = mag[0:nstars-1]
err = err[0:nstars-1]
sky = sky[0:nstars-1]

; Sorting the stars
si = sort(id)
id = id[si]
x = x[si]
y = y[si]
mag = mag[si]
err = err[si]
sky = sky[si]


;###################################
;# Reading in the chi file
num = lonarr(1000)
chi = fltarr(1000)
flag = strarr(1000)


line = ''
count = 0
openr,unit,/get_lun,chifile

; Read in lines until end-of-file
while (~eof(unit)) do begin
  readf,unit,line

  ; Loop through five columns
  for i=0,4 do begin
 
    line1 = strmid(line,i*17,17)
    num1 = strtrim(strmid(line1,0,7),2)
    chi1 = strtrim(strmid(line1,7,7),2)
    flag1 = strtrim(strmid(line1,14,3),2)

    ; We've got a star
    if num1 ne '' then begin
      num[count] = num1
      chi[count] = chi1
      flag[count] = flag1

      count = count + 1
    endif

  end ; five column loop

end  ; reading in lines
close,unit
free_lun,unit

nstars2 = count

; Only keeping the stars read in
num = num[0:nstars2-1]
chi = chi[0:nstars2-1]
flag = flag[0:nstars2-1]

; Sorting the stars
si = sort(num)
num = num[si]
chi = chi[si]
flag = flag[si]


; Making new arrays
id2 = lonarr(nstars2)
x2 = fltarr(nstars2)
y2 = fltarr(nstars2)
mag2 = fltarr(nstars2)
err2 = fltarr(nstars2)
sky2 = fltarr(nstars2)

;###################################
;# Matching the stars
for i=0,nstars2-1 do begin
  g = where(id eq num[i],ng)

  ; Transferring the data
  if (ng gt 0) then begin
    id2[i] = id[g[0]]
    x2[i] = x[g[0]]
    y2[i] = y[g[0]]
    mag2[i] = mag[g[0]]
    err2[i] = err[g[0]]
    sky2[i] = sky[g[0]]
  endif else begin
    id2[i] = -1
  endelse

end


;###################################
;# Only keeping the good ones
gd = where(flag ne '?' and flag ne '*' and chi lt minchi and id2 ge 0,ngd)


;####################################
;# Output the data
openw,unit,/get_lun,outfile
printf,unit,line11
printf,unit,line22
printf,unit,''
form = '(I7,F9.3,F9.3,F9.3,F9.3,F9.3)'
for i=0,ngd-1 do begin
  printf,unit,format=form,id2[gd[i]],x2[gd[i]],y2[gd[i]],mag2[gd[i]],err2[gd[i]],sky2[gd[i]]
end
close,unit
free_lun,unit

print,'You have ',strtrim(ngd,2),' (out of ',strtrim(nstars2,2),') PSF candidates left.

;stop

end
