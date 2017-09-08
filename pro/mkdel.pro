;+
;
; MKDEL
;
; IDL version of Tony's mkdel.f fortran program
;
;  Run this program AFTER running DAOGROW.
;
; The INF file should be in this format:
; obj1108_10a                     1  06 24  1.696    60.000
; obj1108_11a                     1  06 24  1.696    60.000
; obj1108_12a                     1  06 24  1.696    60.000
;
; Each file listed in the INF file must have associated
; .als (created by ALLSTAR) and .tot (created by DAOGROW) files.
;
; INPUTS:
;  inffile  Name of the .inf file.
;
; OUTPUTS:
;  This makes .del files for all the .als files
;  that are in the .inf file.
;
; USAGE:
;  IDL>mkdel,'n1.inf'
;
; By D.Nidever  December 2006 (copy of Tony's fortran code)
;-

pro mkdel,inffile

; No input
if n_elements(inffile) eq 0 then begin
  print,'Syntax - mkdel,inffile'
  return
endif

if file_test(inffile) eq 0 then begin
  print,'FILE ',inffile,' DOES NOT EXIST'
  return
end

; Read in the INF file
READCOL,inffile,infile,format='A',/silent

aalsfile = strtrim(infile,2)+'.als'
atotfile = strtrim(infile,2)+'.tot'
adelfile = strtrim(infile,2)+'.del'

nf = n_elements(infile)

print,strtrim(nf,2),' FILES'

; Loop through the files
for k=0,nf-1 do begin
   
  ; Read the .als file
  READCOL,aalsfile[k],id1,x1,y1,amag,amagerr,sky,it,chi,sharp,format='I,F,F,F,F,F,F,F,F',skipline=3,/silent
  nstars = n_elements(id1)

  ; Read the .tot file
  if file_test(atotfile[k]) eq 1 then begin
    READCOL,atotfile[k],id2,x2,y2,tmag,tmagerr,format='I,F,F,F,F',skipline=3,/silent

    ; Matching them
    MATCH,id1,id2,ind1,ind2,count=nind1
    id3 = id1[ind1]
    amag2 = amag[ind1]
    tmag2 = tmag[ind2]

    ; Calculating del
    del = amag2 - tmag2
    ndel = n_elements(del)

    ; Write the .del file
    openw,unit,/get_lun,adelfile[k]
    for i=0,ndel-1 do $
      printf,unit,format='(A1,I5,F9.4)','',id3[i],del[i]
    close,unit
    free_lun,unit
  endif else begin
    print,atotfile[k],' NOT FOUND'
  endelse
    
endfor


; THIS IS THE FORTRAN VERSION
;
;      PRINT *,'.inf file used in DAOGROW ?'
;      READ(*,'(A30)') inffile
;      OPEN(1,FILE=inffile)
;      DO i = 1, 1000
;         READ(1,'(1X,A11)',END=999) infile
;         aalsfile(i) = infile//'.als'
;         atotfile(i) = infile//'.tot'
;         adelfile(i) = infile//'.del'
;         nf = nf + 1
;         PRINT *,aalsfile(i)
;      ENDDO
;999   CLOSE(1)
;          
;
;      DO k = 1, nf
;         n = 0
;         totac = 0.
;         OPEN(10,FILE=adelfile(k))
;         OPEN(11,FILE=aalsfile(k))
;         READ(11,'(2(/))')
;         DO i = 1, 500
;            READ(11,50,ERR=99) id1,x1,y1,amag,amagerr,sky,it,chi,sharp
;            OPEN(12,FILE=atotfile(k))
;            READ(12,'(2(/))')
;            DO j = 1, 500
;               READ(12,60,ERR=9) id2,x2,y2,tmag,tmagerr
;               IF (id1.EQ.id2 .AND. tmag.LT.40.) THEN
;                  n = n + 1
;                  del = amag - tmag
;                  WRITE(10,'(1X,I5,F9.4)') id2,del
;                  totac = totac + del
;               ENDIF
;            ENDDO
;9           CLOSE(12)
;         ENDDO
;99       CLOSE(11)
;         CLOSE(10)
;      ENDDO
;
;50    FORMAT(I6,3F9.3,F9.4,F9.3,F9.0,2F9.3)
;60    FORMAT(1X,I5,2F9.3,2F9.4)
;

;stop 

end
