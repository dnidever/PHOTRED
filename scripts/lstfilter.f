C######################################################################
C#
C# Below, I include the FORTRAN77 code -- lstfilter.f, in case I might 
C# accidentally delete it.  Take out the preceding "#"s (commented), 
C# copy and paste to a file called "lstfilter.f", and compile the 
C# FORTRAN source using the following command.
C#
C# $ g77 lstfilter.f -o lstfilter
C#
C# Tested under gcc version 3.3.5.
C#
C######################################################################
C#
C#      PROGRAM LSTFILTER
C
C=====================================================================
C
C This program filters out bad sources from the .lst file
C
C=====================================================================
C
      INTEGER nmax
      PARAMETER(nmax=50000)
C
      INTEGER i,j,nstar,nnstar,nleft
      INTEGER id(nmax),cid(nmax)
      REAL x(nmax),y(nmax),mag(nmax),err(nmax),sky(nmax)
      REAL sharp(nmax),round(nmax),dum
      CHARACTER flag(nmax)*1,line*85
      CHARACTER*60 lstfile,coofile,outfile
      CHARACTER*69 line1,line2
      
      PRINT *,'File with PSF candiates (.lst)?'
      READ(*,'(A60)') lstfile
      PRINT *,'File with FIND output (.coo)?'
      READ(*,'(A60)') coofile
      PRINT *,'Output file (.lst1) ?'
      READ(*,'(A60)') outfile
C
C=====================================================================
C Read in .lst file
C=====================================================================
C
      OPEN(11,FILE=lstfile)
      READ(11,'(A69)') line1
      READ(11,'(A69)') line2
      READ(11,*)
      nstar = 0
      DO i = 1, nmax
         READ(11,*,END=9) id(i),x(i),y(i),mag(i),err(i),sky(i)
         nstar = nstar + 1
      ENDDO
9     CLOSE(11)
C           
C=====================================================================
C Read in .coo file
C=====================================================================
C
      OPEN(12,FILE=coofile)
      READ(12,*)
      READ(12,*)
      READ(12,*)
      nnstar = 0
      DO i = 1, nmax
         READ(12,*,END=99) cid(i),dum,dum,dum,sharp(i),round(i),dum
         nnstar = nnstar + 1
      ENDDO
99    CLOSE(12)
C           
C=====================================================================
C Filter out non-stellar sources
C=====================================================================
C
C This only keeps stars with 0.3 < sharp < 1.0
C Changes to  0.2 <= sharp <= 1.0 so it is consistent with DAOPHOT II
C   recommended cuts (DLN 08/05/08)

      nleft = nstar
      OPEN(21,FILE=outfile)
      WRITE(21,'(A69)') line1
      WRITE(21,'(A69)') line2
      WRITE(21,*)
      DO i = 1, nstar
         DO j = 1, nnstar
            IF (id(i).EQ.cid(j)) THEN
               IF (sharp(j).GE.0.2 .AND. sharp(j).LE.1.0) THEN
                  WRITE(21,70) id(i),x(i),y(i),mag(i),err(i),sky(i)
               ELSE
                  nleft = nleft-1
               ENDIF
               GOTO 111
            ENDIF
         ENDDO
111   ENDDO
70    FORMAT(I7,5F9.3)
C      
      WRITE(6,71) nleft,nstar
71    FORMAT(/'You have ',I3,' (out of ',I3,') stars left.'/)
C
      STOP
      END
