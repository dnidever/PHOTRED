      program pullout

      implicit none
      character*40 inpfile(30),tfrfile,outfile,sexinp
      real jnk,chitot,routot,chi(30,60000),rou(30,60000)
      real work(2,14),chis(30),rous(30),Xs(30),Ys(30)
      real x(30,60000),y(30,60000),x1,y1
      real iter(30,60000),sky(30,60000),M(30,60000),Merr(30,60000)
      real T(30,60000),Terr(30,60000),chiavg,rouavg
      integer i,j,k,id(30,60000),count,FMAX,LMAX(30),id1
      integer num(30),num1(30),ids(30),ass,renumber,offset

      write(6,*)'Enter the .tfr file to be read:'
      read(5,*) tfrfile
      open(9,file=tfrfile)

      write(6,*)'Enter the number of files, id offset:'
      read(5,*) FMAX,offset

      write(6,*)'Enter the filename for output:'
      read(5,*) outfile
      open(8,file=outfile)

      write(6,*)'Renumber (1) yes, (2) no?'
      read(5,*) renumber
      

C read in the inpfile names  
      i=11 
      do j=1,FMAX
        read(9,*)inpfile(j)
        open(i,file=inpfile(j))
        i=i+1
      enddo

C  REad in all the data

      do i=1,FMAX

C  Get rid of headers
        do j=1,3
           read(i+10,*)
        enddo 

        j=1
        do while(j.ne.0)

          read(i+10,*,END=300)id(i,j),X(i,j),Y(i,j),m(i,j),merr(i,j),
     &    sky(i,j),iter(i,j),chi(i,j),rou(i,j)

          j=j+1
        enddo
300     LMAX(i)=j-1
        print*, 'Read ',inpfile(i),' which has ',LMAX(i),' lines.'
      enddo
C Get rid of stupid ---- line in .tfr file
      read(9,*)

      k=1  
      do while(k.ne.0)

C        read(9,*,END=400) id1,x,y,num1(1),num1(2),num1(3),num1(4)
C     &  ,num1(5),num1(6),num1(7),num1(8),num1(9)

         read(9,*,END=400) id1,X1,Y1,num(1:FMAX)

C  Compute the mean
        count=1  
        chitot=0.
        routot=0.
        do i=1,FMAX
           if(renumber.eq.2) then
            ids(i)=id1+offset
           else
            ids(i)=k+offset
           endif
           if(num(i).eq.0) then
             work(1,i)=99.9999
             work(2,i)=9.9999
           endif
           if(num(i).ne.0) then
             work(1,i)=m(i,num(i))
             work(2,i)=merr(i,num(i))
             Xs(count)=x1
             Ys(count)=y1
             chitot=chi(i,num(i))+chitot
             routot=rou(i,num(i))+routot
             count=count+1
           endif
C           print*, num1(i),i,work(1:2,i),work(1:2,1),count,
C     &     xs(count),ys(count)
        enddo
C        print*, work(1:2,1:FMAX)
C        pause
        chiavg=chitot/dble(count-1)
        rouavg=routot/dble(count-1)
        write(8,48) ids(1),Xs(1),Ys(1),
     &  work(1:2,1:FMAX),chiavg,rouavg
45      format(1x,I8,4f9.4)
48      format(1x,I8,2f9.3,100f9.4)
        k=k+1
       enddo
400    print*, k-1 
      stop
      end

