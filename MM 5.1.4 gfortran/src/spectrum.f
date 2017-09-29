      SUBROUTINE spectrum(FILEG,W)
      implicit real*8 (a-h,o-z)
      CHARACTER*80 FILEG
      real*4 s
      dimension s(20000)
      COMMON/MAXMIN/EMIN,EMAX
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7,INP8,INP9
500   FORMAT(///25X,' *** GAUSS SPECTRUM *** ',///)
501   FORMAT(/' OUTPUT GAUSS FILE [.GAU]: ',/,A80)
      WRITE(IOUT,500)
      WRITE(IOUT,501)FILEG
      CALL FLUSH(IOUT)
      
C     w=4 
      npoints=10001
C     lr=0.0
C     ur=5000.
      LR=EMIN
      UR=EMAX
      step=(ur-lr)/(npoints-1)
      
      do i=1,npoints
      s(i)=0.0
      enddo

      k=0
C     open (unit=99)
      REWIND INP8

 
C  5  read (unit=99,fmt='(37x,f8.0,42x,e11.0)',end=10) xnu,xint
   5  read (INP8,end=10) xnu,xint
      k=k+1
CC    write (6,*) k,xnu,xint

      
      do i=1,npoints
      z=(i-1)*step+lr
      s(i) = s(i)+xint*exp(-((z-xnu)/w)**2)
      enddo

      goto 5

C  10 close(unit=99)
   10 close(INP8)

      a=0.0
      do i=1,npoints
      if (s(i).gt.a) a=s(i)
      enddo

      do i=1,npoints
      s(i) = s(i)/a
c     s(i) = 1-s(i)
      enddo


      do i=1,npoints
      z=(i-1)*step+lr
C     write (12,'(f14.2,e16.4)') z,s(i)
      write (INP9,*) z,s(i)
      enddo

C     stop
      CLOSE(INP9)
      RETURN
      end
