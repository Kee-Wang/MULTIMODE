C*************************************************************
C*************************************************************
C**MATRIX INVERSION
C*************************************************************
C*************************************************************
      subroutine matinv (a,n,nd,ainv,iwork,ierr)
c
      implicit double precision (a-h,o-z)
c
      dimension a(nd,1),ainv(nd,1),iwork(3,1)
c
      equivalence (irow,jrow), (icolum,jcolum), (amax, t, swap)
c
c     initialization
c
      do 15 i = 1 , n
        do 10 j = 1 , n
          ainv (i,j) = a (i,j)
   10 continue
   15 continue
c
      ierr = 0
      do 20 j=1,n
   20 iwork(3,j)=0
   30 do 550 i=1,n
c
c     search for pivot element
c
   40 amax=0.0
   45 do 105 j=1,n
   50 if (iwork(3,j)-1) 60, 105, 60
   60 do 100 k=1,n
   70 if (iwork(3,k)-1) 80, 100, 740
   80    if ( abs(amax)- abs(a(j,k))) 85, 100, 100
   85 irow=j
   90 icolum=k
   95 amax=a(j,k)
  100 continue
  105 continue
       if ( abs(amax).eq.0.0)    go to 760
  110 iwork(3,icolum)=iwork(3,icolum)+1
c
c     interchange rows to put pivot element on diagonal
c
  130 if (irow-icolum) 140, 260, 140
  140    continue
  300 do 200 l=1,n
  160 swap=a(irow,l)
  170 a(irow,l)=a(icolum,l)
  200 a(icolum,l)=swap
  260 iwork(1,i)=irow
  270 iwork(2,i)=icolum
  310 pivi=a(icolum,icolum)
c
c     divide pivot row by pivot element
c
  330 a(icolum,icolum)=1.0
  340 do 350 l= 1,n
  350 a(icolum,l)=a(icolum,l)/pivi
c
c     reduce non-pivot rows
c
  380 do 550 l1=1,n
  390 if(l1-icolum) 400, 550, 400
  400 t=a(l1,icolum)
  420 a(l1,icolum)=0.0
  430 do 450 l=1,n
  450 a(l1,l)=a(l1,l)-a(icolum,l)*t
  550 continue
c
c     interchange columns
c
  600 do 710 i=1,n
  610 l=n+1-i
  620 if (iwork(1,l)-iwork(2,l)) 630, 710, 630
  630 jrow=iwork(1,l)
  640 jcolum=iwork(2,l)
  650 do 705 k=1,n
  660 swap=a(k,jrow)
  670 a(k,jrow)=a(k,jcolum)
  700 a(k,jcolum)=swap
  705 continue
  710 continue
  740 continue
      do 755 i = 1 , n
        do 750 j = 1 , n
          aux = ainv(i,j)
          ainv(i,j) = a(i,j)
          a(i,j) = aux
  750 continue
  755 continue
      return
  760 ierr = 129
       return
       end
