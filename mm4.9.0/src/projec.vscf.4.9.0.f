C*********************************************************************
C*********************************************************************
        subroutine PROJEC(natom,nmode,xm,x0,omega,xl,xlat,rr,XP,
     1  TEMP,PAROT,ITAU)
c
*************************************************************************
*       natom - no. of atoms
*       nmode - no. of normal modes to return (note nmodes = natoms*3)
*       xm - masses in amu (eg. C = 12.0)
*       xlat - on input equilib geometry corresponding to order of masses in
*            Cartesian coordinates,in au.  Should be in the PA system
*            before the call to 'normals'
*       omega - normal mode frequencies (nmode in ascending order)
*              in atomic units
*       xl - mass-scaled normal mode eigenvectors corresponding to
*            ordering in omega.
C**     XP - PROJECTION MATRIX
*************************************************************************
        implicit real*8(a-h,o-z)
        LOGICAL ABINIT
        parameter (maxmod=180,np=60)
        dimension xm(natom),xlat(natom,3),omega(nmode),
C**NOTE RR IS USED FOR FORCE CONSTANT ARRAY AS WELL AS BOND LENGTHS
C**IT HAS ENOUGH STORAGE FOR RR(3*NATOM,3*NATOM)
     &  xl(natom,nmode,3),x0(natom,3),rr(3*natom,3*natom),
     &  XP(3*NATOM,3*NATOM),TEMP(3*NATOM,3*NATOM),PAROT(3,3)
        dimension coord0(maxmod)
        dimension fijk(6*maxmod,1),derivs(maxmod,maxmod),
     &          potens(-1:1,-1:1)
        dimension rnorms(maxmod,maxmod),coord(maxmod),
     &          freq(maxmod),cmass(np),flag(maxmod),glag(maxmod)
        dimension fv1(maxmod),fv2(maxmod)
      COMMON/FILASS/IOUT
      COMMON/ABINIT/ABINIT
      COMMON/FUNDAM/WAVENM
      COMMON/PRINT/IPRINT,JPRINT
c        bohr = 0.52917706d0
        bohr = 1.0d0
*************************************************************************
*  Xlat is a 60x3 array containing x,y,z coordinates for 60 atoms.  It  *
*  can be made larger, as necessary.  Atomic units are assumed.  For the*
*  moment, the masses are treated indiviually, however, soon this file  *
*  will be changed to include the masses of the atoms.                  *
*************************************************************************
        nmodes=3*natom
        n6=nmodes-nmode
        n7=nmodes-nmode+1
        ISIZE=3*NATOM
  
      IF(ABINIT)THEN
C**READ SECOND DERIVATIVE MATRIX HERE INTO DERIVS
        DO I=1,ISIZE
          READ(17,*)(DERIVS(J,I),J=1,ISIZE)
        END DO
C**USE ROTATION OF ORIGINAL GEOMETRY
        N1=0
        DO I1=1,NATOM
          DO J1=1,3
            N1=N1+1
            N2=0
            DO I2=1,NATOM
              DO J2=1,3
                N2=N2+1
                RNORMS(N2,N1)=0
                DO KK1=1,3
                  K1=(I1-1)*3+KK1
                  DO KK2=1,3
                    K2=(I2-1)*3+KK2
                    RNORMS(N2,N1)=RNORMS(N2,N1)+PAROT(KK2,J2)*
     1              PAROT(KK1,J1)*DERIVS(K2,K1)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        N1=0
        DO I1=1,NATOM
          DO J1=1,3
            N1=N1+1
            N2=0
            DO I2=1,NATOM
              DO J2=1,3
                N2=N2+1
                DERIVS(N2,N1)=RNORMS(N2,N1)
              END DO
            END DO
          END DO
        END DO
      ELSE
        do i=1,natom
          cmass(i) = sqrt(xm(i))
          do j=1,3
            xlat(i,j)=x0(i,j)/CMASS(I)
          end do
          ii = 3*i-2
*************************************************************************
*  Coord0 contains the mass-weighted equilibrium coordinates.   Xlat    *
*  will now be used for other purposes.                                 *
*************************************************************************
          coord0(ii) = xlat(i,1)*cmass(i)
          coord0(ii+1) = xlat(i,2)*cmass(i)
          coord0(ii+2) = xlat(i,3)*cmass(i)
        end do
*************************************************************************
*  Delta is the parameter used in numerical differentiation.            *
*************************************************************************
        delta = 0.003d0
        do i = 1,nmodes
          coord(i) = coord0(i)
        end do
        ntmd = 0
*************************************************************************
*  The xlat array is now used to hold coordinates in a form useful for  *
*  the potential.  The following example has coordinates, un-mass-weigh-*
*  ted, and expressed in bohrs.                                         *
*************************************************************************
        imm = 1
        do mm = 1,nmodes,3
          xlat(imm,1) = coord0(mm)*bohr/cmass(imm)
          xlat(imm,2) = coord0(mm+1)*bohr/cmass(imm)
          xlat(imm,3) = coord0(mm+2)*bohr/cmass(imm)
          imm = imm+1
        end do
*************************************************************************
*  The dpes_shell routine calculates the potential.  Input data are the *
*  coordinates.  Output, here assigned to potens(0,0), is the potential *
*  energy in au at equilibrium.                                         *
*************************************************************************
        call getpot(potens(0,0),natom,xlat,rr)
        do i = 1,nmodes
          do j = i,nmodes
            if (i.ne.j) then
              do ii = -1,1,2
                do jj = -1,1,2
*************************************************************************
*  For each mode, the coordinate is moved by a distance delta.  These   *
*  are put into xlat in potential-friendly units.  Then the potential is*
*  called, and the result is stored in a two dimensional array          *
*  corresponding to the east-west and north-south gradients.  This is   *
*  done when mode i <> mode j.                                          *
*************************************************************************
                  coord(i) = coord0(i) + ii*delta
                  coord(j) = coord0(j) + jj*delta
                  imm = 1
                  do mm = 1,nmodes,3
                    xlat(imm,1) = coord(mm)*bohr/cmass(imm)
                    xlat(imm,2) = coord(mm+1)*bohr/cmass(imm)
                    xlat(imm,3) = coord(mm+2)*bohr/cmass(imm)
                    imm = imm+1
                  end do
                  call getpot(potens(ii,jj),natom,xlat,rr)
                  coord(i) = coord0(i)
                  coord(j) = coord0(j)
                end do
              end do
*************************************************************************
*  The second derivatives are calculated using a four-point formula.    *
*************************************************************************
              derivs(i,j) = (potens(1,1) - potens(-1,1) -
     &          potens(1,-1) + potens(-1,-1))/(8.0d0*delta**2)
              derivs(j,i) = derivs(i,j)
            else
*************************************************************************
*  If i = j, then the derivative is calculated using a  3-point formula *
*  with the (0,0) displacement calculated at the top of the big loop.   *
*  (See above.)                                                         *
*************************************************************************
              do ii = -1,1,2
                coord(i) = coord0(i) + ii*delta
                imm = 1
                do mm = 1,nmodes,3
                  xlat(imm,1) = coord(mm)*bohr/cmass(imm)
                  xlat(imm,2) = coord(mm+1)*bohr/cmass(imm)
                  xlat(imm,3) = coord(mm+2)*bohr/cmass(imm)
                  imm = imm+1
                end do
                ntmd = ntmd + 1
                call getpot(potens(ii,0),natom,xlat,rr)
                coord(i) = coord0(i)
              end do
              derivs(i,i) = (potens(1,0) - 2*potens(0,0) +
     &          potens(-1,0))/(2.0d0*delta**2)
            endif
          enddo
        enddo
      END IF
  
C**NOW PROJECT OUT UNWANTED MODE
      DO I=1,ISIZE
        DO J=1,ISIZE
C**NOTE SCALE OF 10**6 IN GAUSSIAN
          IF(ABINIT)DERIVS(J,I)=DERIVS(J,I)*1.D-6
C**NOTE SCALE OF 10**6 IN GAUSSIAN
          RR(J,I)=DERIVS(J,I)
        END DO
      END DO
      IF(ITAU.GT.1.OR..NOT.ABINIT)THEN
        CALL DGEMM('N','N',ISIZE,ISIZE,ISIZE,-1.0D0,DERIVS,
     &  MAXMOD,XP,ISIZE,1.0D0,RR,ISIZE)
        CALL DGEMM('N','N',ISIZE,ISIZE,ISIZE,-1.0D0,XP,
     &  ISIZE,DERIVS,MAXMOD,1.0D0,RR,ISIZE)
        CALL DGEMM('N','N',ISIZE,ISIZE,ISIZE,1.0D0,DERIVS,
     &  MAXMOD,XP,ISIZE,0.0D0,TEMP,ISIZE)
        CALL DGEMM('N','N',ISIZE,ISIZE,ISIZE,1.0D0,XP,
     &  ISIZE,TEMP,ISIZE,1.0D0,RR,ISIZE)
      END IF
      DO I=1,ISIZE
        DO J=1,ISIZE
          DERIVS(J,I)=RR(J,I)
        END DO
      END DO
  
*************************************************************************
*  The derivs matrix is passed to rs.  The diagonal force constants are*
*  returned in freq.  The eigenvectors are returned as rnorms.  The     *
*  eigenvectors are not stored, but are recalculated from the deriv-    *
*  atives, which are stored.  If desired, file 5 can be used to store   *
*  the eigenvectors.                                                    *
*************************************************************************
      call rs(maxmod,nmodes,derivs,freq,1,rnorms,fv1,fv2,ierr)
c     write(6,10) (freq(k),k=1,nmodes)
10    format(1x,5d13.4)
      INEGS=0
      do i = 1,nmodes
*************************************************************************
*  The frequencies are determined from the diagonal force constants.    *
*************************************************************************
        flag(i)=1.0
        glag(i)=1.0
c  check for imaginary frequencies
        if (freq(i).lt.0.0) flag(i)=-1.0
C**Jelski formulae miss 2 in finite diffs.
        IF(ABINIT)THEN
          freq(i) = flag(i)*sqrt(dabs(freq(i)))
        ELSE
          freq(i) = flag(i)*sqrt(2*dabs(freq(i)))
        END IF
        if (freq(i)*WAVENM.lt.-10.0) glag(i)=-1.0
        IF(GLAG(I).LT.0)INEGS=INEGS+1

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)WRITE(IOUT,*)FLAG(I),GLAG(I),FREQ(I)*WAVENM
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      end do
      IF(ITAU.NE.1)THEN
c  shuffle into position
        N8=8
        DO I=1,NMODES
          IF(GLAG(I).LT.0.0)THEN
            X=FREQ(N8)
            FREQ(N8)=FREQ(I)
            FREQ(I)=X
            IC=0
            DO M=1,NATOM
              DO K=1,3
                IC=IC+1
                X=RNORMS(IC,N8)
                RNORMS(IC,N8)=RNORMS(IC,I)
                RNORMS(IC,I)=X
              END DO
            END DO
            N8=N8+1
          END IF
        END DO
        do i=1,nmode-1
          omega(i)=freq(i+n7)
        enddo
c  stuff rnorms into xl
        do n=1,nmode-1
          ic=0
          do m=1,natom
            do l=1,3
              ic=ic+1
              xl(m,n,l)=rnorms(ic,n+n7)
            enddo
          enddo
        enddo
      ELSE
        IF(ABINIT)THEN
          do i=1,nmode
            omega(i)=freq(i+n6)
          enddo
c  stuff rnorms into xl
          do n=1,nmode
            ic=0
            do m=1,natom
              do l=1,3
                ic=ic+1
                xl(m,n,l)=rnorms(ic,n+n6)
              enddo
            enddo
          enddo
        ELSE
          N6=N6+INEGS
          DO I=1,INEGS
            OMEGA(I)=-FREQ(I)
          END DO
          do i=1,nmode-INEGS
            omega(i+INEGS)=freq(i+n6)
          enddo
c  stuff rnorms into xl
          DO N=1,INEGS
            ic=0
            do m=1,natom
              do l=1,3
                ic=ic+1
                xl(m,N,l)=rnorms(ic,N)
              enddo
            enddo
          END DO
          do n=1,nmode-INEGS
            ic=0
            do m=1,natom
              do l=1,3
                ic=ic+1
                xl(m,n+INEGS,l)=rnorms(ic,n+n6)
              enddo
            enddo
          enddo
        END IF
      END IF
      return
      end
