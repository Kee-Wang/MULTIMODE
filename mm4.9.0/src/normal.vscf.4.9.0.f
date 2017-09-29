C******************************************************************
C******************************************************************
C**NORMAL COORDINATE ANALYSIS
C******************************************************************
C******************************************************************
*  Written by                                                           *
*       Daniel Jelski                                                   *
*       July, 1995                                                      *
*  Modified by Joel Bowman, Sept. 1997
*  Interface version to be used with Carter VSCF code
*************************************************************************
c  needs to be linked to getpot(v,natom,xlat,rr) routine
c
C       CALL NORMAL(NATOM,NMODE,XM,X0,OMEGA,XL,XX,RR)

        subroutine normal(natom,nmode,xm,x0,omega,xl,xlat,rr)
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
*************************************************************************
        implicit real*8(a-h,o-z)
        parameter (maxmod=180,np=60)
        dimension xm(natom),xlat(natom,3),omega(nmode),
     &  xl(natom,nmode,3),x0(natom,3),rr(natom,natom)
        dimension coord0(maxmod)
        dimension fijk(6*maxmod,1),derivs(maxmod,maxmod),
     &          potens(-1:1,-1:1)
        dimension rnorms(maxmod,maxmod),coord(maxmod),
     &          freq(maxmod),cmass(np),flag(maxmod)
        dimension fv1(maxmod),fv2(maxmod)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/VMIN/VMIN
      COMMON/SADDLE/JNORM
      COMMON/FILASS/IOUT
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
          do i=1,natom
           do j=1,3
             xlat(i,j)=x0(i,j)
           end do
          cmass(i) = sqrt(xm(i))
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
C         delta = 0.005d0
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
        VMIN=POTENS(0,0)
        WRITE(IOUT,*)
        WRITE(IOUT,*)'V AT EQU ',POTENS(0,0)*WAVENM
        WRITE(IOUT,*)
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
     &            potens(1,-1) + potens(-1,-1))/(8.0d0*delta**2)
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
     &            potens(-1,0))/(2.0d0*delta**2)
              endif
            enddo
          enddo
*************************************************************************
*  The derivs matrix is passed to rs.  The diagonal force constants are*
*  returned in freq.  The eigenvectors are returned as rnorms.  The     *
*  eigenvectors are not stored, but are recalculated from the deriv-    *
*  atives, which are stored.  If desired, file 5 can be used to store   *
*  the eigenvectors.                                                    *
*************************************************************************
        call rs(maxmod,nmodes,derivs,freq,1,rnorms,fv1,fv2,ierr)
c       write(6,10) (freq(k),k=1,nmodes)
10      format(1x,5d13.4)
      WRITE(2,*)
      WRITE(2,444)
444   FORMAT(10X,'FLAG',22X,'FREQ',/)
      INEGS=0
      do i = 1,nmodes
*************************************************************************
*  The frequencies are determined from the diagonal force constants.    *
*************************************************************************
        flag(i)=1.0
c  check for imaginary frequencies (<-10cm-1)
        if (freq(i).lt.0.0) flag(i)=-1.0
        freq(i) = flag(i)*sqrt(2*dabs(freq(i)))
        IF(FREQ(I)*WAVENM.LT.-10)INEGS=INEGS+1
        WRITE(2,*)FLAG(I),FREQ(I)*WAVENM
      end do
C**IMAGINARY FREQUENCY WILL ALWAYS BE MODE 1 (BUT USE JNORM)
      IF(JNORM.GT.0)THEN
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
      return
      end
C*********************************************************************
C*********************************************************************

      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a REAL SYMMETRIC matrix.
C
C     On Input
C
C        NM  must be set to the row dimension of the two-dimensional
C        array parameters as declared in the calling program
C        dimension statement.
C
C        N  is the order of the matrix  A.
C
C        A  contains the real symmetric matrix.
C
C        MATZ  is an integer variable set equal to zero if
C        only eigenvalues are desired.  Otherwise it is set to
C        any non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        W  contains the eigenvalues in ascending order.
C
C        Z  contains the eigenvectors if MATZ is not zero.
C
C        IERR  is an integer output variable set equal to an
C        error completion code described in section 2B of the
C        documentation.  The normal completion code is zero.
C
C        FV1  and  FV2  are temporary storage arrays.
C
      INTEGER N,NM,IERR,MATZ,I,J,NN1,L,K
      REAL*8 A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N),A1
C
      IF (N .LE. NM) GO TO 20
      IERR = 10 * N
      GO TO 50
C
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
C****************************************S.C.
   20 CALL TRID(N,A,NM,FV1,W)
      DO I=1,N
        DO J=1,N
          Z(J,I)=A(J,I)
        END DO
      END DO
      NN1=N-1
      DO 333 I=1,NN1
      L=I+1
      DO 33 J=L,N
      IF(W(I)-W(J))3,3,1
   1  CONTINUE
      A1=W(I)
      W(I)=W(J)
      W(J)=A1
      DO 2 K=1,N
      A1=Z(K,I)
      Z(K,I)=Z(K,J)
      Z(K,J)=A1
   2  CONTINUE
   3  CONTINUE
  33  CONTINUE
 333  CONTINUE
C****************************************S.C.
   50 RETURN
      END
