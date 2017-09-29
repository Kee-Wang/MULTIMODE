C**************************************************************
C**************************************************************
      SUBROUTINE USERIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FILASS/IOUT,INP
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      RETURN
      END
C******************************************
C******************************************
      SUBROUTINE GETPOT
      RETURN
      END
C******************************************
C******************************************
      SUBROUTINE GETQPT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETAPT(VDP,NATOM,XX,RR,QQ,NMODE,ITAU)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION QQ(NMODE),XX(NATOM,3),RR(NATOM,NATOM)
      COMMON/DAVE/VE,V(722),WE(12),OMEGA(11,722),DERIV(11,722)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      IF(ITAU .EQ. 0) THEN
        VDP=0
        DO I = 1,12
          VDP = VDP + 0.5D0*WE(I)**2*QQ(I)**2
        ENDDO
      ELSE
        VDP=V(ITAU)-VE
        DO I = 1,11
          VDP = VDP + 0.5D0*OMEGA(I,ITAU)**2*QQ(I)**2
     1                     +DERIV(I,ITAU)*QQ(I)
        ENDDO
      ENDIF
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDIP
      RETURN
      END
C******************************************
C******************************************
      SUBROUTINE GETQDT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE MINPOT
      RETURN
      END
