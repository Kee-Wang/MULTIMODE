C*******************************************************
C*******************************************************
C**TIMER ROUTINES
C*******************************************************
C*******************************************************
C######################################################################
C     SUBROUTINE TIMIT(INDEX)
C     ************** TIMER ROUTINE FOR THE CRAY ************
C     COMMON/XIMES/TI,TX,TIM,TTGO,TSTART,TTGOO
C     COMMON/OPTNS/TIMLIM
C     COMMON/FILASS/IOUT
C     GOTO (1,2,3,1),INDEX
C2     TIM=0.0E0
C     TX=0.0E0
C     TI=0.0E0
C     DAY=DATE()
C     HOUR=CLOCK()
C     CALL TREMAIN(TTGO)
C     CALL SECOND(TSTART)
C     TIMLIM=TTGO
C     TTGOO=TIMLIM
C     WRITE(IOUT,99) DAY,HOUR,TIMLIM
C99    FORMAT(1X,'JOB STARTED ON ',A10,' AT TIME ',A10
C    1   //1X,' TIME LIMIT  ',F12.2,' SECONDS ')
C     RETURN
C1     CALL SECOND(TTT)
C     TIM=TTT-TSTART
C     TX=TIM-TI
C     TI=TIM
C     RETURN
C3     CALL SECOND(TTT)
C     TIM=TTT-TSTART
C     TX=TIM-TI
C     TI=TIM
C     TTGO=TIMLIM-TIM
C4     IF(INDEX.EQ.4) RETURN
C     WRITE(IOUT,1000) TX,TTGO
C1000  FORMAT(
C    &1X,'Last step ',F8.2,' secs : remaining ',F8.2,' secs')
C     RETURN
C     END
C######################################################################
      SUBROUTINE TIMIT(INDEX)
C**************** Timer for SGI *******************
C      EXTERNAL FDATE,DTIME
C      REAL*4 DTIME
      REAL*4 TARRAY(2),ETIME,TTT
      character*24 fdate
      COMMON/XIMES/TI,TX,TIM,TTGO,TSTART,TTGOO
      COMMON/OPTNS/TIMLIM
      COMMON/FILASS/IOUT
      GOTO (1,2,3,4),INDEX
 1     CONTINUE
      call DTIME(TARRAY,TTT)
      RETURN
 2     WRITE(IOUT,99) FDATE()
 99    FORMAT(1X,'JOB START AT:',A24)
      RETURN
 3     CONTINUE
      call DTIME(TARRAY,TTT)
      WRITE(IOUT,1000) TTT,TARRAY(1),TARRAY(2)
      WRITE(IOUT,1001) ETIME(TARRAY)
 1000  FORMAT(1X,'Last step ',F10.2,
     &       ' secs  ( user:',F10.2,' s and sys:',F10.2,' s)')
 1001  FORMAT(1X,'Total elapsed time', F10.2,' s.')
      RETURN
 4     call DTIME(TARRAY,TTT)
      WRITE(IOUT,1001) ETIME(TARRAY)
      RETURN
      END
C######################################################################
C     SUBROUTINE TIMIT(INDEX)
C**************** Timer for IBM *******************
C     CHARACTER*26 STR
C     REAL*4 DELTA
C     REAL*4 ELAPSED
C     TYPE TB_TYPE
C       SEQUENCE
C       REAL*4 USRTIME
C       REAL*4 SYSTIME
C     END TYPE
C     TYPE(TB_TYPE) DTIME_STRUCT
C     TYPE(TB_TYPE) ETIME_STRUCT
C     COMMON/FILASS/IOUT
C     GOTO (1,2,3,4),INDEX
C2     CALL fdate_(STR)
C     WRITE(IOUT,99) STR
C99    FORMAT(1X,'JOB START AT:',A26)
C     RETURN
C1     DELTA=dtime_(DTIME_STRUCT)
C     RETURN
C3     DELTA=dtime_(DTIME_STRUCT)
C     WRITE(IOUT,1000) DELTA
C1000  FORMAT(1X,'Last step ',F8.2,' secs',/)
C     RETURN
C4     ELAPSED=etime_(ETIME_STRUCT)
C     WRITE(IOUT,1001) ELAPSED
C1001  FORMAT(1X,'Total elapsed time', F8.2,' secs',/)
C     RETURN
C     END
