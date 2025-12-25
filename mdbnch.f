C***********************************************************************
C     MDBNCH - Molecular Dynamics Benchmark
C     A simple molecular dynamics simulation benchmark
C     Written in Fortran 77
C***********************************************************************
      PROGRAM MDBNCH
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=864)
      COMMON /PARTCL/ X(NMAX),Y(NMAX),Z(NMAX)
      COMMON /VLCTY/  VX(NMAX),VY(NMAX),VZ(NMAX)
      COMMON /FORCES/ FX(NMAX),FY(NMAX),FZ(NMAX)
      COMMON /PARAMS/ N,DT,SIDE,RCOFF
      
      N = 864
      DT = 0.001D0
      SIDE = 6.8D0
      RCOFF = 2.5D0
      
      WRITE(*,*) 'MOLECULAR DYNAMICS BENCHMARK'
      WRITE(*,*) 'Number of particles: ', N
      WRITE(*,*) 'Timestep: ', DT
      WRITE(*,*) 'Box side: ', SIDE
      WRITE(*,*) 'Cutoff radius: ', RCOFF
      WRITE(*,*)
      
      CALL INITPO
      CALL INITVL
      
      T1 = SECOND()
      
      DO 100 ISTEP = 1, 100
        CALL CALCFO
        CALL MOVEA
        CALL MOVEB
 100  CONTINUE
      
      T2 = SECOND()
      
      EKIN = VKCALC()
      TEMP = EKIN * 2.0D0 / (3.0D0 * DBLE(N))
      
      WRITE(*,*) 'Simulation completed'
      WRITE(*,*) 'Final kinetic energy: ', EKIN
      WRITE(*,*) 'Final temperature: ', TEMP
      WRITE(*,*) 'Elapsed time: ', T2-T1, ' seconds'
      
      STOP
      END

C***********************************************************************
C     Initialize positions on FCC lattice
C***********************************************************************
      SUBROUTINE INITPO
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=864)
      COMMON /PARTCL/ X(NMAX),Y(NMAX),Z(NMAX)
      COMMON /PARAMS/ N,DT,SIDE,RCOFF
      
      NCELL = INT(DBLE(N/4)**(1.0D0/3.0D0) + 0.5D0)
      DCELL = SIDE / DBLE(NCELL)
      
      ITEL = 0
      DO 30 LZ = 1, NCELL
        DO 30 LY = 1, NCELL
          DO 30 LX = 1, NCELL
            X0 = DBLE(LX-1) * DCELL
            Y0 = DBLE(LY-1) * DCELL
            Z0 = DBLE(LZ-1) * DCELL
            
            ITEL = ITEL + 1
            IF (ITEL .GT. N) GOTO 40
            X(ITEL) = X0
            Y(ITEL) = Y0
            Z(ITEL) = Z0
            
            ITEL = ITEL + 1
            IF (ITEL .GT. N) GOTO 40
            X(ITEL) = X0 + 0.5D0*DCELL
            Y(ITEL) = Y0 + 0.5D0*DCELL
            Z(ITEL) = Z0
            
            ITEL = ITEL + 1
            IF (ITEL .GT. N) GOTO 40
            X(ITEL) = X0 + 0.5D0*DCELL
            Y(ITEL) = Y0
            Z(ITEL) = Z0 + 0.5D0*DCELL
            
            ITEL = ITEL + 1
            IF (ITEL .GT. N) GOTO 40
            X(ITEL) = X0
            Y(ITEL) = Y0 + 0.5D0*DCELL
            Z(ITEL) = Z0 + 0.5D0*DCELL
 30   CONTINUE
      
 40   CONTINUE
      RETURN
      END

C***********************************************************************
C     Initialize velocities
C***********************************************************************
      SUBROUTINE INITVL
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=864)
      COMMON /VLCTY/  VX(NMAX),VY(NMAX),VZ(NMAX)
      COMMON /PARAMS/ N,DT,SIDE,RCOFF
      
      SUMVX = 0.0D0
      SUMVY = 0.0D0
      SUMVZ = 0.0D0
      
      DO 10 I = 1, N
        VX(I) = RANF() - 0.5D0
        VY(I) = RANF() - 0.5D0
        VZ(I) = RANF() - 0.5D0
        SUMVX = SUMVX + VX(I)
        SUMVY = SUMVY + VY(I)
        SUMVZ = SUMVZ + VZ(I)
 10   CONTINUE
      
      SUMVX = SUMVX / DBLE(N)
      SUMVY = SUMVY / DBLE(N)
      SUMVZ = SUMVZ / DBLE(N)
      
      DO 20 I = 1, N
        VX(I) = VX(I) - SUMVX
        VY(I) = VY(I) - SUMVY
        VZ(I) = VZ(I) - SUMVZ
 20   CONTINUE
      
      RETURN
      END

C***********************************************************************
C     Calculate forces using Lennard-Jones potential
C***********************************************************************
      SUBROUTINE CALCFO
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=864)
      COMMON /PARTCL/ X(NMAX),Y(NMAX),Z(NMAX)
      COMMON /FORCES/ FX(NMAX),FY(NMAX),FZ(NMAX)
      COMMON /PARAMS/ N,DT,SIDE,RCOFF
      
      RCOFFS = RCOFF * RCOFF
      
      DO 10 I = 1, N
        FX(I) = 0.0D0
        FY(I) = 0.0D0
        FZ(I) = 0.0D0
 10   CONTINUE
      
      DO 30 I = 1, N-1
        DO 20 J = I+1, N
          XX = X(I) - X(J)
          YY = Y(I) - Y(J)
          ZZ = Z(I) - Z(J)
          
          XX = XX - SIDE * DNINT(XX/SIDE)
          YY = YY - SIDE * DNINT(YY/SIDE)
          ZZ = ZZ - SIDE * DNINT(ZZ/SIDE)
          
          RD = XX*XX + YY*YY + ZZ*ZZ
          
          IF (RD .LE. RCOFFS) THEN
            R2I = 1.0D0 / RD
            R6I = R2I * R2I * R2I
            FF = 48.0D0 * R2I * R6I * (R6I - 0.5D0)
            
            FX(I) = FX(I) + FF * XX
            FY(I) = FY(I) + FF * YY
            FZ(I) = FZ(I) + FF * ZZ
            
            FX(J) = FX(J) - FF * XX
            FY(J) = FY(J) - FF * YY
            FZ(J) = FZ(J) - FF * ZZ
          ENDIF
 20     CONTINUE
 30   CONTINUE
      
      RETURN
      END

C***********************************************************************
C     First half of velocity Verlet integration
C***********************************************************************
      SUBROUTINE MOVEA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=864)
      COMMON /PARTCL/ X(NMAX),Y(NMAX),Z(NMAX)
      COMMON /VLCTY/  VX(NMAX),VY(NMAX),VZ(NMAX)
      COMMON /FORCES/ FX(NMAX),FY(NMAX),FZ(NMAX)
      COMMON /PARAMS/ N,DT,SIDE,RCOFF
      
      DO 10 I = 1, N
        VX(I) = VX(I) + 0.5D0 * DT * FX(I)
        VY(I) = VY(I) + 0.5D0 * DT * FY(I)
        VZ(I) = VZ(I) + 0.5D0 * DT * FZ(I)
        
        X(I) = X(I) + DT * VX(I)
        Y(I) = Y(I) + DT * VY(I)
        Z(I) = Z(I) + DT * VZ(I)
        
        X(I) = X(I) - SIDE * DNINT(X(I)/SIDE)
        Y(I) = Y(I) - SIDE * DNINT(Y(I)/SIDE)
        Z(I) = Z(I) - SIDE * DNINT(Z(I)/SIDE)
 10   CONTINUE
      
      RETURN
      END

C***********************************************************************
C     Second half of velocity Verlet integration
C***********************************************************************
      SUBROUTINE MOVEB
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=864)
      COMMON /VLCTY/  VX(NMAX),VY(NMAX),VZ(NMAX)
      COMMON /FORCES/ FX(NMAX),FY(NMAX),FZ(NMAX)
      COMMON /PARAMS/ N,DT,SIDE,RCOFF
      
      DO 10 I = 1, N
        VX(I) = VX(I) + 0.5D0 * DT * FX(I)
        VY(I) = VY(I) + 0.5D0 * DT * FY(I)
        VZ(I) = VZ(I) + 0.5D0 * DT * FZ(I)
 10   CONTINUE
      
      RETURN
      END

C***********************************************************************
C     Calculate kinetic energy
C***********************************************************************
      FUNCTION VKCALC()
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=864)
      COMMON /VLCTY/  VX(NMAX),VY(NMAX),VZ(NMAX)
      COMMON /PARAMS/ N,DT,SIDE,RCOFF
      
      SUMV = 0.0D0
      DO 10 I = 1, N
        SUMV = SUMV + VX(I)*VX(I) + VY(I)*VY(I) + VZ(I)*VZ(I)
 10   CONTINUE
      
      VKCALC = 0.5D0 * SUMV
      RETURN
      END

C***********************************************************************
C     Random number generator (0 to 1)
C***********************************************************************
      FUNCTION RANF()
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE ISEED
      DATA ISEED /12345/
      
      ISEED = MOD(ISEED * 1103515245 + 12345, 2147483647)
      RANF = DBLE(ISEED) / 2147483647.0D0
      
      RETURN
      END

C***********************************************************************
C     Timer function (returns CPU time in seconds)
C***********************************************************************
      FUNCTION SECOND()
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 ETIME, TARRAY(2)
      
      SECOND = DBLE(ETIME(TARRAY))
      
      RETURN
      END
