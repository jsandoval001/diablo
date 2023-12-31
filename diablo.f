C******************************************************************************|
C diablo.f -> DNS In A Box, Laptop Optimized                       VERSION 0.9
C
C This Fortran 77 code computes incompressible flow in a box.
C
C Primative variables (u,v,w,p) are used, and continuity is enforced with a
C fractional step algorithm.
C
C SPATIAL DERIVATIVES:
C   0, 1, 2, or 3 directions are taken to be periodic and handled spectrally
C   (these cases are referred to as the "periodic", "channel", "duct", and
C    "cavity" cases respectively).
C   The remaining directions are taken to be bounded by walls and handled with
C   momentum- and energy-conserving second-order central finite differences.
C
C TIME ADVANCEMENT
C   Two main approaches are implemented:
C     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
C     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
C
C The emphasis in this introductory code is on code simplicity:
C   -> All variables are in core.
C   -> The code is not explicitly designed for use with either MPI or SMP.
C   -> Overindexing is not used.
C A few simple high-performance programming constructs are used:
C   -> The inner 2 loops are broken out in such a way as to enable out-of-order
C      execution of these loops as much as possible, thereby leveraging
C      vector and superscalar CPU architectures.
C   -> The outer loops are fairly long (including as many operations as
C      possible inside on a single J plane of data) in order to make effective
C      use of cache.
C Multiple time advancement algorithms are implemented for the periodic,
C channel, duct, and cavity cases in order to compare their efficiency for
C various flows on different computational architectures.  In a few of the
C algorithms, the code is broken out fully into a PASS1/PASS2 architecture
C to maximize the efficient use of cache.
C
C This code was developed as a joint project in MAE 223 (CFD), taught by
C Thomas Bewley, at UC San Diego (spring 2001, spring 2005).
C Primary contributions follow:
C Thomas Bewley was the chief software architect
C John R. Taylor wrote the channel flow solvers
C******************************************************************************|
C
C This code is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version. This code is distributed in the hope that it
C will be useful, but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C GNU General Public License for more details. You should have received a
C copy of the GNU General Public License along with this code; if not,
C write to the Free Software Foundation, Inc., 59 Temple Place - Suite
C 330, Boston, MA 02111-1307, USA.
C
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      PROGRAM DIABLO
!      implicit none
      INCLUDE 'header'
      INTEGER N
      LOGICAL FLAG
      

      WRITE(6,*) 
      WRITE(6,*) '             ****** WELCOME TO DIABLO ******'
      WRITE(6,*) '                      VERTICAL CODE'
      WRITE(6,*)

!      OPEN(11,FILE='/media/data-2/stokes_edge_tracking/wall/Re_1000/'//
!     &     'first_guests_2023_11_03/005_1e-08/path.dat'
!     &     ,FORM='FORMATTED',
!     &     STATUS='OLD')
      OPEN(11,FILE='path.dat',FORM='FORMATTED', STATUS='OLD')
      READ(11,'(A)') INPATH
      READ(11,'(A)') SAVPATH

      CLOSE(11)
      WRITE(*,*) 'INPATH: ',INPATH
      WRITE(*,*) 'OUTPATH: ',SAVPATH
      
      LIP = LEN_TRIM(INPATH)
      LSP = LEN_TRIM(SAVPATH)

      CALL INITIALIZE

! COMPUTE TIME TAKEN
      CALL WALL_TIME(START_TIME)

C A flag to determine if we are considering the first time-step
      FIRST_TIME=.TRUE.

      CALL SAVE_FLOW(.FALSE.)

      DO TIME_STEP = TIME_STEP+1, TIME_STEP+N_TIME_STEPS
      IF (RANK_G.EQ.0) THEN
        WRITE(6,*) 'Now beginning TIME_STEP = ',TIME_STEP
      END IF

!        CALL GET_ENERGY(.FALSE.)

        DO RK_STEP=1,3                            
            CALL RK_CHAN_1         
        END DO

        IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
          CALL GET_DISSIPATION(.FALSE.)      
          CALL GET_ENERGY(.FALSE.)
        END IF
               
        TIME=TIME+DELTA_T
        FIRST_TIME=.FALSE.
! Save statistics to an output file
!        IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
!            CALL SAVE_STATS(.FALSE.)
!        END IF
        
! Save the flow to a restart file
        IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
          CALL SAVE_FLOW(.FALSE.)
!          CALL GET_ENERGY(.FALSE.)
        END IF
        
! Save the full three-dimensional fields in NETCDF format to vis.nc
! (optional)
        IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
          CALL NETCDF_OPEN_VIS
          CALL NETCDF_WRITE_VIS
          CALL NETCDF_CLOSE_VIS
        END IF

! Filter the scalar field
!      DO N = 1,N_TH
!        IF (FILTER_TH(N).AND.(MOD(TIME_STEP,FILTER_INT(N)).EQ.0)) THEN
!          IF (RANK_G.EQ.0) THEN    
!          write(*,*) 'Filtering...'
!          END IF
!          CALL FILTER_CHAN
!        END IF 
!      END DO

        IF (USE_MPI) THEN
          CALL END_RUN_MPI(FLAG)
        ELSE
          CALL END_RUN(FLAG)
        END IF
        IF (FLAG) THEN
          EXIT
        END IF

      END DO

! COMPUTE TIME TAKEN
      CALL GET_DISSIPATION(.TRUE.)
      CALL GET_ENERGY(.FALSE.)
      CALL WALL_TIME(END_TIME)
     
      IF (RANK_G.EQ.0) THEN
      WRITE(*,*) 'Elapsed Time (sec): ',END_TIME-START_TIME
      WRITE(*,*) 'Seconds per Iteration: '
     &     ,(END_TIME-START_TIME)/N_TIME_STEPS
      END IF
      
      IF (RANK_G.EQ.0) THEN
!        WRITE(*,*) 'DISSCOUNT ',DISSCOUNT
!        DO N = 1,DISSCOUNT
!          WRITE(*,*) 'DISSIPATIONS ',DISSIPATIONS(N)
!        END DO
!        WRITE(*,*) 'ENERGYCOUNT ',ENERGYCOUNT
!        DO N = 1,ENERGYCOUNT
!          WRITE(*,*) 'POTS ',POTS(N)
!          WRITE(*,*) 'KINS ',KINS(N)
!        END DO
      END IF
      
      IF (USE_MPI) THEN
        CALL WAIT
      END IF
      
      TIME_STEP=TIME_STEP-1
      CALL SAVE_FLOW(.TRUE.)
      CALL SAVE_STATS(.TRUE.)
      WRITE(6,*)
      WRITE(6,*) '        ****** Hello world!  Have a nice day! ******'
      WRITE(6,*)
      END


