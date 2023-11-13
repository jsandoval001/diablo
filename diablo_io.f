C This file contains subroutines for inputting and outputting data in
C Diablo as well as all subroutines called directly from diablo.f
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INITIALIZE
      INCLUDE 'header'

      REAL    VERSION, CURRENT_VERSION
      logical RESET_TIME
      INTEGER I, J, K, N

      WRITE(*,*) 'here'

      CURRENTFILE = INPATH(:LIP)//'input.dat'
      WRITE(*,*) 'file ',CURRENTFILE
      OPEN (11,file=CURRENTFILE,form='formatted',status='old')      
      WRITE(*,*) 'INPUT.DAT OPENED'

C Read input file.
C   (Note - if you change the following section of code, update the
C    CURRENT_VERSION number to make obsolete previous input files!)

      CURRENT_VERSION=2.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) FLAVOR,   VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input data format.'
      READ(11,*)
      READ(11,*) USE_MPI,    LES,    READ_HDF5,    SAVE_HDF5
      READ(11,*)
      READ(11,*) RE, LX, LY, LZ, TIME_LIMIT
      READ(11,*)
      READ(11,*) NUM_PER_DIR, CREATE_NEW_FLOW, RHO_TYPE, TIME_AVERAGE
      READ(11,*)
      READ(11,*) N_TIME_STEPS, DELTA_T, RESET_TIME, VARIABLE_DT, CFL
     &            , UPDATE_DT
      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, MOVIE, INIT_E
     &     , T0
      READ(11,*)
! Read in the parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) FILTER_TH(N), FILTER_INT(N)
        READ(11,*)
        READ(11,*) RI_TAU(N), PR(N), REACTION(N)
      END DO

	NU = 1.0d0/RE

C If we are using MPI, then Initialize the MPI Variables
        WRITE(*,*) 'ABOUT TO INITIALIZE MPI'
      IF (USE_MPI) THEN
        CALL INIT_MPI
      ELSE
        RANK_G = 0
      END IF
      WRITE(*,*) 'MPI INITIALIZED'

      IF (RANK_G.EQ.0) THEN
          WRITE(*,*) 'RHO_TYPE: ',RHO_TYPE
      END IF
	
C Initialize channel package.

        CALL INPUT_CHAN
        CALL CREATE_GRID_CHAN
        IF (USE_MPI) THEN
          CALL INIT_CHAN_MPI
        ELSE 
          CALL INIT_CHAN
        END IF
     

C Initialize grid
      IF (FLAVOR.NE.'Ensemble') THEN
      IF (RANK_G.EQ.0) THEN
	WRITE(6,*) 'Note that this code is distributed under the ',
     *           'GNU General Public License.'
      WRITE(6,*) 'No warranty is expressed or implied.'
      WRITE(6,*)
      write(*,*) 'Flavor: ',FLAVOR
      WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
      WRITE(6,*) 'RE: ',RE,', NU: ',NU
      DO N=1,N_TH
        WRITE(6,*) 'Scalar number: ',N
        WRITE(6,*) '  Richardson number: ',RI_TAU(N)
        WRITE(6,*) '  Prandlt number: ',PR(N)
      END DO
      END IF
	END IF
      NXM=NX-1
      NYM=NY-1
      NZM=NZ-1

C Initialize storage arrays.
      DO K=0,NZ+1
        DO I=0,NX+1 
          DO J=0,NY+1
            U1(I,K,J)=0.d0
            U3(I,K,J)=0.d0
            U2(I,K,J)=0.d0
            P (I,K,J)=0.d0
            R1(I,K,J)=0.d0
            R2(I,K,J)=0.d0
            R3(I,K,J)=0.d0
            F1(I,K,J)=0.d0
            F2(I,K,J)=0.d0
            F3(I,K,J)=0.d0
            BU1(I,K,J)=0.d0
            DBU1(I,K,J) = 0.D0
            BU1F(I,K,J)=0.d0
            DO N = 1,N_TH
            BTH(I,K,J,N)=0.d0
            END DO
          END DO
        END DO
      END DO

C Initialize FFT package (includes defining the wavenumber vectors).
      CALL INIT_FFT
      
      ENERGIES(:) = 0.D0
      DISSIPATIONS(:) = 0.D0
      ENERGYCOUNT = 1
      DISSCOUNT = 1
	
C Initialize RKW3 parameters.
      H_BAR(1)=DELTA_T*(8.0d0/15.0d0)
      H_BAR(2)=DELTA_T*(2.0d0/15.0d0)
      H_BAR(3)=DELTA_T*(5.0d0/15.0d0)
      BETA_BAR(1)=1.0d0
      BETA_BAR(2)=25.0d0/8.0d0
      BETA_BAR(3)=9.0d0/4.0d0
      ZETA_BAR(1)=0.0d0
      ZETA_BAR(2)=-17.0d0/8.0d0
      ZETA_BAR(3)=-5.0d0/4.0d0
      H_ALPHA(1)=DELTA_T*(0.D0)
      H_ALPHA(2)=DELTA_T*(8.D0/15.D0)
      H_ALPHA(3)=DELTA_T*(2.D0/3.D0)
      DISSIPATION = 0.D0
      NSAMPLES = 0

C Initialize values for reading of scalars
      NUM_READ_TH=0
      DO N=1,N_TH
        IF (CREATE_NEW_TH(N)) THEN
          NUM_READ_TH=NUM_READ_TH 
        ELSE
          NUM_READ_TH=NUM_READ_TH + 1
          READ_TH_INDEX(NUM_READ_TH)=N
        END IF
      END DO
	
      CALL CREATE_BACKGROUND_TH_CHAN
      CALL CREATE_TH_CHAN 
	
C initialise adjoint system (to be done)

!	  CALL INIT_ADJOINT

C Create flow.
C Initialize flow.
      IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
        PREVIOUS_TIME_STEP=0
        TIME_STEP=0
        TIME=T0
      END IF

      TIME_TIME_STEP = TIME
	CALL CREATE_BACKGROUND_FLOW_CHAN
	
      IF (CREATE_NEW_FLOW) THEN
	  
          CALL CREATE_FLOW_CHAN
     
        IF (FLAVOR.NE.'Ensemble') THEN
        IF (RANK_G.EQ.0) THEN
	    write(*,*) 'A new flowfield has been created'
	  END IF
	  END IF
        IF (FLAVOR.EQ.'Basic') THEN
	    CALL SAVE_STATS(.FALSE.)
!	    CALL SAVE_FLOW(.FALSE.)
	  END IF
      ELSE
      IF (RANK_G.EQ.0) THEN
        write(*,*) 'Reading flow...'
      END IF
        CALL READ_FLOW
!        CALL PERTURB
      IF (RANK_G.EQ.0) THEN
        write(*,*) 'Done reading flow'
      END IF
      END IF
 
C Initialize flow.
      IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
        PREVIOUS_TIME_STEP=0
        TIME_STEP=0
        TIME=T0
      END IF


      IF (RANK_G.EQ.0) THEN
        WRITE(*,*) 'INIT_E', INIT_E
      END IF
      CALL GET_ENERGY(.TRUE.)
      CALL GET_ENERGY(.FALSE.)
      
        CALL SAVE_STATS(.FALSE.)

          CALL POISSON_P_CHAN

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      LOGICAL FINAL

        CALL SAVE_STATS_CHAN(FINAL)          
   
	
	if (flavor.eq.'Basic') then
	IF (RANK_G.EQ.0) THEN
        write(*,*) 'done save_stats diablo'
      END IF
	end if

      IF (FINAL) THEN
     
          CALL VIS_FLOW_CHAN         

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*95 FNAME
      CHARACTER*95 FNAME_TH(N_TH)
      INTEGER I, J, K, N, NUM_PER_DIR_T

#ifdef HDF5
      IF (READ_HDF5) THEN
        FNAME = INPATH(:LIP)//'start.h5'
      ELSE IF (USE_MPI) THEN
        FNAME = INPATH(:LIP)//'diablo_'//MPI_IO_NUM//'.start'
      ELSE
        FNAME=INPATH(:LIP)//'diablo.start'
      END IF
#else
      IF (USE_MPI) THEN
        FNAME = INPATH(:LIP)//'diablo_'//MPI_IO_NUM//'.start'
      ELSE
        FNAME=INPATH(:LIP)//'diablo.start'
      END IF
#endif

      IF (RANK_G.EQ.0) THEN
        WRITE(6,*) 'Reading flow from ',FNAME
      END IF
      IF (FNAME(LEN_TRIM(FNAME)-2:LEN_TRIM(FNAME)).EQ.".h5") THEN
#ifdef HDF5      
      CALL READHDF5(FNAME)
#else
         write(*,*) ' **** ERROR ******************************'
         write(*,*) ' Program not compiled with HDF5 libraries.'
         stop       
#endif      
      
      ELSE

      DO N=1,N_TH
        IF (USE_MPI) THEN
                FNAME_TH(N)=INPATH(:LIP)//'diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) //'_'//MPI_IO_NUM// '.start'
        ELSE
        FNAME_TH(N)=INPATH(:LIP)//'diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.start'
        END IF
      END DO
      WRITE(6,*)   'Reading flow from ',FNAME

      OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
      READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP


      write(*,*) 'NX_T, NY_T, NZ_T: ',NX_T,NY_T,NZ_T

! TEMPORARY, UNCOMMENT!!!
!      IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T))
!     *     STOP 'Error: old flowfield wrong dimensions. '
!      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T)
!     *     STOP 'Error: old flowfield wrong NUM_PER_DIR. '

      write(*,*) 'READING FLOW'

        READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CP(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1)

        DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
          READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX),K=0,TNKZ),J=0,NY+1)
         CLOSE(11)
        END DO
 
      CLOSE(10)
      CLOSE(11)
      END IF

C Apply initial boundary conditions, set ghost cells
      IF (USE_MPI) THEN
        call APPLY_BC_VEL_MPI
      ELSE
        call APPLY_BC_VEL_LOWER
        call APPLY_BC_VEL_UPPER
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*95 FNAME
      CHARACTER*95 FNAME_TH(N_TH)
      INTEGER      I, J, K, N
      LOGICAL      FINAL, FLAG
      INTEGER ID
      CHARACTER*4 STR

      IF (FINAL) THEN
#ifdef HDF5
            IF(SAVE_HDF5) THEN
              FNAME = SAVPATH(:LSP)//'end.h5'
            ELSE IF (USE_MPI) THEN
              FNAME = SAVPATH(:LSP)//'diablo_'//MPI_IO_NUM//'.res'
            ELSE
              FNAME=SAVPATH(:LSP)//'diablo.res'
            END IF
#else      
            IF (USE_MPI) THEN
               FNAME = SAVPATH(:LSP)//
     &              'diablo_'//MPI_IO_NUM//'.res'
            ELSE
               FNAME=SAVPATH(:LSP)//'diablo.res'
            END IF
#endif            
         
      ELSE
#ifdef HDF5
        IF (SAVE_HDF5) THEN
          FNAME = SAVPATH(:LSP)//'out.h5'                
        ELSE IF (USE_MPI) THEN
          FNAME = SAVPATH(:LSP)//
     &          'diablo_'//MPI_IO_NUM//'.saved'
        ELSE
          FNAME=SAVPATH(:LSP)//'diablo.saved'
        END IF
#else
        IF (USE_MPI) THEN
          FNAME = SAVPATH(:LSP)//
     &          'diablo_'//MPI_IO_NUM//'.saved'
        ELSE
          FNAME=SAVPATH(:LSP)//'diablo.saved'
        END IF            
#endif        
         
      END IF
      
      IF (FNAME(LEN_TRIM(FNAME)-2:LEN_TRIM(FNAME)).EQ.".h5") THEN
      
#ifdef HDF5
         inquire(file=SAVPATH(:LSP)//'save_flow.ind',exist=flag)
         if (flag) then
            open(unit=500,file=SAVPATH(:LSP)//
     &           'save_flow.ind',status='old'
     &           ,form='formatted')
            read(500,'(1I10)') id
            close(unit=500) 
            write(str,'(1I0.4)') id
            FNAME=SAVPATH(:LSP)//str//'_'//'out.h5'
         end if
            call WriteHDF5(FNAME)
         if (flag) then
            open(unit=500,file=SAVPATH(:LSP)//
     &           'save_flow.ind',form='formatted')
            write(500,'(1I10)') id+1
            close(unit=500)
         end if
#else
         write(*,*) ' **** ERROR ******************************'
         write(*,*) ' Program not compiled with HDF5 libraries.'
         stop
#endif      
      
      ELSE
       IF (FINAL) THEN
         DO N=1,N_TH
           IF (USE_MPI) THEN
             FNAME_TH(N) = SAVPATH(:LSP)//'diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) //'_'//MPI_IO_NUM// '.res'             
           ELSE
             FNAME_TH(N)=SAVPATH(:LSP)//'diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.res'
           END IF

         END DO      
      ELSE
      
         DO N=1,N_TH
           IF (USE_MPI) THEN
             FNAME_TH(N) = SAVPATH(:LSP)//'diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) //'_'//MPI_IO_NUM// '.saved'
           ELSE
             FNAME_TH(N)=SAVPATH(:LSP)//'diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.saved'
           END IF


         END DO      
      
      END IF
      
      
      IF (RANK_G.EQ.0) THEN
      WRITE(6,*) 'Writing flow to ',FNAME
      END IF

      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
      WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP



        WRITE (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CP(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,NY+1)
          CLOSE(11)
        END DO

      CLOSE(10)
      CLOSE(11)
      END IF

      RETURN
      END
      
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE GET_ENERGY(SET_ENERGY)
C----*|--.---------.---------.---------.---------.---------.---------.-|--      

      INCLUDE 'header'
      REAL*8 TEMP1(0:NY+1)
      REAL*8 TEMP2(0:NY+1)
      REAL*8 ENERGY,ANS,KINENERGY,POTENERGY
      LOGICAL, INTENT(IN) :: SET_ENERGY            
      INTEGER I,J,K,N
      COMPLEX*16 CTHTEMP(0:NKX,0:TNKZ,1:NY,1:N_TH)
      CHARACTER*150 FNAMET,FNAMEK,FNAMEP

      FNAMET = SAVPATH(:LSP)//'./tot.txt'     
      FNAMEK = SAVPATH(:LSP)//'./kin.txt'     
      FNAMEP = SAVPATH(:LSP)//'./pot.txt'     

      OPEN(UNIT=20,FILE=FNAMET,ACCESS='APPEND',STATUS='OLD')
      OPEN(UNIT=40,FILE=FNAMEK,ACCESS='APPEND',STATUS='OLD')
      OPEN(UNIT=60,FILE=FNAMEP,ACCESS='APPEND',STATUS='OLD')

      TEMP1(:) = 0.d0
      TEMP2(:) = 0.D0
      ENERGY = 0.d0
      ANS = 0.d0

      DO J = 1,NY
        DO K = 0,TNKZ
          TEMP1(J) = TEMP1(J) + LX*LZ
     &                  * CDABS(CU1(0,K,J))**2.0d0
          TEMP1(J) = TEMP1(J) + LX*LZ 
     &                  * CDABS(CU3(0,K,J))**2.0d0
          DO N = 1,N_TH
              TEMP2(J) = TEMP2(J) + RI_TAU(N) * LX*LZ 
     &         * ABS(CTHTEMP(0,K,J,N))**2.0d0
          END DO
          DO I = 1,NKX
            TEMP1(J) = TEMP1(J) + 2.0d0*LX*LZ
     &                  * CDABS(CU1(I,K,J))**2.0d0
            TEMP1(J) = TEMP1(J) + 2.0d0*LX*LZ
     &                  * CDABS(CU3(I,K,J))**2.0d0
            DO N = 1,N_TH
              TEMP2(J) = TEMP2(J) + RI_TAU(N) * 2.0d0*LX*LZ
     &         * CDABS(CTH(I,K,J,N))**2.0d0
            END DO
          END DO
        END DO
      END DO
      CALL TRAPZ(TEMP1,KINENERGY,1)
      CALL TRAPZ(TEMP2,POTENERGY,1)

      TEMP1(:) = 0.d0
      
      DO J = 1,NY
        DO K = 0,TNKZ
          TEMP1(J) = TEMP1(J) + LX*LZ
     &                  * CDABS(CU2(0,K,J))**2.0d0
          DO I = 1,NKX
            TEMP1(J) = TEMP1(J) + 2.0d0 * LX*LZ 
     &                  * CDABS(CU2(I,K,J))**2.0d0
          END DO
        END DO
      END DO
      CALL TRAPZ(TEMP1,ANS,2)
      
      KINENERGY = KINENERGY + ANS
      
      KINENERGY = KINENERGY/2.0d0   
      POTENERGY = POTENERGY/2.0D0
      
      KINENERGY = KINENERGY/(LX*LZ*LY)                                      
      POTENERGY = POTENERGY/(LX*LZ*LY)
             
      CALL GATHER_SCALAR_MPI(KINENERGY,KINENERGY)
      CALL GATHER_SCALAR_MPI(POTENERGY,POTENERGY)
      ENERGY = KINENERGY + POTENERGY
      
      IF (RANK_G.EQ.0) THEN
        WRITE(*,*) 'ENERGY:', ENERGY
        WRITE(20,*) ENERGY
        WRITE(40,*) KINENERGY
        WRITE(60,*) POTENERGY                       
      END IF

      CLOSE(20)
      CLOSE(40)
      CLOSE(60)
      
      IF (SET_ENERGY) THEN      
      
        CU1(:,:,:) = CU1(:,:,:)*DSQRT(INIT_E/ENERGY)
        CU2(:,:,:) = CU2(:,:,:)*DSQRT(INIT_E/ENERGY)
        CU3(:,:,:) = CU3(:,:,:)*DSQRT(INIT_E/ENERGY)
        CTH(:,:,:,:) = CTH(:,:,:,:)*DSQRT(INIT_E/ENERGY)      
                         
        IF (USE_MPI) THEN
          CALL GHOST_CHAN_MPI
        END IF                           
                         
        CALL POISSON_P_CHAN
      
        IF (USE_MPI) THEN
          CALL GHOST_CHAN_MPI
        END IF                         
        
      ELSE               
        
        POTS(ENERGYCOUNT) = POTENERGY
        KINS(ENERGYCOUNT) = KINENERGY
        ENERGIES(ENERGYCOUNT) = ENERGY
        ENERGYCOUNT = ENERGYCOUNT + 1
                         
      END IF
                               
      RETURN
      END 
      
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE GET_DISSIPATION(LAST_TIME)
C----*|--.---------.---------.---------.---------.---------.---------.-|--

! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

      include 'header'
      character*150 FNAME
      LOGICAL LAST_TIME      
      REAL*8 NEW_DISS
      REAL*8 ANS
      REAL*8 VECTOR(0:NY+1)

      integer I,J,K,N

      IF (LAST_TIME) THEN                                          
                       
      IF (RANK_G.EQ.0) THEN     
        WRITE(*,*) 'Dissipation = ',DISSIPATION        
      END IF
            
      ELSE

      FNAME = SAVPATH(:LSP)//'diss.txt'

      OPEN(UNIT=80,FILE=FNAME,ACCESS='APPEND',STATUS='OLD')

      NEW_DISS = 0.d0
      ANS = 0.d0
      VECTOR(:) = 0.d0
      NSAMPLES = NSAMPLES + 1

! DU1DX
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J)+NU*LX*LZ*KX2(0)*CDABS(CU1(0,K,J))**2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KX2(I)*CDABS(CU1(I,K,J))**2.0d0
          END DO
        END DO
      END DO
      CALL TRAPZ(VECTOR,ANS,1)
      NEW_DISS = NEW_DISS + ANS
      
! DU1DY
      VECTOR(:) = 0.d0      
      
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ
     &            *(CDABS(CU1(0,K,J)-CU1(0,K,J-1))**2.0d0)/(DY(J)*DY(J))
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 2.0d0*NU*LX*LZ
     &            *(CDABS(CU1(I,K,J)-CU1(I,K,J-1))**2.0d0)/(DY(J)*DY(J))                    
          END DO
        END DO            
      END DO       
      CALL TRAPZ(VECTOR,ANS,2)    
      NEW_DISS = NEW_DISS + ANS      
      
! DU1DZ
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ*KZ2(K)*CDABS(CU1(0,K,J))
     &          **2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KZ2(K)*CDABS(CU1(I,K,J))**2.0d0
          END DO         
        END DO
      END DO 
      CALL TRAPZ(VECTOR,ANS,1)
      NEW_DISS = NEW_DISS + ANS
      
! DU3DX
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ*KX2(0)*CDABS(CU3(0,K,J))
     &          **2.0d0
          DO I = 1,NKX    
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KX2(I)*CDABS(CU3(I,K,J))**2.0d0
          END DO          
        END DO
      END DO
      CALL TRAPZ(VECTOR,ANS,1)
      NEW_DISS = NEW_DISS + ANS
      
! DU3DY
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ
     &            *(CDABS(CU3(0,K,J)-CU3(0,K,J-1))**2.0d0)/(DY(J)*DY(J))
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 2.0d0*NU*LX*LZ
     &            *(CDABS(CU3(I,K,J)-CU3(I,K,J-1))**2.0d0)/(DY(J)*DY(J))
          END DO         
        END DO
      END DO  
      CALL TRAPZ(VECTOR,ANS,2)      
      NEW_DISS = NEW_DISS + ANS    
      
! DU3DZ
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ*KZ2(K)*CDABS(CU3(0,K,J))
     &          **2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KZ2(K)*CDABS(CU3(I,K,J))**2.0d0
          END DO          
        END DO
      END DO  
      CALL TRAPZ(VECTOR,ANS,1)
      NEW_DISS = NEW_DISS + ANS
      
! DU2DX
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + 
     &             NU*LX*LZ*KX2(0)*CDABS(CU2(0,K,J))**2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KX2(I)*CDABS(CU2(I,K,J))**2.0d0
          END DO                    
        END DO
      END DO
      CALL TRAPZ(VECTOR,ANS,2)
      NEW_DISS = NEW_DISS + ANS
      
! DU2DY
      DO J = 1,NY
      VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + NU*LX*LZ
     &            *(CDABS(CU2(0,K,J+1)-CU2(0,K,J))**2.0d0)/
     &     (DYF(J)*DYF(J))
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 2.0d0*NU*LX*LZ
     &            *(CDABS(CU2(I,K,J+1)-CU2(I,K,J))**2.0d0)/
     &     (DYF(J)*DYF(J))
          END DO         
        END DO
      END DO           
      CALL TRAPZ(VECTOR,ANS,1) 
      NEW_DISS = NEW_DISS + ANS
      
! DU2DZ
      DO J = 1,NY
        VECTOR(J) = 0.d0
        DO K = 0,TNKZ
          VECTOR(J) = VECTOR(J) + 
     &             NU*LX*LZ*KZ2(K)*CDABS(CU2(0,K,J))**2.0d0
          DO I = 1,NKX
            VECTOR(J) = VECTOR(J) + 
     &             2.0d0*NU*LX*LZ*KZ2(K)*CDABS(CU2(I,K,J))**2.0d0
          END DO          
        END DO
      END DO  
      CALL TRAPZ(VECTOR,ANS,2)   
      NEW_DISS = NEW_DISS + ANS                            
      
! DTHDX
!      DO N = 1,N_TH
!      DO J = 1,NY
!        VECTOR(J) = 0.d0
!        DO K = 0,TNKZ
!          VECTOR(J) = VECTOR(J) + 
!     &       NU*LX*LZ*KX2(0)*CDABS(CTH(0,K,J,N))**2.0d0
!     &            *RI_TAU(N)/PR(N)
!          DO I = 1,NKX
!            VECTOR(J) = VECTOR(J) + 
!     &       2.0d0*NU*LX*LZ*KX2(I)*CDABS(CTH(I,K,J,N))**2.0d0
!     &            *RI_TAU(N)/PR(N)
!          END DO         
!        END DO
!      END DO   
!      CALL TRAPZ(VECTOR,ANS,1) 
!      NEW_DISS = NEW_DISS + ANS
!      END DO
      
! DTHDY
!      DO N = 1,N_TH
!      DO J = 1,NY
!        VECTOR(J) = 0.d0
!        DO K = 0,TNKZ
!          VECTOR(J) = VECTOR(J) + 
!     &             NU*LX*LZ
!     &          *(CDABS(CTH(0,K,J,N)-CTH(0,K,J-1,N))**2.0d0)/
!     &           (DY(J)*DY(J))
!     &            *RI_TAU(N)/PR(N)
!          DO I = 1,NKX
!            VECTOR(J) = VECTOR(J) + 
!     &             2.0d0*NU*LX*LZ
!     &          *(CDABS(CTH(I,K,J,N)-CTH(I,K,J-1,N))**2.0d0)/
!     &           (DY(J)*DY(J))
!     &            *RI_TAU(N)/PR(N)
!          END DO          
!        END DO
!      END DO  
!      CALL TRAPZ(VECTOR,ANS,2)   
!      NEW_DISS = NEW_DISS + ANS 
!      END DO
      
! DTHDZ
!      DO N = 1,N_TH
!      DO J = 1,NY
!        VECTOR(J) = 0.d0
!        DO K = 0,TNKZ
!          VECTOR(J) = VECTOR(J) + 
!     &             NU*LX*LZ*KZ2(K)*CDABS(CTH(0,K,J,N))**2.0d0
!     &            *RI_TAU(N)/PR(N)
!          DO I = 1,NKX
!            VECTOR(J) = VECTOR(J) + 
!     &             2.0d0*NU*LX*LZ*KZ2(K)*CDABS(CTH(I,K,J,N))**2.0d0
!     &            *RI_TAU(N)/PR(N)
!          END DO          
!        END DO
!      END DO    
!      CALL TRAPZ(VECTOR,ANS,1)
!      NEW_DISS = NEW_DISS + ANS
!      END DO      
      
      NEW_DISS = NEW_DISS / (LX*LY*LZ)
      
      IF (USE_MPI) THEN
        CALL GATHER_SCALAR_MPI(NEW_DISS,NEW_DISS)
      END IF      
      
      IF (RANK_G.EQ.0) THEN
        WRITE(80,*) NEW_DISS
      END IF
      DISSIPATIONS(DISSCOUNT) = NEW_DISS
      DISSCOUNT = DISSCOUNT + 1
      
      CLOSE(80)

! TRAPEZOIDAL RULE IN TIME      
      IF (TIME_AVERAGE) THEN
      IF (NSAMPLES.EQ.1) THEN
      DISSIPATION = DISSIPATION + 0.5d0 * NEW_DISS / FLOAT(N_TIME_STEPS)
      ELSEIF (NSAMPLES.EQ.N_TIME_STEPS) THEN
      DISSIPATION = DISSIPATION + 0.5d0 * NEW_DISS / FLOAT(N_TIME_STEPS)
      ELSE
        DISSIPATION = DISSIPATION + NEW_DISS / FLOAT(N_TIME_STEPS)
      END IF
      ELSE
        DISSIPATION = NEW_DISS
      END IF                  
      
! TIME AVERAGED
!      DISSIPATION = (1.d0/FLOAT(NSAMPLES)*NEW_DISS)
!     &      +((FLOAT(NSAMPLES)-1.d0)/FLOAT(NSAMPLES))*DISSIPATION
      
      END IF
      
      RETURN
      END   
      
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|      
      SUBROUTINE END_RUN(FLAG)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
      INCLUDE 'header'

      LOGICAL FLAG,FILE_EXISTS
      
      FLAG=.FALSE.
      ! Check for the time
      call WALL_TIME(END_TIME)
      if (END_TIME-START_TIME.gt.TIME_LIMIT) THEN
         write(*,*) ' STOP beacuse of wall-time hit!'
         FLAG=.TRUE.
      END IF
      
      INQUIRE(FILE="stop.now", EXIST=FILE_EXISTS)
      IF ( FILE_EXISTS ) THEN
         write(*,*) ' STOP beacuse of stop.now file!'
         FLAG=.TRUE.
      END IF
      
      RETURN
      END 


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine wall_time(wt)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c
c     Return wall-clock time as seconds after Jan. 1, 2014.
c     Support for leap year is not included anymore.
c
c     Next leap year is 2016!
c
c     By using a 'save' statement, the wall-time after the first
c     call to the subroutine could be computed, but that is not
c     intended with the present subroutine (e.g. the history file)
c
      implicit none

      real*8 wt
      integer val(8),i,shift,day

      integer mon(12,2)
      data mon /
     &     31,28,31,30,31,30,31,31,30,31,30,31,
     &     31,29,31,30,31,30,31,31,30,31,30,31/
c
c     Get current date and time
c     val(1) : year
c     val(2) : month
c     val(3) : day
c     val(4) : difference to GMT
c     val(5) : hour
c     val(6) : minute
c     val(7) : second
c     val(8) : 1/1000 second
c
      call date_and_time(values=val)
c
c     Determine leap year
c
      if (mod(val(1),4).eq.0) then
         if (mod(val(1),100).eq.0) then
            if (mod(val(1),400).eq.0) then
               shift=2
            else
               shift=1
            end if
         else
            shift=2
         end if
      else
         shift = 1
      end if
c
c     Construct day of the year
c
      day = val(3)-1
      do i=1,val(2)-1
         day=day+mon(i,shift)
      end do
c
c     And compute wall-clock time
c
      wt = (val(1)-2014)*365*86400+
     &     day*86400+val(5)*3600+val(6)*60+val(7)+val(8)/1000.

      end                
