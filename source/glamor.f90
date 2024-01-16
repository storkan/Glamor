PROGRAM glamoR
  USE declarations
  USE gyro_mag_mod
  USE ReadControlFile_MOD
  USE ReadXYZFile_MOD
  USE ReadTimeFile_MOD
  USE POWDER_TYPE_DEC, ONLY: powder_type, CO_type
  USE LEB3J005_000084_MOD, ONLY:LEB3J005_000084
  USE LEB3J011_000600_MOD, ONLY:LEB3J011_000600
  USE LEB3J017_001980_MOD, ONLY:LEB3J017_001980
  USE LEB3J023_004656_MOD, ONLY:LEB3J023_004656
  USE LEB3J029_009060_MOD, ONLY:LEB3J029_009060
  USE LEB3J035_015624_MOD, ONLY:LEB3J035_015624
  USE LEB3J041_024780_MOD, ONLY:LEB3J041_024780

  IMPLICIT NONE
  !deklarationsteil


  ! variables in ReadControlFile_MOD
  CHARACTER(len=5):: recoupledisotope, observedisotope, observedatomtype, recoupledatomtype
  CHARACTER(len=30):: powderspec
  DOUBLE PRECISION:: convergence_criterion
  CHARACTER(len=30):: XYZfile, DYNfile

  ! program related variable, e.g variables that are necessary for reading and checking inputfiles etc.
  CHARACTER(len=30):: line
  INTEGER:: ios24=0
  INTEGER:: counter1=0,counter5=0 
  INTEGER:: filecheck, dynamiccheck
  INTEGER:: o, r, i !(schon in timesequence-loop und dann später in pow-av-loops) !, powder
  INTEGER:: j
  DOUBLE PRECISION:: RedCosinusTerm
  !########
  INTEGER :: DIM1 
  INTEGER, DIMENSION(:), POINTER :: tmpAtomNrlist
  DOUBLE PRECISION, DIMENSION(:), POINTER :: tmpWeightlist

  TYPE(powder_type), DIMENSION(:),POINTER :: powderav
  TYPE(xyz_static), DIMENSION(:), POINTER:: obs
  TYPE(xyz_dynamic), DIMENSION(:), POINTER:: recoup,voids
  !########

  ! calculation related variables
  ! FAIL

  DOUBLE PRECISION:: delta_x, delta_y, delta_z, gyro_ri, gyro_oi
  DOUBLE PRECISION, DIMENSION(:), POINTER::grid_array
  DOUBLE PRECISION :: D_k,vector_length,beta_pc_rad, gamma_pc_rad
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793d0, h_quer = 1.054571628d-34
  integer:: user=0

  !##########################
  ! variables needed in the kernel
  INTEGER :: powderstart, powderend, t
  DOUBLE PRECISION :: SinCosBeta_cr,SinGamma_cr,CosBeta_cr,SinBeta_cr,TwoSqrtTwoDwTOvThree,CosGamma_cr
  INTEGER :: p, powdervariable
  DOUBLE PRECISION :: OneminTwoCosBeta_crTwo,  Signal_Sk_In, Signal_S_In, Tau_k, C_1k, S_1k
  DOUBLE PRECISION :: Signal_S_I, Signal_pwdr_avrg
  DOUBLE PRECISION, DIMENSION (5,128) :: Signal_Sk_In_tau
  DOUBLE PRECISION:: kappa= 0.367553d0 
  DOUBLE PRECISION:: alp, bet, gam, wei
  Double precision:: convergencesum=0, convergencepoint
  Type DynamicAtom_TYPE
    Integer,DIMENSION(:),Pointer :: Pos
    INTEGER,DIMENSION(:),Pointer :: weight
!    DOUBLE PRECISION,DIMENSION(:),POINTER :: weight     
  END TYPE DynamicAtom_TYPE
  !###############################

  ! Programminternerteil
  ALLOCATE(powderav(8))

  DIM1=UBOUND(LEB3J005_000084,DIM=1)
  ALLOCATE(powderav(1)%CO(DIM1))
  powderav(1)%name="LEB3J005_000084"
  DO i=1,DIM1
     powderav(1)%CO(i)%a=LEB3J005_000084(i)%a 
     powderav(1)%CO(i)%b=LEB3J005_000084(i)%b 
     powderav(1)%CO(i)%g=LEB3J005_000084(i)%g 
     powderav(1)%CO(i)%w=LEB3J005_000084(i)%w 
  ENDDO

  DIM1=UBOUND(LEB3J011_000600,DIM=1)
  ALLOCATE(powderav(2)%CO(DIM1))
  powderav(2)%name="LEB3J011_000600"
  DO i=1,DIM1
     powderav(2)%CO(i)%a=LEB3J011_000600(i)%a 
     powderav(2)%CO(i)%b=LEB3J011_000600(i)%b 
     powderav(2)%CO(i)%g=LEB3J011_000600(i)%g 
     powderav(2)%CO(i)%w=LEB3J011_000600(i)%w 
  ENDDO

  DIM1=UBOUND(LEB3J017_001980,DIM=1)
  ALLOCATE(powderav(3)%CO(DIM1))
  powderav(3)%name="LEB3J017_001980"
  DO i=1,DIM1
     powderav(3)%CO(i)%a=LEB3J017_001980(i)%a 
     powderav(3)%CO(i)%b=LEB3J017_001980(i)%b 
     powderav(3)%CO(i)%g=LEB3J017_001980(i)%g 
     powderav(3)%CO(i)%w=LEB3J017_001980(i)%w 
  ENDDO

  DIM1=UBOUND(LEB3J023_004656,DIM=1)
  ALLOCATE(powderav(4)%CO(DIM1))
  powderav(4)%name="LEB3J023_004656"
  DO i=1,DIM1
     powderav(4)%CO(i)%a=LEB3J023_004656(i)%a 
     powderav(4)%CO(i)%b=LEB3J023_004656(i)%b 
     powderav(4)%CO(i)%g=LEB3J023_004656(i)%g 
     powderav(4)%CO(i)%w=LEB3J023_004656(i)%w 
  ENDDO

  DIM1=UBOUND(LEB3J029_009060,DIM=1)
  ALLOCATE(powderav(5)%CO(DIM1))
  powderav(5)%name="LEB3J029_009060"
  DO i=1,DIM1
     powderav(5)%CO(i)%a=LEB3J029_009060(i)%a 
     powderav(5)%CO(i)%b=LEB3J029_009060(i)%b 
     powderav(5)%CO(i)%g=LEB3J029_009060(i)%g 
     powderav(5)%CO(i)%w=LEB3J029_009060(i)%w 
  ENDDO

  DIM1=UBOUND(LEB3J035_015624,DIM=1)
  ALLOCATE(powderav(6)%CO(DIM1))
  powderav(6)%name="LEB3J035_015624"
  DO i=1,DIM1
     powderav(6)%CO(i)%a=LEB3J035_015624(i)%a 
     powderav(6)%CO(i)%b=LEB3J035_015624(i)%b 
     powderav(6)%CO(i)%g=LEB3J035_015624(i)%g 
     powderav(6)%CO(i)%w=LEB3J035_015624(i)%w 
  ENDDO

  DIM1=UBOUND(LEB3J041_024780,DIM=1)
  ALLOCATE(powderav(7)%CO(DIM1))
  powderav(7)%name="LEB3J041_024780"
  DO i=1,DIM1
     powderav(7)%CO(i)%a=LEB3J041_024780(i)%a 
     powderav(7)%CO(i)%b=LEB3J041_024780(i)%b 
     powderav(7)%CO(i)%g=LEB3J041_024780(i)%g 
     powderav(7)%CO(i)%w=LEB3J041_024780(i)%w 
  ENDDO


  !##################################################################
  !Read Control file
  Call ReadControlFile('control-sequence.dat',recoupledisotope,observedisotope,&
       &recoupledatomtype,observedatomtype,powderspec,&
       &convergence_criterion,DYNfile,dynamiccheck,XYZfile,filecheck)
  !##################################################################
  ! Read XYZ-file
  CALL ReadXYZfile(filecheck,XYZfile,recoupledatomtype,observedatomtype,&
       &dynamiccheck,DYNfile,obs,voids,recoup)
  !## uncomment the following to check ReadXYZfile
  write(*,*) "Dimension obs:",UBOUND(obs,Dim=1)
  write(*,*) "Dimension voids:",UBOUND(voids,Dim=1)
  write(*,*) "Dimension recoup:",UBOUND(recoup,Dim=1)
  Do i=1,UBOUND(obs,Dim=1)
     write(*,*) obs(i)%atomtype_xyz,obs(i)%x_xyz,obs(i)%weight 
  ENDDO
!  Do i=1,UBOUND(voids,Dim=1)
!     write(*,*) voids(i)%atomtype_xyz,voids(i)%x_xyz,voids(i)%weight 
!  ENDDO
  Do i=1,UBOUND(recoup,Dim=1)
     write(*,*) recoup(i)%atomtype_xyz !,obs(i)%x_xyz,weight 
     write(*,*) recoup(i)%x_xyz, recoup(i)%y_xyz,recoup(i)%z_xyz
     Do j=1,UBOUND(recoup(i)%AtomNr,Dim=1)
        write(*,*) recoup(i)%AtomNr(j),recoup(i)%weight(j)  
        IF(Associated(recoup(i)%voidNr)) THEN
           write(*,*) recoup(i)%voidNr(j),recoup(i)%void_weight(j)
        ENDIF
!        Do k=1, UBOUND(recoup(i)%AtomNr(j)%weight,Dim=1)
!           write(*,*) recoup(i)%AtomNr(j)%weight(k)
!        ENDDO
     ENDDO
  ENDDO
  !##################################################################

  !################ TIME SCALE ########################################## 
  CALL ReadTimeFile('time-sequence.dat',grid_array)

  !lies die versch. Zeiten ein und speichere sie in ein array, damit sie später mit eingerechnet werden können.


  !#######################Read Procedures Finished #########################

  !###############'Selection Part'##########################

  Do i=1, 23
     IF(recoupledisotope==gyro_array(i)%isotope) THEN
        gyro_ri = gyro_array(i)%gyrovalue
     END IF

     IF (observedisotope==gyro_array(i)%isotope) THEN
        gyro_oi = gyro_array(i)%gyrovalue
     END IF
  END DO


  !#################################################################
  ! Depending on the powderspec read from control-sequencefile 
  ! we now select our powder averages:
  SELECT CASE(powderspec)
  CASE ('auto')
     powderstart=1
     powderend=7
  CASE ('userspecified.cry')
     user=1
     OPEN (24, FILE='userspecified.cry', ACCESS='SEQUENTIAL', STATUS='OLD')
     ! COUNT LINES TO ALLOCATE THE POWDERARRAY
     DO WHILE(ios24==0)
        READ(24,FMT=*,IOSTAT=ios24) line
        IF(ios24==0) counter5=counter5+1
     END DO
     !Kein Powderarray mehr! siehe die 8.spalte in dem neuen array!
     CLOSE(24)
     powderstart=8
     powderend=8
  CASE ('LEB3J005_000084')
     powderstart=1
     powderend=1
  CASE ('LEB3J011_000600')
     powderstart=2
     powderend=2
  CASE ('LEB3J017_001980')
     powderstart=3
     powderend=3
  CASE ('LEB3J023_004656')
     powderstart=4
     powderend=4
  CASE ('LEB3J029_009060')
     powderstart=5
     powderend=5
  CASE ('LEB3J035_015624')
     powderstart=6
     powderend=6
  CASE ('LEB3J041_024780')
     powderstart=7
     powderend=7
  CASE DEFAULT
     OPEN(50,FILE='powder_error.dat',STATUS='replace',ACTION='write')
     WRITE(50,*)'The powderfile which you want to use does not exist!'
     WRITE(50,*)'You can either create a powder file with the ending .cry'
     WRITE(50,*)'or choose one from the following list.'
     WRITE(50,*)'  LEB3J005_000084'
     WRITE(50,*)'  LEB3J011_000600'
     WRITE(50,*)'  LEB3J017_001980'
     WRITE(50,*)'  LEB3J023_004656'
     WRITE(50,*)'  LEB3J029_009060'
     WRITE(50,*)'  LEB3J035_015624'
     WRITE(50,*)'  LEB3J041_024780'
     CLOSE(50)
     WRITE(*,*)''
     WRITE(*,*)'###########################################################################'
     WRITE(*,*)'You are trying to use a powderfile which is not implemented in this program.' 
     WRITE(*,*)'Please check the file powder_error.dat !'
     WRITE(*,*)'###########################################################################'
     WRITE(*,*)''
     STOP
  END SELECT


  IF (user==1) THEN
     OPEN (24, FILE='userspecified.cry', ACCESS='SEQUENTIAL', STATUS='OLD')
     ALLOCATE(powderav(8)%CO(counter5))
     DO i=1, counter5
        READ(24,FMT=*,IOSTAT=ios24)alp, bet, gam, wei
        powderav(8)%CO(i)%a=alp 
        powderav(8)%CO(i)%b=bet 
        powderav(8)%CO(i)%g=gam 
        powderav(8)%CO(i)%w=wei 
        WRITE(*,*) powderav(8)%CO(i)%a
        WRITE(*,*) powderav(8)%CO(i)%b
        WRITE(*,*) powderav(8)%CO(i)%g
        WRITE(*,*) powderav(8)%CO(i)%w
     END DO
     CLOSE(24)
  END IF


  !################ KERNEL 


  !Sk_In for different powder averages:
  DO powdervariable=powderstart , powderend
     WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*,*)''
     IF(user==0)THEN
        WRITE(*,*)'calculation with:', powderav(powdervariable)%name  
        WRITE(*,*)''
        WRITE(*,*)'#tau_uds/s    DeltaS/S0'  
     ELSE
        WRITE(*,*) 'calculation with: userspecified.cry'
     END IF
     DIM1=UBOUND(powderav(powdervariable)%CO,DIM=1)

     !tau-loop:
     DO t=1, UBOUND(grid_array,DIM=1)
        Signal_Sk_In=0.0d0    
        TwoSqrtTwoDwTOvThree=((2.0d0*sqrt(2.0d0)*grid_array(t))/3.0d0)
        ! loop over all S spins:
!dir$ loop count min(1)
        DO o=1, UBOUND(obs,DIM=1)
           Signal_pwdr_avrg = 0.0d0
           !powder average-loop:
           DO p=1, DIM1
              !loop over all I-spins:
              OneminTwoCosBeta_crTwo=1.0d0-2.0d0*cos(powderav(powdervariable)%CO(p)%b)*cos(powderav(powdervariable)%CO(p)%b)
              SinCosBeta_cr=sin(powderav(powdervariable)%CO(p)%b)*cos(powderav(powdervariable)%CO(p)%b)      
              SinGamma_cr=sin(powderav(powdervariable)%CO(p)%g)
              CosGamma_cr=cos(powderav(powdervariable)%CO(p)%g)
              CosBeta_cr=cos(powderav(powdervariable)%CO(p)%b)
              SinBeta_cr=sin(powderav(powdervariable)%CO(p)%b)
              Signal_S_In=1.0d0
              DO r=1,UBOUND(recoup,DIM=1)
                 redCosinusTerm=0
                 !     WRITE(*,*) "r=",r,"UBOUND(recoup(r)%AtomNr,DIM=1)",UBOUND(recoup(r)%AtomNr,DIM=1)

                 Do j=1,UBOUND(recoup(r)%AtomNr,DIM=1)
                    !     WRITE(*,*) "J=",J
                    !                    IF(r==9) then
                    !                    WRITE(*,*) recoup(r)%AtomNr(j)
                    !                    ENDIF
                    !                    WRITE(*,*) recoup(recoup(r)%AtomNr(j))%AtomNr(1),recoup(recoup(r)%AtomNr(j))%x_xyz, &
                    !                               recoup(recoup(r)%AtomNr(j))%y_xyz,&
                    !                               recoup(recoup(r)%AtomNr(j))%z_xyz
                    delta_x=recoup(recoup(r)%AtomNr(j))%x_xyz-obs(o)%x_xyz
                    delta_y=recoup(recoup(r)%AtomNr(j))%y_xyz-obs(o)%y_xyz
                    delta_z=recoup(recoup(r)%AtomNr(j))%z_xyz-obs(o)%z_xyz
                    vector_length= 1.0D-10*sqrt(delta_x**2.0d0 +delta_y**2.0d0 + delta_z**2.0d0)
                    D_k = -(gyro_oi*gyro_ri*h_quer*1.0d-7)/(2.0d0*PI*(vector_length**3.0d0))
                    beta_pc_rad = acos(delta_z*1.0d-10/vector_length)
                    gamma_pc_rad= modulo(atan2(delta_y,delta_x),2.0d0*pi)
                    Tau_k = powderav(powdervariable)%CO(p)%a + gamma_pc_rad

                    C_1k = 3.0d0*OneminTwoCosBeta_crTwo*sin(2.0d0*beta_pc_rad)*cos(Tau_k) &
                         + 3.0d0*(sin(beta_pc_rad)*sin(beta_pc_rad)*cos(2.0d0*Tau_k) &
                         - 3.0d0*cos(beta_pc_rad)*cos(beta_pc_rad)+1.0d0)*SinCosBeta_cr
                    S_1k = 3.0d0*(sin(2.0d0*beta_pc_rad)*sin(Tau_k)*CosBeta_cr-sin(beta_pc_rad) &
                         *sin(beta_pc_rad)*sin(2.0d0*Tau_k)*SinBeta_cr)
                    !                    redCosinusTerm=recoup(r)%weight(j)*redCosinusTerm+D_k*(C_1k*SinGamma_cr-S_1k*CosGamma_cr) !sum over k
                    redCosinusTerm=redCosinusTerm+recoup(r)%weight(j)*D_k*(C_1k*SinGamma_cr-S_1k*CosGamma_cr) !sum over k
                 ENDDO
                 IF(Associated(recoup(r)%VoidNr)) THEN
                    Do j=1,UBOUND(recoup(r)%VoidNr,DIM=1)
                       !     WRITE(*,*) "J=",J
                       !                    IF(r==9) then
                       !                    WRITE(*,*) recoup(r)%AtomNr(j)
                       !                    ENDIF
                       !                    WRITE(*,*) recoup(recoup(r)%AtomNr(j))%AtomNr(1),recoup(recoup(r)%AtomNr(j))%x_xyz, &
                       !                               recoup(recoup(r)%AtomNr(j))%y_xyz,&
                       !                               recoup(recoup(r)%AtomNr(j))%z_xyz
                       delta_x=voids(recoup(r)%VoidNr(j))%x_xyz-obs(o)%x_xyz
                       delta_y=voids(recoup(r)%VoidNr(j))%y_xyz-obs(o)%y_xyz
                       delta_z=voids(recoup(r)%VoidNr(j))%z_xyz-obs(o)%z_xyz
                       vector_length= 1.0D-10*sqrt(delta_x**2.0d0 +delta_y**2.0d0 + delta_z**2.0d0)
                       D_k = -(gyro_oi*gyro_ri*h_quer*1.0d-7)/(2.0d0*PI*(vector_length**3.0d0))
                       beta_pc_rad = acos(delta_z*1.0d-10/vector_length)
                       gamma_pc_rad= modulo(atan2(delta_y,delta_x),2.0d0*pi)
                       Tau_k = powderav(powdervariable)%CO(p)%a + gamma_pc_rad

                       C_1k = 3.0d0*OneminTwoCosBeta_crTwo*sin(2.0d0*beta_pc_rad)*cos(Tau_k) &
                            + 3.0d0*(sin(beta_pc_rad)*sin(beta_pc_rad)*cos(2.0d0*Tau_k) &
                            - 3.0d0*cos(beta_pc_rad)*cos(beta_pc_rad)+1.0d0)*SinCosBeta_cr
                       S_1k = 3.0d0*(sin(2.0d0*beta_pc_rad)*sin(Tau_k)*CosBeta_cr-sin(beta_pc_rad) &
                            *sin(beta_pc_rad)*sin(2.0d0*Tau_k)*SinBeta_cr)
                       !                    redCosinusTerm=recoup(r)%weight(j)*redCosinusTerm+D_k*(C_1k*SinGamma_cr-S_1k*CosGamma_cr) !sum over k
                       redCosinusTerm=redCosinusTerm+recoup(r)%void_weight(j)*D_k*(C_1k*SinGamma_cr-S_1k*CosGamma_cr) !sum over k
                    ENDDO
                 ENDIF
                 Signal_S_I= cos(TwoSqrtTwoDwTOvThree*redCosinusTerm) ! cos-term of eq (13)
                 Signal_S_In=Signal_S_In*Signal_S_I
              END DO

              Signal_pwdr_avrg = Signal_S_In*(powderav(powdervariable)%CO(p)%w) + Signal_pwdr_avrg
           END DO
           Signal_Sk_In= Signal_pwdr_avrg+Signal_Sk_In
        END DO
        Signal_Sk_In_tau(powdervariable,t)=1.0d0-Signal_Sk_In/DBLE(UBOUND(obs,DIM=1))
        WRITE(*,*) kappa*1.0D0*grid_array(t), Signal_Sk_In_tau(powdervariable,t)
     END DO

     !convergence criterion
     IF (powdervariable>1 .AND. powderstart<8) THEN
        WRITE(*,*)''
        WRITE(*,*) 'Checking for convergence'
        WRITE(*,*)''
        convergencesum=0
        DO t=1, UBOUND(grid_array,DIM=1)
           IF (powdervariable==1) THEN
              EXIT !no convergence possible if only 1 powder is calculated
           END IF
           convergencesum=convergencesum+((Signal_Sk_In_tau(powdervariable,t))-(Signal_Sk_In_tau(powdervariable-1,t)))**2
        END DO
        convergencepoint=sqrt(convergencesum/(UBOUND(grid_array,DIM=1)-1))
        WRITE(*,*) "Number of points: ",UBOUND(grid_array,DIM=1)
        WRITE(*,*) "Convergence criterion: ",convergence_criterion
        WRITE(*,*) "Squared difference to previous: ",convergencepoint
        WRITE(*,*)''
        IF (convergencepoint<convergence_criterion) THEN
           WRITE(*,*)''
           WRITE(*,*) 'converged with: ', powderav(powdervariable-1)%name
           WRITE(*,*)''
           EXIT !STOP
        END IF
     END IF
  END DO
  WRITE(*,*)'tau_uds: dephasing duration on the universal desphasing scale'
  WRITE(*,*)'         http://dx.doi.org/10.1021/cm800805f'
  WRITE(*,*)''
  WRITE(*,*)'Glamor reference: http://dx.doi.org/10.1016/j.ssnmr.2012.10.001'
  WRITE(*,*)''
  !deallocate all pointers
  do i=1,8
     IF(associated(powderav(i)%CO)) DEALLOCATE(powderav(i)%CO)
  enddo
  DEALLOCATE(powderav)
  DEALLOCATE(obs)
  do i=1,counter1
     IF(associated(recoup(i)%AtomNr)) DEALLOCATE(recoup(i)%AtomNr)
     IF(associated(recoup(i)%VoidNr)) DEALLOCATE(recoup(i)%VoidNr)
     IF(associated(recoup(i)%weight)) DEALLOCATE(recoup(i)%weight)
     IF(associated(recoup(i)%void_weight)) DEALLOCATE(recoup(i)%void_weight)
  enddo
  deallocate(recoup)
  IF(associated(voids)) deallocate(voids)
  IF(associated(tmpAtomNrlist)) deallocate(tmpAtomNrlist)
  IF(associated(tmpWeightlist)) deallocate(tmpWeightlist)
  deallocate(grid_array)

END PROGRAM glamoR
