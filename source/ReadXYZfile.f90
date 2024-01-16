MODULE ReadXYZfile_MOD
  private
  public :: ReadXYZfile
contains
  SUBROUTINE ReadXYZfile(filecheck,XYZfile,recoupledatomtype,observedatomtype,dynamiccheck,DYNfile,obs,voids,recoup)
    USE declarations
    USE shift_mod
    !Declare Variables
    IMPLICIT NONE
    INTEGER :: filecheck,ios22=0,ios26=0,counter1=0,counter2=0,counter3=0,counter4
    INTEGER :: dynamiccheck,linecounter,randomnumber
    INTEGER :: voidstart,voidend,voidcounter,NrPrimaryLines
    INTEGER :: tempAtomNr,OtherPositions,OtherVoids
    INTEGER :: i,j,k,index,v_index,Pos,voidPos,tmpDIM
    CHARACTER(len=30):: XYZfile, DYNfile,a,line
    CHARACTER(len=5)::  observedatomtype, recoupledatomtype
    CHARACTER(len=3)::tempAtomName
    DOUBLE PRECISION :: weight,tempX,tempY,tempZ,voidWeight
    TYPE(xyz_static), DIMENSION(:), POINTER:: obs
    TYPE(xyz_dynamic), DIMENSION(:), POINTER:: recoup,voids
    INTEGER, DIMENSION(:), POINTER :: tmpAtomNrlist
    DOUBLE PRECISION, DIMENSION(:), POINTER :: tmpWeightlist
    ! Read XYZ-file
    ! check if input files (which one?...) are correct
    counter4=0
    IF(filecheck==0) THEN
       ! CHECK FOR OBS AND REC
       OPEN(22, FILE=XYZfile, ACCESS='SEQUENTIAL', STATUS='OLD')
       ios22=0
       READ (22,FMT=*,IOSTAT=ios22) a
       READ (22,FMT=*,IOSTAT=ios22) a
       DO WHILE(ios22==0)    
          READ (22,FMT=*,IOSTAT=ios22) a
          ! Determine COUNTERs FOR ALLOCATION OF ARRAYS
          IF (a==recoupledatomtype .OR. a==observedatomtype) THEN
             IF (a==recoupledatomtype .AND. ios22==0) THEN
                counter1=counter1+1 
             END IF
             IF (a==observedatomtype .AND. ios22==0) THEN
                counter2=counter2+1
             END IF
             CONTINUE
          ELSE
             WRITE(*,*)''
             WRITE(*,*)'########################################################################'
             WRITE(*,*) 'There is something wrong with your atomtypes, please check your file!'
             WRITE(*,*)'########################################################################'
             WRITE(*,*)''
             STOP
          END IF
       END DO
       REWIND(22)

       ! ALLOCATION OF THE 2 ARRAYS
       ALLOCATE(obs(counter2))
       ALLOCATE(recoup(counter1))

       ! Determination of the number of voids in order to allocate the array (recoup)
       IF(dynamiccheck==0) THEN
           DO i=1,UBOUND(recoup,DIM=1)
!             IF (.NOT.ASSOCIATED(recoup(i)%AtomNr)) THEN
                ALLOCATE(recoup(i)%AtomNr(1))
                ALLOCATE(recoup(i)%weight(1))
                recoup(i)%AtomNr(1)=i
                recoup(i)%weight(1)=1.0D0
!             ENDIF
          ENDDO
       ELSE !  dynamiccheck==1, i.e. we calculate with nuclear dynamics
          OPEN(26, FILE=DYNfile, ACCESS='SEQUENTIAL', STATUS='OLD')
          ios26=0
          linecounter=0
          DO WHILE(ios26==0)    
             READ (26,FMT=*,IOSTAT=ios26) a
             linecounter=linecounter+1
             ! Determine the number of voids through a counter        
             IF (a=='%voids') THEN
                voidstart=linecounter+1
             ELSEIF(a=='%end_voids') THEN
                voidend=linecounter
                voidcounter=voidend-voidstart
                ALLOCATE(voids(voidcounter))
                write (*,*) voidcounter, voidend, voidstart
             ELSEIF(a=='%dynamic') THEN
                READ(26,FMT=*,IOSTAT=ios26) NrPrimaryLines
                DO i=1,NrPrimaryLines
                   READ(26,FMT=*,IOSTAT=ios26) tempAtomNr,OtherPositions,OtherVoids,weight

                   ALLOCATE(recoup(tempAtomNr)%AtomNr(1+OtherPositions)) ! allocate memory for AtomNr-list of primary atom
                   ALLOCATE(recoup(tempAtomNr)%weight(1+OtherPositions))
                   ALLOCATE(recoup(tempAtomNr)%VoidNr(OtherVoids))       ! allocate memory for VoidNr-list of primary atom
                   ALLOCATE(recoup(tempAtomNr)%void_weight(OtherVoids))
                   recoup(tempAtomNr)%AtomNr(1)=tempAtomNr               ! set first item in AtomNr-list of primary atom to the own atoms number
                   recoup(tempAtomNr)%weight(1)=weight
                   index=1
                   v_index=0
                   Do j=1,OtherPositions
                      READ(26,FMT=*,IOSTAT=ios26) Pos,weight
                      index=index+1
                      recoup(tempAtomNr)%AtomNr(index)=Pos     !append secondary atom to AtomNr-list of primary atom
                      recoup(tempAtomNr)%weight(index)=weight  !append secondary atom weight to weight-list of primary atom
                      ALLOCATE(recoup(Pos)%AtomNr(1+OtherPositions)) ! allocate memory for AtomNr-list of secondary atom
                      ALLOCATE(recoup(Pos)%weight(1+OtherPositions)) ! allocate memory for weight-list of secondary atom
                      ALLOCATE(recoup(Pos)%VoidNr(OtherPositions))     ! allocate memory for VoidNr-list of secondary atom
                      ALLOCATE(recoup(Pos)%void_weight(OtherPositions))! allocate memory for void_weight-list of secondary atom
                   ENDDO

                   Do j=1,OtherVoids
                      READ(26,FMT=*,IOSTAT=ios26) voidPos,voidWeight
                      v_index=v_index+1
                      recoup(tempAtomNr)%VoidNr(v_index)=voidPos           ! append voidNr to VoidNr-list of primary atom
                      recoup(tempAtomNr)%void_weight(v_index)=voidWeight   ! append voidWeight to void_weight-list of primary atom
                   ENDDO
                   !all secondary atoms and voids connected with tempAtomNr have been read in at this point
                   !now update fields of the secondary atoms
                   IF (ubound(recoup(tempAtomNr)%AtomNr,DIM=1)>1) THEN !set the AtomNr-list entries for the *secondary* atoms
                      tmpDIM=ubound(recoup(tempAtomNr)%AtomNr,DIM=1)   ! tmpDIM = No of entries in AtomNr-list of primary atom 
                      Allocate(tmpAtomNrlist(tmpDIM))                  
                      Allocate(tmpWeightlist(tmpDIM))                  
!                      tmpDIM2=ubound(recoup(tempAtomNr)%VoidNr,DIM=1)   ! tmpDIM = No of entries in AtomNr-list of primary atom 
!                      Allocate(tmpVoidNrlist(tmpDIM2))                  
!                      Allocate(tmpVoid_weightlist(tmpDIM2))                  
                      tmpAtomNrlist=recoup(tempAtomNr)%AtomNr  ! vector aqssignment, copies AtomNr-list of primary atom in a temporary List
                      tmpWeightlist=recoup(tempAtomNr)%weight  ! vector aqssignment, copies Weight-list of primary atom in a temporary List
                      Do k=2,ubound(recoup(tempAtomNr)%AtomNr,DIM=1)
                         CALL shift(tmpAtomNrlist,tmpWeightlist)  !shift subroutine shifts values in tmpAtomNrlist to the left and 
                         !puts the element at index one to the end
                         recoup(tmpAtomNrlist(k))%AtomNr=tmpAtomNrlist
                         recoup(tmpAtomNrlist(k))%weight=tmpWeightlist !JW: replaced %AtomNr by %weight
                         recoup(tmpAtomNrlist(k))%VoidNr=recoup(tempAtomNr)%VoidNr !vector aqssignment, copies VoidNr-list of primary atom 
                                                                                   !to VoidNr-list of secondary atom
                         recoup(tmpAtomNrlist(k))%Void_weight=recoup(tempAtomNr)%Void_weight !vector aqssignment, copies Void_weight-list of primary atom 
                                                                                   !to Void_weight-list of secondary atom
                         WRITE(*,*) k, tmpAtomNrlist, tmpWeightlist
                      ENDDO
                      
                   ENDIF
                ENDDO
             ENDIF
          END DO

          REWIND(26)
          ios26=0

          DO WHILE(ios26==0) 
             READ (26,FMT=*,IOSTAT=ios26) a
             IF (a=='%voids') THEN  
                DO i=1, voidcounter 
                   READ (26,FMT=*,IOSTAT=ios26) randomnumber,tempX,tempY,tempZ
write (*,*) 'blaaaa'
                   voids(i)%atomtype_xyz='BQ'
                   voids(i)%x_xyz=tempX
                   voids(i)%y_xyz=tempY
write (*,*) 'blaaaa1'
                   voids(i)%z_xyz=tempZ
write (*,*) 'blaaaa2'
!                   voids(i)%weight=0.1D0
write (*,*) 'blaaaaaaaaaa'
                END DO
             END IF
          END DO


          CLOSE(26)
          DO i=1,UBOUND(recoup,DIM=1)
             IF (.NOT.ASSOCIATED(recoup(i)%AtomNr)) THEN
                ALLOCATE(recoup(i)%AtomNr(1))
                ALLOCATE(recoup(i)%weight(1))
                recoup(i)%AtomNr(1)=i
                recoup(i)%weight(1)=1.0D0
             ENDIF
          ENDDO
       END IF

       ! write data about the atoms into the arrays recoup and obs
       READ(22,FMT=*,IOSTAT=ios22) line
       READ(22,FMT=*,IOSTAT=ios22) line
       DO WHILE(ios22==0)
          READ(22,FMT=*,IOSTAT=ios22) tempAtomName,tempX,tempY,tempZ
          IF (tempAtomName==recoupledatomtype .AND. ios22==0) THEN
             counter3=counter3+1     
             recoup(counter3)%atomtype_xyz=tempAtomName
             recoup(counter3)%x_xyz=tempX
             recoup(counter3)%y_xyz=tempY
             recoup(counter3)%z_xyz=tempZ   
             !WRITE(*,*)counter3        
             !write(*,FMT='(a f10.5 f10.5 f10.5)')  recoup(counter3)%atomtype_xyz,&
             !                        recoup(counter3)%x_xyz,&
             !                        recoup(counter3)%y_xyz,&
             !                        recoup(counter3)%z_xyz
          ELSEIF(tempAtomName==observedatomtype .AND. ios22==0) THEN
             counter4=counter4+1
             obs(counter4)%atomtype_xyz=tempAtomName
             obs(counter4)%x_xyz=tempX
             obs(counter4)%y_xyz=tempY
             obs(counter4)%z_xyz=tempZ
             obs(counter4)%weight=1.0D0
             !write(*,FMT='(a f10.5 f10.5 f10.5)') obs(counter4)%atomtype_xyz,&
             !                         obs(counter4)%x_xyz,&
             !                         obs(counter4)%y_xyz,&
             !                         obs(counter4)%z_xyz
          END IF
       END DO
       CLOSE(22)
    ENDIF
    !STOP
    IF(filecheck==1) THEN
       WRITE(*,*)'I havent been working on cif-files yet. Please be patient.'
       STOP
    END IF

    IF(filecheck==2) THEN
       WRITE(*,*)''
       WRITE(*,*)'################################################################################################################'
       WRITE(*,*)'You are trying to read an inputfiletype that has not been intended for this program. Please choose -cif or -xyz.'
       WRITE(*,*)'################################################################################################################'
       WRITE(*,*)''
       STOP
    END IF
  END SUBROUTINE ReadXYZfile
END MODULE ReadXYZfile_MOD
