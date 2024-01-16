MODULE ReadControlFile_MOD
private
public :: ReadControlFile
contains
SUBROUTINE ReadControlFile(filename,recoupledisotope,observedisotope,&
                    &recoupledatomtype,observedatomtype,powderspec,&
                    &convergence_criterion,DYNfile,dynamiccheck,XYZfile,filecheck)
!Variable Declarations
IMPLICIT NONE
CHARACTER(len=*):: filename ! name of inputfile to open
CHARACTER(len=30):: a, b, b2 ! trashvariable only useful to read data into other variables
CHARACTER(len=5):: recoupledisotope, observedisotope, observedatomtype, recoupledatomtype
CHARACTER(len=30):: powderspec ! specifies what kind of powdraveraging shall be used
CHARACTER(len=30):: XYZfile, DYNfile
DOUBLE PRECISION:: convergence_criterion
INTEGER:: filecheck, dynamiccheck  !logical parameters
INTEGER:: ios23=0

  ! Read control-sequencefile 
  ! -read set recoupled and observed atoms
  ! -read the information about the powderaverages
  dynamiccheck=0 ! per defualt we calculate statically
  OPEN(23, FILE=filename, ACCESS='SEQUENTIAL', STATUS='OLD')
  ios23=0
  DO WHILE(ios23==0)
     READ(23,FMT=*,IOSTAT=ios23) a
     IF (a=='staticfile') THEN ! read recoupled and observed atoms
        READ (23,FMT=*,IOSTAT=ios23) a, recoupledisotope
        READ (23,FMT=*,IOSTAT=ios23) a, observedisotope
        READ (23,FMT=*,IOSTAT=ios23) a, recoupledatomtype
        READ (23,FMT=*,IOSTAT=ios23) a, observedatomtype
     ELSEIF (a=='%powderaverage') THEN ! read powderaverage specification
        READ (23,FMT=*,IOSTAT=ios23) a, powderspec
        READ (23,FMT=*,IOSTAT=ios23) a, convergence_criterion
     ELSEIF (a=='%dynamicsfile') THEN ! locate if the dynamics file is needed
        READ (23,FMT=*,IOSTAT=ios23) a, DYNfile
        dynamiccheck=1
     END IF
  END DO
  REWIND(23)
  ! locate which files are needed
  ios23=0
  DO WHILE(ios23==0)
     READ(23,FMT=*,IOSTAT=ios23) a
     IF (a=='%structureinput') THEN
        READ (23,FMT=*,IOSTAT=ios23) a, b
        IF(a=='inputtype' .AND. b=='xyz') THEN
           filecheck=0
           READ (23,FMT=*,IOSTAT=ios23) a, XYZfile
        ELSE
           BACKSPACE(23, IOSTAT=ios23)
           READ (23,FMT=*,IOSTAT=ios23) a, b2
           IF(a=='inputtype' .AND. b=='cif') THEN
              filecheck=1
           ELSE
              filecheck=2
           END IF
           EXIT
        END IF
     END IF
  END DO
  CLOSE(23)
END SUBROUTINE ReadControlFile
END MODULE ReadControlFile_MOD
