MODULE ReadTimeFile_MOD
  private
  public :: ReadTimeFile
contains
  SUBROUTINE ReadTimeFile(timefilename,grid_array)
    IMPLICIT NONE
    !Variable Declarations
    CHARACTER(len=*):: timefilename ! name of inputfile to open
    DOUBLE PRECISION::startvalue, increment
    INTEGER:: gridnumber,ios25
    INTEGER:: z, n, counter6
    CHARACTER(len=30):: a

    DOUBLE PRECISION, DIMENSION(:), POINTER::grid_array


    OPEN(25, FILE='time-sequence.dat', ACCESS='SEQUENTIAL', STATUS='OLD')
    READ(25,FMT=*,IOSTAT=ios25) a
    IF (a=='%timetype') THEN
       READ(25,FMT=*,IOSTAT=ios25) a
    END IF

    SELECT CASE(a)
    CASE ('linear')
       READ(25,FMT=*,IOSTAT=ios25) a
       READ(25,FMT=*,IOSTAT=ios25)startvalue, increment, gridnumber
       ALLOCATE(grid_array(gridnumber))
       DO z=0, gridnumber-1
          startvalue=z*increment
          grid_array(z+1)=startvalue
       END DO
    CASE ('grid')        
       BACKSPACE(25, IOSTAT=ios25)
       READ(25,FMT=*,IOSTAT=ios25)a, gridnumber
       READ(25,FMT=*,IOSTAT=ios25)a
       ALLOCATE(grid_array(gridnumber))
       DO n=1, gridnumber
          counter6=counter6+1
          READ(25,FMT=*,IOSTAT=ios25) grid_array(counter6)  
       END DO
    CASE DEFAULT
       WRITE(*,*)''
       WRITE(*,*)'########################################################################################################'
       WRITE(*,*)'You have selected an undefined timetype. Please select *linear* or *grid* in the file time_sequence.dat.'
       WRITE(*,*)'#########################################################################################################'
       WRITE(*,*)''
       STOP      
    END SELECT
    CLOSE(25)

  END SUBROUTINE ReadTimeFile
END MODULE ReadTimeFile_MOD
