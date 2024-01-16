MODULE POWDER_TYPE_DEC
IMPLICIT NONE
PRIVATE
PUBLIC::powder_type, CO_type !, powderav

TYPE CO_type
  DOUBLE PRECISION :: a, b, g, w
END TYPE CO_type

TYPE powder_type
  CHARACTER(len=20) :: name
  TYPE (CO_type), DIMENSION(:), POINTER :: CO 
END TYPE powder_type
! TYPE(powder_type), DIMENSION(7) :: powderav

CONTAINS
END MODULE POWDER_TYPE_DEC






