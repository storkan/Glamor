MODULE declarations
  IMPLICIT NONE
  TYPE :: xyz_static 
    CHARACTER( len =3):: atomtype_xyz            ! ch. symbol of the element
    DOUBLE PRECISION:: x_xyz, y_xyz, z_xyz       ! coorinates of the static xyz file
    DOUBLE PRECISION:: weight                    ! weight of the occupied position
  END TYPE xyz_static

  TYPE :: xyz_dynamic 
    CHARACTER( len =3):: atomtype_xyz            ! ch. symbol of the element
    DOUBLE PRECISION :: x_xyz, y_xyz, z_xyz      ! coorinates of the static xyz file
    INTEGER, DIMENSION(:),Pointer :: AtomNr
    INTEGER, DIMENSION(:),Pointer :: VoidNr
    DOUBLE PRECISION,DIMENSION(:),POINTER :: weight               ! weight of the occupied position
    DOUBLE PRECISION,DIMENSION(:),POINTER :: void_weight          ! weight of the void position
  END TYPE xyz_dynamic
END MODULE declarations

