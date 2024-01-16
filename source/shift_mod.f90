module shift_mod
contains
SUBROUTINE shift(tmpAtomNrlist,tmpWeightlist)   !shift subroutine shifts values in tmpAtomNrlist and tmpWeightlist to the left 
INTEGER, DIMENSION(:),POINTER :: tmpAtomNrlist  !and puts the element at index one to the end
INTEGER, DIMENSION(:),POINTER :: int_Vector
INTEGER :: listDim
DOUBLE PRECISION, DIMENSION(:),POINTER :: tmpWeightlist
DOUBLE PRECISION, DIMENSION(:),POINTER :: double_Vector
listDim=ubound(tmpAtomNrlist,DIM=1)
ALLOCATE(int_Vector(listDim))
int_Vector(1:listDim-1)=tmpAtomNrlist(2:listDim)
int_Vector(listDim)=tmpAtomNrlist(1)
tmpAtomNrlist=int_Vector

ALLOCATE(double_Vector(listDim))
double_vector(1:listDim-1)=tmpWeightlist(2:listDim)
double_vector(listDim)=tmpWeightlist(1)
tmpWeightlist=double_Vector

DEALLOCATE(int_Vector)
DEALLOCATE(double_Vector)
END SUBROUTINE shift
end module shift_mod
