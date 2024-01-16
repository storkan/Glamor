#############################
#for gnu fortran:
F95 = gfortran
FFLAGS = -O3 -g -Wall

############################
#for intel fortran:
#F95 = ifort
#FFLAGS = -parallel -par-threshold0


#------------------------------------------------------

ALL: glamor
allways_younger:

shift_mod.o: shift_mod.f90
	$(F95) $(FFLAGS) -c shift_mod.f90

gyro_mag_mod.o: gyro_mag_mod.f90
	$(F95) $(FFLAGS) -c gyro_mag_mod.f90

powder_type_dec.o: powder_type_dec.f90
	$(F95) $(FFLAGS) -c powder_type_dec.f90

LEB3J005_000084.o: LEB3J005_000084.f90
	$(F95) $(FFLAGS) -c LEB3J005_000084.f90

LEB3J011_000600.o: LEB3J011_000600.f90
	$(F95) $(FFLAGS) -c LEB3J011_000600.f90

LEB3J017_001980.o: LEB3J017_001980.f90
	$(F95) $(FFLAGS) -c LEB3J017_001980.f90

LEB3J023_004656.o: LEB3J023_004656.f90
	$(F95) $(FFLAGS) -c LEB3J023_004656.f90

LEB3J029_009060.o: LEB3J029_009060.f90
	$(F95) $(FFLAGS) -c LEB3J029_009060.f90

LEB3J035_015624.o: LEB3J035_015624.f90
	$(F95) $(FFLAGS) -c LEB3J035_015624.f90

LEB3J041_024780.o: LEB3J041_024780.f90
	$(F95) $(FFLAGS) -c LEB3J041_024780.f90

declarations.o: declarations.f90
	$(F95) $(FFLAGS) -c declarations.f90

ReadControlFile.o: ReadControlFile.f90
	$(F95) $(FFLAGS) -c ReadControlFile.f90

ReadXYZfile.o: ReadXYZfile.f90
	$(F95) $(FFLAGS) -c ReadXYZfile.f90

ReadTimeFile.o: ReadTimeFile.f90
	$(F95) $(FFLAGS) -c ReadTimeFile.f90

objects = gyro_mag_mod.o powder_type_dec.o LEB3J005_000084.o LEB3J011_000600.o LEB3J017_001980.o \
          LEB3J023_004656.o LEB3J029_009060.o LEB3J035_015624.o LEB3J041_024780.o declarations.o shift_mod.o\
          ReadControlFile.o ReadXYZfile.o ReadTimeFile.o


glamor: glamor.f90 $(objects)
	$(F95) $(FFLAGS) -c glamor.f90
	$(F95) $(FFLAGS) -o glamor glamor.o $(objects)

#----------------------------------------------------
clean:
	$(RM) *.o *~ *.mod glamor
