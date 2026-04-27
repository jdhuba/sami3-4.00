OBJ	= \
	parameter_mod.o\
	parameter_apex_mod.o\
	namelist_mod.o\
	variable_mod.o\
	message_passing_mod.o\
	chemistry_mod.o\
	time_mod.o\
	atomic_mod.o\
	conductance_mod.o\
	exb_mod.o\
	hardy_mod.o\
	misc_mod.o\
	photo_production_mod.o\
	grid_mod.o\
	msis2.o\
	gtd8d.o\
	grid3-4.00.o\
	madala_sevp_dp.o\
	sami3-4.00.o\
	sevp13_dp.o\
	hwm14.o\
	heating.o\
	solvers.o\
	photo_production.o\
	output.o\
	chemistry.o\
	potential.o\
	exb_transport.o\
	misc.o\
	neutral.o\
	hardy_sami3.o\
	weimer_omni.o\
	weimer_grid.o\
	weimer_phi.o\
	apexsh.o\
        magfld.o


# intel old ifort

  f90 = mpif90 -O2 -save -mcmodel=large -shared-intel 

# debug

#  f90 = mpif90 -C -traceback -shared-intel -mcmodel=large

.SUFFIXES : .o .f90 .f .F90 .F


%.o:%.mod

.f.o:
	$(f90) -c $*.f 

.f90.o:
	$(f90) -c $*.f90

.F90.o:
	$(f90) -c $*.F90

.F.o:
	$(f90) -c $*.F

clean:
	rm *.o *.x *.mod

sami3.x:	$(OBJ)
	$(f90) -o sami3.x $(OBJ)

#$(OBJ): com3_mpi-2.00.inc
#$(OBJ): param3_mpi-2.00.inc


