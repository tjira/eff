# Makefile for eff2 -- Julius Su

CXX = icc
CXXFLAGS = -O2 -xW -axW -iop -static
CC = icc
CFLAGS = -O2 -xW -axW -iop -static

LIBS = -lm

OBJECTS = eff_access.o eff_create.o eff_driver.o eff_dynamics.o eff_erf.o eff_global.o eff_update.o eff_util.o eff_shanno.o eff_minimize.o eff_makemol.o eff_output.o eff_initial.o eff_input.o eff_params.o eff_bounds.o eff_ewald.o eff_timing.o eff_pressure.o eff_efield.o eff_collapse.o eff_properties.o eff_cores.o eff_eigen.o

OBJECTS2 = makemol.o

eff: $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

eff_access.o: eff_access.h eff_global.h
eff_create.o: eff_create.h eff_global.h
eff_driver.o: eff_driver.h
eff_dynamics.o: eff_dynamics.h
eff_erf.o: eff_erf.h
eff_global.o: eff_global.h
eff_update.o: eff_update.h
eff_util.o: eff_util.h
eff_shanno.o: eff_shanno.h
eff_minimize.o: eff_minimize.h
eff_makemol.o: eff_makemol.h
eff_output.o: eff_output.h
eff_initial.o: eff_initial.h
eff_input.o: eff_input.h
eff_params.o: eff_params.h
eff_bounds.o: eff_bounds.h
eff_ewald.o: eff_ewald.h
eff_timing.o: eff_timing.h
eff_efield.o: eff_efield.h
eff_collapse.o: eff_collapse.h
eff_properties.o: eff_properties.h
eff_cores.o: eff_cores.h
eff_eigen.o: eff_eigen.h
