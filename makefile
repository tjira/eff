INCLUDE := -I include
FLAGS := -O3
LIBS := -lm

OBJECTS = eff_access.o eff_create.o eff_driver.o eff_dynamics.o eff_erf.o eff_global.o eff_update.o eff_util.o eff_shanno.o eff_minimize.o eff_makemol.o eff_output.o eff_initial.o eff_input.o eff_params.o eff_bounds.o eff_ewald.o eff_timing.o eff_pressure.o eff_efield.o eff_collapse.o eff_properties.o eff_cores.o eff_eigen.o

all: bin build bin/eff

bin/eff: $(addprefix build/, $(OBJECTS))
	gcc $(FLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

$(addprefix build/, $(OBJECTS)): build/%.o: src/%.c
	gcc -c $(FLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

bin build:
	@mkdir -p $@

clean:
	rm -rf bin build
