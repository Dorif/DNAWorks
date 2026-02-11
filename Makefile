# compiler
FC = gfortran

# compile flags (release by default; use `make debug` for debug build)
FCFLAGS = -Ofast -march=native -flto -funroll-loops

# program name
PROGRAM = dnaworks

# required objects
objects = dnaworks.o dnaworks_data.o dnaworks_test.o \
	control_func.o encoding.o input.o misc_func.o \
	mutate.o output.o overlaps.o scores.o str_func.o time_func.o

# required modules
modules = dnaworks_data.mod dnaworks_test.mod

# the main linking step
$(PROGRAM): $(objects)
	$(FC) $(FCFLAGS) -o $(PROGRAM) $(objects)

# specific requirements for each object
$(objects): $(modules)

# compile recipe for modules
%.mod: %.f90
	$(FC) $(FCFLAGS) -c $<

# compile recipe for objects
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# debug build with bounds checking and warnings
.PHONY: debug
debug: FCFLAGS = -g -O0 -fbounds-check -Wall -Wextra -fbacktrace -static-libgcc
debug: clean $(PROGRAM)

# extra rules
.PHONY: clean
clean:
	rm -f $(objects) $(modules) $(PROGRAM)
