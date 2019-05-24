#
# This is a simple makefile for building spin-squared calculation code.
#
LIBS         = -llapack -lblas
#RunF         = pgfortran -i8 -r8 -Mallocatable=03
RunF         = pgfortran -i8 -r8
#
#
# The 'all' rule.
#
all: ci_singles_alpha_beta_strings.exe

#
# Generic rules for building module (*.mod) and object (*.o) files.
#
#%.mod: %.f90
#	$(RunF) -c $*.f90

%.o: %.f90
	$(RunF) -I$(MODS) -c $*.f90

#
# Generic rule for building general executable program (*.exe) from a standard
# f90 source (*.f90) file.
#
%.exe: %.f90
	$(RunF) $(LIBS) -o $*.exe $*.f90
