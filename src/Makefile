# Set build type

#BUILD = debug
#BUILD = default
BUILD = production
#BUILD = openmp

# Compiler options for different machines

HOST = $(shell hostname)
Exe = ../mdsbi

# Mac Pro

ifeq ($(findstring civet,$(HOST)),civet)
 F95 = openmpif90 -fdefault-real-8 -fdefault-double-8
 ifeq ($(BUILD),debug)
  F95FLAGS = -g -Wall -Wextra -Wconversion -fbounds-check -fbacktrace \
	-fimplicit-none -std=f2003
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -O5 -Wuninitialized
 endif
 FFTW2 = -L/opt/local/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
 FFTW3 = -L/opt/local/lib -lfftw3
 LIBS0 = -framework vecLib
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw2.f90),)
	LIBS = $(FFTW2) $(LIBS0)
 endif
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw2c.f90),)
	LIBS = $(FFTW2) $(LIBS0)
 endif
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw3.f90),)
	LIBS = $(FFTW3) $(LIBS0)
 endif
endif

# Stanford CEES

ifeq ($(findstring cees,$(HOST)),cees)
 Exe = /data/dunham1/edunham/mlaunch/mdsbi
 F95 = mpif90 -r8 -i4 -i_dynamic
 ifeq ($(BUILD),debug)
  F95FLAGS = -g -check all -check noarg_temp_created -warn all \
	-traceback -ftrapuv -fpe0 -fp-stack-check -fltconsistency -std03
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -O2
 endif
 ifeq ($(BUILD),openmp)
  F95FLAGS = -O2 -openmp -parallel
 endif
 FFTW2 = 
 MKLPATH = /usr/local/INTEL/Compiler/11.1/069/mkl
 FFTW3 = -L$(MKLPATH)/include/fftw -lfftw3
 LIBS0 = -L$(MKLPATH) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential \
         -lmkl_core -Wl,--end-group -lpthread
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw2.f90),)
	LIBS = $(FFTW2) $(LIBS0)
 endif
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw2c.f90),)
	LIBS = $(FFTW2) $(LIBS0)
 endif
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw3.f90),)
	LIBS = $(FFTW3) $(LIBS0)
 endif
endif

# BlueGene/L
ifeq ($(ARCH),ppc64)
 make-on-BGL:
	@echo 'Running make -f Makefile.bgl (makefile for BlueGene/L)' 1>&2
	make -f Makefile.bgl
endif

ifneq ($(shell test -f fft_routines.f90 && echo yes),yes)
 abort-on-fft_routines:
	@echo 'Link appropriate fft_routines.f90 file:' 1>&2
	@echo 'ln -s fft_routines_fftw2.f90 fft_routines.f90 (for FFTW-2.x MPI, real)' 1>&2
	@echo 'ln -s fft_routines_fftw2c.f90 fft_routines.f90 (for FFTW-2.x MPI, complex)' 1>&2
	@echo 'ln -s fft_routines_fftw3.f90 fft_routines.f90 (for FFTW-3.x non-MPI)' 1>&2
	@exit 1
endif

# Files

Files = asperity.f90 constants.f90 convolution.f90 convolution_bm.f90 \
	convolution_im.f90 convolution_routines.f90 fault.f90 \
	fft_routines.f90 fields.f90 fourier.f90 friction.f90 \
	friction_bm.f90 friction_im.f90 friction_routines.f90 \
	friction_static.f90 front.f90 history.f90 init.f90 integration.f90 \
	io.f90 kernel.f90 load.f90 main.f90 mesh.f90 model.f90 \
	mpi_routines.f90 preslip.f90 problem.f90 problem_routines.f90 \
	rates.f90 ratestate.f90 rk.f90 slipweak.f90 static.f90 substep.f90 \
	test_utilities.f90 thermpres.f90 time_step.f90 utilities.f90

Obs = 	$(Files:.f90=.o)

$(Exe): $(Obs)
	$(F95) $(F95FLAGS) -o $(Exe) \
	$(Obs) $(LIBS) $(INCL)

%.o : %.f90
	$(F95) $(F95FLAGS) -c $< -o $@ $(INCL)

clean:
	rm -f *.o *.mod $(Exe)

.SUFFIXES: .o .f90

# Dependencies
# DO NOT DELETE THIS LINE - used by make depend
asperity.o: constants.o io.o


convolution.o: constants.o history.o io.o kernel.o model.o mpi_routines.o

convolution_bm.o: constants.o convolution.o history.o kernel.o model.o
convolution_bm.o: utilities.o

convolution_im.o: constants.o convolution.o history.o kernel.o model.o
convolution_im.o: utilities.o

convolution_routines.o: convolution.o convolution_bm.o convolution_im.o
convolution_routines.o: kernel.o model.o

fault.o: asperity.o constants.o io.o model.o mpi_routines.o preslip.o

fft_routines.o: constants.o convolution.o io.o model.o mpi_routines.o

fft_routines_fftw2.o: constants.o convolution.o io.o model.o mpi_routines.o

fft_routines_fftw2c.o: constants.o convolution.o io.o model.o mpi_routines.o

fft_routines_fftw3.o: constants.o convolution.o io.o model.o mpi_routines.o

fields.o: constants.o fault.o io.o load.o model.o thermpres.o

fourier.o: constants.o convolution.o fault.o fft_routines_fftw3.o io.o kernel.o
fourier.o: model.o

friction.o: constants.o fault.o io.o model.o mpi_routines.o ratestate.o
friction.o: slipweak.o

friction_bm.o: constants.o friction.o utilities.o

friction_im.o: constants.o friction.o utilities.o

friction_routines.o: constants.o convolution.o fault.o fourier.o friction.o
friction_routines.o: friction_bm.o friction_im.o friction_static.o io.o model.o
friction_routines.o: ratestate.o slipweak.o utilities.o

friction_solver.o: constants.o friction.o model.o utilities.o

friction_static.o: constants.o

front.o: constants.o fft_routines_fftw3.o fields.o io.o model.o mpi_routines.o

history.o: constants.o model.o

init.o: constants.o convolution.o fft_routines_fftw3.o fields.o friction.o
init.o: front.o io.o kernel.o mesh.o model.o mpi_routines.o problem.o

integration.o: constants.o fault.o io.o model.o mpi_routines.o

io.o: constants.o mpi_routines.o

kernel.o: constants.o io.o model.o mpi_routines.o utilities.o

load.o: asperity.o constants.o fault.o io.o model.o mpi_routines.o

main.o: constants.o io.o mpi_routines.o problem.o problem_routines.o
main.o: test_utilities.o

mesh.o: constants.o fft_routines_fftw3.o fields.o io.o model.o mpi_routines.o
mesh.o: utilities.o

model.o: constants.o io.o mpi_routines.o rk.o utilities.o

mpi_routines.o: constants.o

preslip.o: asperity.o constants.o io.o model.o mpi_routines.o utilities.o

problem.o: convolution.o fft_routines_fftw3.o fields.o friction.o front.o
problem.o: kernel.o mesh.o model.o

problem_routines.o: constants.o front.o init.o integration.o io.o mesh.o
problem_routines.o: mpi_routines.o problem.o rates.o static.o substep.o
problem_routines.o: thermpres.o time_step.o

rates.o: constants.o convolution_routines.o fault.o fourier.o
rates.o: friction_routines.o io.o load.o model.o problem.o thermpres.o

ratestate.o: asperity.o constants.o io.o model.o mpi_routines.o utilities.o

rk.o: constants.o io.o

slipweak.o: asperity.o constants.o io.o model.o mpi_routines.o

static.o: constants.o fault.o friction_routines.o io.o model.o problem.o
static.o: rates.o

substep.o: constants.o integration.o mpi_routines.o problem.o time_step.o

test_utilities.o: constants.o integration.o utilities.o

thermpres.o: asperity.o constants.o fault.o io.o model.o mpi_routines.o

time_step.o: constants.o fault.o fourier.o friction_routines.o integration.o
time_step.o: io.o load.o mpi_routines.o problem.o rates.o thermpres.o

utilities.o: constants.o io.o

fft_routines.mod: fft_routines_fftw3.o
runge_kutta.mod: rk.o
