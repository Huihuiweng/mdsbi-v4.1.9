# Set build type

#BUILD = debug
#BUILD = default
BUILD = production
#BUILD = openmp

# Compiler options for different machines

Exe = ../mdsbi

F95 = mpiifort
 ifeq ($(BUILD),debug)
  F95FLAGS = -g -check all -check noarg_temp_created -warn all -traceback -ftrapuv -fpe0
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -O2
 endif
 ifeq ($(BUILD),openmp)
   F95FLAGS = -o2 -OPENMP -PARALLEL
 endif
 FFTW2 = -L/home/weng/Works/Softwares/fftw-3.3.10/lib/ -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw \
	 -L/home/weng/Works/Softwares/lapack-3.9.0 -llapack
 FFTW3 = -L/home/weng/Works/Softwares/fftw-3.3.10/lib/ -lfftw3 \
	  -llapack
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw2.f90),)
	LIBS = $(FFTW2) 
 endif
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw2c.f90),)
	LIBS = $(FFTW2)
 endif
 ifeq ($(shell diff fft_routines.f90 fft_routines_fftw3.f90),)
	LIBS = $(FFTW3)
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
