FC = gfortran
FFLAGS = -O2 -fopenmp

ifeq ($(CRAY_PRGENVCRAY), loaded)
FC = ftn
FFLAGS = -h omp -O2 -openmp
else ifeq ($(CRAY_PRGENVINTEL), loaded)
FC = ftn
FFLAGS = -O2 -openmp
else ifeq ($(CRAY_PRGENVGNU), loaded)
ifeq ($(shell expr $(GCC_VERSION) '<' 5.0), 1)
$(error Unsupported GCC version, use at least v5.x (module swap gcc gcc/5.1.0))
endif
FC = ftn
FFLAGS = -O2 -fopenmp
endif

SRC = vtk_export.f90 shwater2d.f90
SRC = vtk_export.f90 shwater2d_naive.f90
SRC = vtk_export.f90 shwater2d_opt_transpose.f90
SRC = vtk_export.f90 shwater2d_opt_transpose_and_numerics.f90

OBJS = ${SRC:.f90=.o}
DEST = shwater2d

all: $(DEST)

$(DEST): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@	

clean:
	rm -f $(DEST) *.mod *.MOD *.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

