CXX = g++
CXXFLAGS = -std=c++14 -Wall -O2 -fopenmp

ifeq ($(CRAY_PRGENVCRAY), loaded)
   CXX = CC
   CXXFLAGS = -std=c++14 -Wall -O2 -openmp
else ifeq ($(CRAY_PRGENVINTEL), loaded)
   CXX = CC
   CXXFLAGS = -std=c++14 -Wall -O2 -openmp -D_Float128=__float128
else ifeq ($(CRAY_PRGENVGNU), loaded)
   CXX = CC
   CXXFLAGS = -std=c++14 -Wall -O2 -fopenmp
endif

ifeq (, $(shell which srun))
   SRUN = 
else
   SRUN = srun -n 1   
endif

all: shwater2d_ser shwater2d_par shwater2d_opt
	@$(MAKE) -C orig
	@#$(MAKE) -C f90

shwater2d_ser : shwater2d_ser.cpp
	$(CXX) $(CXXFLAGS) shwater2d_ser.cpp -o shwater2d_ser

shwater2d_par : shwater2d_par.cpp
	$(CXX) $(CXXFLAGS) shwater2d_par.cpp -o shwater2d_par

shwater2d_opt : shwater2d_opt.cpp
	$(CXX) $(CXXFLAGS) shwater2d_opt.cpp -o shwater2d_opt

clean:
	rm -f shwater2d_ser shwater2d_par shwater2d_opt *.o
	rm -f bin obj *.vtk
	@$(MAKE) -C orig clean
	@$(MAKE) -C f90 clean

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

test-par: shwater2d
	$(SRUN) ./shwater2d 0_par 1000 1000 0.1 1

test-ser: shwater2d_opt
	$(SRUN) ./shwater2d_ser 0 1000 1000 0.1 1

test-opt: shwater2d_opt
	$(SRUN) ./shwater2d_opt 0 1000 1000 0.1 1

test-ref:
	@$(MAKE) -C orig all
	( cd orig; $(SRUN) ./shwater2d_naive )

diff-ref:
	diff result.vtk orig/result.vtk

