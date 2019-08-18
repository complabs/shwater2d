CXX = g++
CXXFLAGS = -O2 -fopenmp
CPPFLAGS = -DEXPORT_VTK

ifeq ($(CRAY_PRGENVCRAY), loaded)
CXX = CC
CXXFLAGS = -O2 -openmp
else ifeq ($(CRAY_PRGENVINTEL), loaded)
CXX = CC
CXXFLAGS = -O2 -openmp -D_Float128=__float128
else ifeq ($(CRAY_PRGENVGNU), loaded)
CXX = CC
CXXFLAGS = -O2 -fopenmp
endif

OBJS = ${SRC:.cpp=.o}

SRC = shwater2d.cpp
DEST = shwater2d

all: $(DEST)

$(DEST): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ -lm	

clean:
	rm -f $(DEST) *.o
	rm -f bin obj *.vtk

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

