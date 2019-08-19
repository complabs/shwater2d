CC = gcc
CFLAGS = -O2 -fopenmp

ifeq ($(CRAY_PRGENVCRAY), loaded)
CFLAGS = -O2 -openmp
else ifeq ($(CRAY_PRGENVINTEL), loaded)
CFLAGS = -O2 -openmp -D_Float128=__float128
else ifeq ($(CRAY_PRGENVGNU), loaded)
CFLAGS = -O2 -fopenmp
endif

SRC = vtk_export.c shwater2d.c
DEST = shwater2d
OBJS = ${SRC:.c=.o}

SRC2 = vtk_export.c shwater2d_naive.c
DEST2 = shwater2d_naive
OBJS2= ${SRC2:.c=.o}

all: $(DEST) $(DEST2)

$(DEST): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm	

$(DEST2): $(OBJS2)
	$(CC) $(CFLAGS) $(OBJS2) -o $@ -lm	

clean:
	rm -f $(DEST) $(DEST2) *.o *.vtk

%.o: %.c
	$(CC) $(CFLAGS) -c $<

