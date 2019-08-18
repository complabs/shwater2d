CC = gcc
CFLAGS = -O2 -fopenmp

ifeq ($(CRAY_PRGENVCRAY), loaded)
CFLAGS = -O2 -openmp
else ifeq ($(CRAY_PRGENVINTEL), loaded)
CFLAGS = -O2 -openmp -D_Float128=__float128
else ifeq ($(CRAY_PRGENVGNU), loaded)
CFLAGS = -O2 -fopenmp
endif

OBJS = ${SRC:.c=.o}

SRC = vtk_export.c shwater2d_naive.c
DEST = shwater2d_naive

SRC = vtk_export.c shwater2d.c
DEST = shwater2d

all: $(DEST)

$(DEST): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm	

clean:
	rm -f $(DEST) *.o

%.o: %.c
	$(CC) $(CFLAGS) -c $<

