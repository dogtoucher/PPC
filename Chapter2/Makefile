CC = g++
CFLAGS= -g -O3 -march=native -std=c++17
OMP_FLAG = -fopenmp
MAIN=shortcut.o
OBJ=v0.o v0-omp.o v1.o v1-omp.o v2.o
OUTS=v0.out v0-omp.out v1.out v1-omp.out v2.out

all: $(MAIN) $(OBJ) $(OUTS)

%.out: %.o
	$(CC) shortcut.o $< -o $@ $(CFLAGS) $(OMP_FLAG)

%.o: %.cc
	$(CC) -c -o $@ $< $(CFLAGS)

v0: $(MAIN)
	$(CC) shortcut.o v0.o $(CFLAGS) -o v0.out

clean:
	rm ./*.o ./*.out