
CC = gcc
COM =  mrc.o
OBJECTS = $(COM)
HEADER = mrc.h
PG = Vectorize
Alias = Sample

OPT = -O3 -msse4.2 -msse -msse2 -msse3 -mmmx -ftree-vectorize -ftree-vectorizer-verbose=3 -ffast-math -march=native -funroll-loops -fopenmp

$(PG) : main.c $(OBJECTS) $(HEADER) $(COM)
	$(CC) -o $(Alias) main.c $(OBJECTS) -lm $(OPT) -lfftw3_threads -lfftw3f

mrc.o : mrc.c $(HEADER)
	$(CC) -O3 -c mrc.c -lm $(OPT)

