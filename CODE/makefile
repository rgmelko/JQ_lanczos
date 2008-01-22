OBJS = ED_Lan_1107.o GenHam.o  Lanczos_07.o lapack.o Lattice_16B.cpp Lattice_16A.cpp Lattice_18A.cpp Lattice_20A.cpp  Lattice_24A.cpp 
CC = g++
CFLAGS = -O3
LIBS = -lm -framework veclib

a.out: $(OBJS)
	$(CC) $(OBJS) -o a.out $(LIBS)

ED_Lan_1107.o: ED_Lan_1107.cpp GenHam.h Lanczos_07.h lapack.h
	$(CC) $(CFLAGS) -c ED_Lan_1107.cpp

GenHam.o: GenHam.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c GenHam.cpp

Lattice_20A.o: Lattice_20A.cpp GenHam.h 
	$(CC) $(CFLAGS) -c Lattice_20A.cpp

Lanczos_07.o: Lanczos_07.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c Lanczos_07.cpp

lapack.o: lapack.cpp lapack.h 
	$(CC) $(CFLAGS) -c lapack.cpp
clean :
	rm *.o
