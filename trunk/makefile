OBJS = ED_Lan_1107.o GenHam.o  Lanczos_07.o lapack.o Lattice_30R.cpp
CC = g++
#CFLAGS = -O2 
CFLAGS = -O2 -arch x86_64
#LIBS = -lm -framework veclib
LIBS = -framework Accelerate

a.out: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o a.out $(LIBS)

ED_Lan_1107.o: ED_Lan_1107.cpp GenHam.h Lanczos_07.h lapack.h simparam.h
	$(CC) $(CFLAGS) -c ED_Lan_1107.cpp

GenHam.o: GenHam.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c GenHam.cpp

Lattice_30R.o: Lattice_30R.cpp GenHam.h 
	$(CC) $(CFLAGS) -c Lattice_30R.cpp

Lanczos_07.o: Lanczos_07.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c Lanczos_07.cpp

lapack.o: lapack.cpp lapack.h 
	$(CC) $(CFLAGS) -c lapack.cpp
clean :
	rm *.o
