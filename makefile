OBJS = ED_Lan_FFhoney.o GenHam.o  Lanczos_07.o lapack.o Lattice_18_PBC.cpp
CC = g++
#CFLAGS = -g
CFLAGS = -O2 -arch x86_64
#LIBS = -lm -framework veclib
LIBS = -framework Accelerate

a.out: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o a.out $(LIBS)

ED_Lan_FFhoney.o: ED_Lan_FFhoney.cpp GenHam.h Lanczos_07.h lapack.h simparam.h
	$(CC) $(CFLAGS) -c ED_Lan_FFhoney.cpp

GenHam.o: GenHam.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c GenHam.cpp

Lattice_18_PBC.o: Lattice_18_PBC.cpp GenHam.h 
	$(CC) $(CFLAGS) -c Lattice_18_PBC.cpp

Lanczos_07.o: Lanczos_07.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c Lanczos_07.cpp

lapack.o: lapack.cpp lapack.h 
	$(CC) $(CFLAGS) -c lapack.cpp
clean :
	rm *.o
