CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
RND = ../../RandomGen/

all: simulator.exe

simulator.exe : system.o NSL_SIMULATOR.o particle.o $(RND)random.o ../../library.o
	$(CC) system.o particle.o $(RND)random.o NSL_SIMULATOR.o -o simulator.exe ../../library.o $(AFLAGS)
NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.hpp particle.hpp
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS)
system.o : system.cpp system.hpp $(RND)random.hpp particle.hpp
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.hpp
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) -c $(RND)random.cpp -o $(RND)random.o $(CFLAGS)
../../library.o : ../../library.cpp ../../library.hpp
	$(CC) -c ../../library.cpp -o ../../library.o $(CFLAGS)

clean :
	rm *.o rm $(RND)*.o simulator.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
