CC = g++
CFLAGS = -Wall -O3 --std=c++20
AFLAGS = -larmadillo -rpath /opt/homebrew/Cellar/armadillo-12.8.3
SIM = ../../Libraries/Simulator/SOURCE/
RND = ../../Libraries/RandomGen/
LIB = ../../Libraries/


main : $(RND)random.o $(LIB)blockaverage.o $(LIB)metropolis.o $(LIB)library.o main.o
	$(CC) $(CFLAGS) $(AFLAGS) $(RND)random.o $(LIB)blockaverage.o $(LIB)metropolis.o $(LIB)library.o main.o -o main

main.o : main.cpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)blockaverage.o : $(LIB)blockaverage.cpp $(LIB)blockaverage.hpp
	$(CC) $(CFLAGS) -c $(LIB)blockaverage.cpp -o $(LIB)blockaverage.o

$(LIB)metropolis.o : $(LIB)metropolis.cpp $(LIB)metropolis.hpp
	$(CC) $(CFLAGS) -c $(LIB)metropolis.cpp -o $(LIB)metropolis.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o

exe :
	./main
all : main exe
clean :
	rm *.o main $(LIB)*.o $(RND)*.o