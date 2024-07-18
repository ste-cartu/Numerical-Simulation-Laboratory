CC = g++
CFLAGS = -Wall -std=c++20 -O3
AFLAGS = -larmadillo -rpath /opt/homebrew/Cellar/armadillo-12.8.3
LIB = ../../Libraries/
RND = ../../Libraries/RandomGen/


main : main.o $(RND)random.o $(LIB)blockaverage.o $(LIB)metropolis.o $(LIB)library.o
	$(CC) $(CFLAGS) $(AFLAGS) $(RND)random.o $(LIB)blockaverage.o $(LIB)metropolis.o $(LIB)library.o main.o -o main

main.o : main.cpp $(RND)random.hpp $(LIB)blockaverage.hpp $(LIB)library.hpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)blockaverage.o : $(LIB)blockaverage.cpp $(LIB)blockaverage.hpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)blockaverage.cpp -o $(LIB)blockaverage.o

$(LIB)metropolis.o : $(LIB)metropolis.hpp $(LIB)metropolis.cpp $(RND)random.hpp $(LIB)functions.hpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c $(LIB)metropolis.cpp -o $(LIB)metropolis.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o


exe :
	./main 1e6 1e2
clean :
	rm *.o $(RND)*.o $(LIB)*.o main
all : main exe
