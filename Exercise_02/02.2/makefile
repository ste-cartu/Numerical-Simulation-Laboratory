CC = g++
CFLAGS = -Wall -std=c++20 -O3
AFLAGS = -larmadillo -rpath /opt/homebrew/Cellar/armadillo-12.8.3
LIB = ../../Libraries/
RND = ../../Libraries/RandomGen/


main : $(RND)random.o $(LIB)randomwalk.o $(LIB)library.o main.o
	$(CC) $(CFLAGS) $(AFLAGS) $(RND)random.o $(LIB)randomwalk.o $(LIB)library.o main.o -o main

main.o : main.cpp $(RND)random.hpp $(LIB)blockaverage.hpp $(LIB)library.hpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(RND)random.o : $(RND)random.hpp $(RND)random.cpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)randomwalk.o : $(LIB)blockaverage.hpp $(LIB)randomwalk.hpp $(LIB)randomwalk.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)randomwalk.cpp -o $(LIB)randomwalk.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o


exe :
	./main 1 1e3 1e5 1e2
clean :
	rm *.o $(RND)*.o $(LIB)*.o main
all : main exe
