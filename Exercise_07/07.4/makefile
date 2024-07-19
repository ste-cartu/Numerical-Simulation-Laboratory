CC = g++
CFLAGS = -Wall -O3 --std=c++20
AFLAGS = -larmadillo -rpath /opt/homebrew/Cellar/armadillo-12.8.3
SIM = ../../Libraries/Simulator/SOURCE/
RND = ../../Libraries/RandomGen/
LIB = ../../Libraries/


main : $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o main.o
	$(CC) $(CFLAGS) $(AFLAGS) $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o main.o -o main

main.o : main.cpp $(SIM)system.hpp $(SIM)particle.hpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(SIM)system.o : $(SIM)system.cpp $(SIM)system.hpp $(RND)random.hpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(SIM)system.cpp -o $(SIM)system.o

$(SIM)particle.o : $(SIM)particle.cpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(SIM)particle.cpp -o $(SIM)particle.o

$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o

clean :
	rm *.o $(RND)*.o $(SIM)*.o $(LIB)*.o main