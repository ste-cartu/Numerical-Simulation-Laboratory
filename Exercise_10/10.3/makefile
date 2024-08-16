CC = mpic++
CFLAGS = -Wall -O3 --std=c++20
AFLAGS = -larmadillo -rpath /opt/homebrew/Cellar/armadillo-12.8.3
SIM = ../../Libraries/Simulator/SOURCE/
RND = ../../Libraries/RandomGen/
LIB = ../../Libraries/
NP = 8


main : $(RND)random.o $(LIB)geneticalg.o $(LIB)library.o main.o
	$(CC) $(CFLAGS) $(AFLAGS) $(RND)random.o $(LIB)geneticalg.o $(LIB)library.o main.o -o main

main.o : main.cpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)geneticalg.o : $(LIB)geneticalg.cpp $(LIB)geneticalg.hpp
	$(CC) $(CFLAGS) -c $(LIB)geneticalg.cpp -o $(LIB)geneticalg.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o

exe :
	mpiexec -np $(NP) ./main
all : main exe
clean :
	rm *.o main $(LIB)*.o $(RND)*.o main
remove :
	rm Output/*.*