CC = g++ 
CFLAGS = -Wall -std=c++20 -O3
LIB = ../../Libraries/
RND = ../../Libraries/RandomGen/


main : $(RND)random.o classes.o $(LIB)library.o main.o
	$(CC) $(CFLAGS) $(RND)random.o classes.o $(LIB)library.o main.o -o main

main.o : main.cpp $(RND)random.hpp $(LIB)classes.hpp $(LIB)library.hpp $(LIB)fmtlib.h
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(RND)random.o : $(RND)random.hpp $(RND)random.cpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

classes.o : $(LIB)classes.hpp classes.hpp classes.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c classes.cpp -o classes.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o


exe :
	./main 1 1e4 1e6 1e2
clean :
	rm *.o $(RND)*.o $(LIB)*.o
cleanall :
	rm *.o $(RND)*.o $(LIB)*.o main
all : main exe
