CC = g++ 
CFLAGS = -Wall -std=c++20 -O3
LIB = ../../Libraries/
RND = ../../Libraries/RandomGen/


main : main.o $(RND)random.o $(LIB)classes.o $(LIB)library.o
	$(CC) $(CFLAGS) $(RND)random.o $(LIB)classes.o $(LIB)library.o main.o -o main

main.o : main.cpp $(RND)random.hpp $(LIB)classes.hpp $(LIB)library.hpp $(LIB)fmtlib.h
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)classes.o : $(LIB)classes.cpp $(LIB)classes.hpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)classes.cpp -o $(LIB)classes.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o


exe :
	./main 1e6 1e2
clean :
	rm *.o $(RND)*.o $(LIB)*.o
cleanall :
	rm *.o $(RND)*.o $(LIB)*.o main
all : main exe
