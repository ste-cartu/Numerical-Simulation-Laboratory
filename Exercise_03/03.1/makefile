CC = g++ 
CFLAGS = -Wall -std=c++20 -O3
LIB = ../../Libraries/
RND = ../../Libraries/RandomGen/


main : $(RND)random.o $(LIB)blockaverage.o $(LIB)library.o main.o
	$(CC) $(CFLAGS) $(RND)random.o $(LIB)blockaverage.o $(LIB)library.o main.o -o main

main.o : main.cpp $(RND)random.hpp $(LIB)blockaverage.hpp $(LIB)library.hpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(RND)random.o : $(RND)random.hpp $(RND)random.cpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)blockaverage.o : $(LIB)blockaverage.cpp $(LIB)blockaverage.hpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)blockaverage.cpp -o $(LIB)blockaverage.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o


exe :
	./main
clean :
	rm *.o $(RND)*.o $(LIB)*.o
cleanall :
	rm *.o $(RND)*.o $(LIB)*.o main
all : main exe
