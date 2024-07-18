CC = g++
CFLAGS = -Wall -O3 --std=c++20
AFLAGS = -larmadillo -rpath /opt/homebrew/Cellar/armadillo-12.8.3
SIM = ../../../Libraries/Simulator/SOURCE/
RND = ../../../Libraries/RandomGen/
LIB = ../../../Libraries/
E = Equilibration/
S = Simulation/


liquid : $(E)liquid_eq $(S)liquid_sim
eq : $(E)liquid_eq
sim : $(S)liquid_sim


$(E)liquid_eq : $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(E)liquid_eq.o
	$(CC) $(CFLAGS) $(AFLAGS) $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(E)liquid_eq.o -o $(E)liquid_eq

$(E)liquid_eq.o : $(E)liquid_eq.cpp $(SIM)system.hpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(E)liquid_eq.cpp -o $(E)liquid_eq.o

$(S)liquid_sim : $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(S)liquid_sim.o
	$(CC) $(CFLAGS) $(AFLAGS) $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(S)liquid_sim.o -o $(S)liquid_sim

$(S)liquid_sim.o : $(S)liquid_sim.cpp $(SIM)system.hpp $(SIM)particle.hpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c $(S)liquid_sim.cpp -o $(S)liquid_sim.o


$(SIM)system.o : $(SIM)system.cpp $(SIM)system.hpp $(RND)random.hpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(SIM)system.cpp -o $(SIM)system.o

$(SIM)particle.o : $(SIM)particle.cpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(SIM)particle.cpp -o $(SIM)particle.o

$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o


doeq : $(E)liquid_eq
	./$(E)liquid_eq
dosim : $(S)liquid_sim
	./$(S)liquid_sim
clean :
	rm $(E)*.o $(S)*.o $(SIM)*.o $(RND)*.o $(E)liquid_eq $(S)liquid_sim
rmeq :
	rm $(E)OUTPUT/*.*
	rm $(E)OUTPUT/CONFIG/*.*
rmsim :
	rm $(S)OUTPUT/*.*
	rm $(S)OUTPUT/CONFIG/*.*
rmall : rmeq rmsim