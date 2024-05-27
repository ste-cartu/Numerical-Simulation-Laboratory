CC = g++
CFLAGS = -O3 --std=c++20
AFLAGS = -larmadillo -rpath /opt/homebrew/Cellar/armadillo-12.8.3
SIM = ../../../Libraries/Simulator/SOURCE/
RND = ../../../Libraries/RandomGen/
LIB = ../../../Libraries/
E = Equilibration/
S = Simulation/


gas : $(E)gas_eq $(S)gas_sim


$(E)gas_eq : $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(E)gas_eq.o
	$(CC) $(CFLAGS) $(AFLAGS) $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(E)gas_eq.o -o $(E)gas_eq

$(E)gas_eq.o : $(E)gas_eq.cpp $(SIM)system.hpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(E)gas_eq.cpp -o $(E)gas_eq.o

$(S)gas_sim : $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(S)gas_sim.o
	$(CC) $(CFLAGS) $(AFLAGS) $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(S)gas_sim.o -o $(S)gas_sim

$(S)gas_sim.o : $(S)gas_sim.cpp $(SIM)system.hpp $(SIM)particle.hpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c $(S)gas_sim.cpp -o $(S)gas_sim.o


$(SIM)system.o : $(SIM)system.cpp $(SIM)system.hpp $(RND)random.hpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(SIM)system.cpp -o $(SIM)system.o

$(SIM)particle.o : $(SIM)particle.cpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(SIM)particle.cpp -o $(SIM)particle.o

$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o


doeq : $(E)gas_eq
	./$(E)gas_eq
dosim : $(S)gas_sim
	./$(S)gas_sim
clean :
	rm $(E)*.o $(S)*.o $(SIM)*.o $(RND)*.o
cleanall :
	rm $(SIM)*.o $(RND)*.o $(E)*.o $(S)*.o $(E)gas_eq $(S)gas_sim
rmeq :
	rm $(E)OUTPUT/*.*
	rm $(E)OUTPUT/CONFIG/*.*
rmsim :
	rm $(S)OUTPUT/*.*
	rm $(S)OUTPUT/CONFIG/*.*
rmall : rmeq rmsim