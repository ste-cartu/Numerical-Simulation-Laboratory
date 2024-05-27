CC = g++
CFLAGS = -O3 --std=c++20
AFLAGS = -larmadillo -rpath /opt/homebrew/Cellar/armadillo-12.8.3
SIM = ../../../Libraries/Simulator/SOURCE/
RND = ../../../Libraries/RandomGen/
LIB = ../../../Libraries/
E = Equilibration/
S = Simulation/


solid : $(E)solid_eq $(S)solid_sim
eq : $(E)solid_eq
sim : $(S)solid_sim


$(E)solid_eq : $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(E)solid_eq.o
	$(CC) $(CFLAGS) $(AFLAGS) $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(E)solid_eq.o -o $(E)solid_eq

$(E)solid_eq.o : $(E)solid_eq.cpp $(SIM)system.hpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(E)solid_eq.cpp -o $(E)solid_eq.o

$(S)solid_sim : $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(S)solid_sim.o
	$(CC) $(CFLAGS) $(AFLAGS) $(SIM)system.o $(SIM)particle.o $(RND)random.o $(LIB)library.o $(S)solid_sim.o -o $(S)solid_sim

$(S)solid_sim.o : $(S)solid_sim.cpp $(SIM)system.hpp $(SIM)particle.hpp $(LIB)fmtlib.hpp
	$(CC) $(CFLAGS) -c $(S)solid_sim.cpp -o $(S)solid_sim.o


$(SIM)system.o : $(SIM)system.cpp $(SIM)system.hpp $(RND)random.hpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(SIM)system.cpp -o $(SIM)system.o

$(SIM)particle.o : $(SIM)particle.cpp $(SIM)particle.hpp
	$(CC) $(CFLAGS) -c $(SIM)particle.cpp -o $(SIM)particle.o

$(RND)random.o : $(RND)random.cpp $(RND)random.hpp
	$(CC) $(CFLAGS) -c $(RND)random.cpp -o $(RND)random.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.hpp
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o


doeq : $(E)solid_eq
	./$(E)solid_eq
dosim : $(S)solid_sim
	./$(S)solid_sim
clean :
	rm $(E)*.o $(S)*.o $(SIM)*.o $(RND)*.o
cleanall :
	rm $(SIM)*.o $(RND)*.o $(E)*.o $(S)*.o $(E)solid_eq $(S)solid_sim
rmeq :
	rm $(E)OUTPUT/*.*
	rm $(E)OUTPUT/CONFIG/*.*
rmsim :
	rm $(S)OUTPUT/*.*
	rm $(S)OUTPUT/CONFIG/*.*
rmall : rmeq rmsim