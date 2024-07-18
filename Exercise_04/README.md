# Exercise 04

The exercise **04.2** can be executed in two ways.

## First method: equilibration and simulation
The first way to compute the exercise is divided in two steps: **equilibration** and **simulation**, performed by two different source file for each matter phase: solid, liquid or gas. This exploits the `RESTART` option of the simulator.

To run the exercise you have move to the directory corresponding the phase you want to simulate. Then you have several options:
- `make eq`: compile the **equilibration** source file
- `make sim`: compile the **simulation** source file
- `make` or `make <phase>`: compile both the **equilibration** and the **simulation** source files
- `make doeq`: perform the **equilibration**
- `make dosim`: perform the **simulation**
- `make clean`: delete all objects and executable files
- `make rmeq`: delete all the outputs of the **equilibration**
- `make rmsim`: delete all the outputs of the **simulation**
- `make rmall`: delete all the outputs of the both the **equilibration** and the **simulation**

**ATTENTION** if you do this way you must manually set the input temperature in file `<phase>/Simulation/INPUT/input.dat` to the correct value, obtained as an output in the file `<phase>/Equilibration/OUTPUT/equilibration.dat`. You must also copy the file `velocities.dat` from directory `<phase>/Equilibration/OUTPUT/CONFIG` to `<phase>/Simulation/INPUT/CONFIG`.


## Second method: all at once
If you want to perform **equilibration** and **simulation** straight away you can just launch the commands in the directory **04.2**. Here there is a single source file and the possible shell commands are:
- `make`: compile the source file
- `make sol` or `./main Solid` perform **equilibration** and **simulation** for the solid phase
- `make liq` or `./main Liquid` perform **equilibration** and **simulation** for the liquid phase
- `make gas` or `./main Gas` perform **equilibration** and **simulation** for the gas phase
- `make all` perform **equilibration** and **simulation** for all phases
- `make clean`: delete all objects and executable files


## Note: equilibration outputs
Every time you perform an equilibration, the program will store all physical quantities you measure during this phase in the folder `OUTPUT/EQUILIBRATION`.