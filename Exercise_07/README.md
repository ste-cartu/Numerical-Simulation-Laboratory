# Exercise 07

## 07.2 - Usage
To run exercise **07.2** you can use the following commands:
- `make`: compile the source file
- `make sol` or `./main Solid` perform **equilibration** and **simulation** for the solid phase
- `make liq` or `./main Liquid` perform **equilibration** and **simulation** for the liquid phase
- `make gas` or `./main Gas` perform **equilibration** and **simulation** for the gas phase
- `make all` perform **equilibration** and **simulation** for all phases
- `make clean`: delete all objects and executable files
- `make rmsol`: delete all output files relative to the solid phase
- `make rmliq`: delete all output files relative to the liquid phase
- `make rmgas`: delete all output files relative to the gas phase


## 07.4 - Usage
To run exercise **07.4** you can use the following commands:
- `make`: compile the source file 
- `make clean`: delete all objects and executable files
- `./run.sh NVE`: perform **NVE** equilibration and simulation over all the phases: solid, liquid and gas
- `./run.sh NVT`: perform **NVT** equilibration and simulation over all the phases: solid, liquid and gas
- `./run.sh NVE NVT`: perform both **NVE** and **NVT** equilibration and simulation over all the phases