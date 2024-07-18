# Exercise 06

## Usage
To run exercise **06.1** you can move in its directory and then use the following commands:
- `make`: compile the source file
- `make metro`: execution with Mteropolis sampling at a fixed temperature
- `make gibbs`: execution with Gibbs sampling at a fixed temperature
- `make all`: execution with Mteropolis and then with Gibbs sampling at a fixed temperature
- `make clean`: delete all objects and executable files
- `make remove`: delete all the outputs of both the sampling methods
- `./run.sh metro`: perform simulations at all temperatures with Metropolis sampling
- `./run.sh gibbs`: perform simulations at all temperatures with Gibbs sampling
- `./run.sh metro gibbs`: perform simulations at all temperatures with Metropolis and Gibbs sampling

## Output
The equilibration outputs are stored in the `OUTPUT/EQUILIBRATION` directory, while the simulation outputs are in `OUTPUT/SIMULATION`. There are also two file `summary_H=0.00.csv` and `summary_H=0.02.csv` that contain the summary of all simulations.