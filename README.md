# Numerical Simulation Laboratory
#### _Stefano Carturan - A.A. 2023-2024_
---

## Structure
In this repository it is stored the code I wrote for the Numerical Simulation Laboratory project. There is a folder for every exercise, and also one directory containing all the assignments and one for the libraries I used. Every exercise-folder contains a Jupyter Notebook showing the results and one or more subdirectories for different assignments (e.g. in the directory Exercise_01 there are the following subdirectories: 01.1, 01.2, 01.3). Typically in eanch subdirectory there is a single source file `main.cpp` and a `makefile`.

## Prerequisites
To run the exercises contained in this repository you must have:
- A `C++` compiler supporting `C++20`
- `Armadillo` library: please note that in the `makefiles` contained here this library is linked following the path `/opt/homebrew/Cellar/armadillo-12.8.3`. If you have it installed anywhere else on your computer you should change this path.
- `Openmp` librray (only for exercise 10)

## Usage
Unless specified otherwise, these are the makefile commands you can use to run all exercises:
- `make`: compile the source file
- `make exe` or `./main`: run the executable file
- `make all`: compile and run all at once
- `make clean`: delete all object and executable files
- `make remove`: delete all output files of a simulation