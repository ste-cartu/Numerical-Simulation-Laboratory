# Numerical Simulation Laboratory
#### _Stefano Carturan - A.A. 2023-2024_
---

## Structure
In this repository it is stored the code I wrote for the Numerical Simulation Laboratory project. There is a folder for every exercise, and also one directory containing all the assignments and one for the libraries I used. Every exercise-folder contains a Jupyter Notebook showing the results and one or more subdirectories for different assignments (e.g. in the directory Exercise_01 there are the following subdirectories: 01.1, 01.2, 01.3).

## Prerequisites
To run the exercises contained in this repository you must have:
- A `C++` compiler supporting `C++20`
- `Armadillo` library
- `Openmp` librray (only for exercise 10)

## Usage
Unless specified otherwise, every exercise can be compiled typing `make` in the subdirecory corresponding to the assignment, and can be executed typing `make exe` or just `./main`. To compile and execute in one command you can type `make all` and to delete object and executable files the command is `make clean`.