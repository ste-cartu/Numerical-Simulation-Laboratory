# Exercise 10

Both exercises **10.1** and **10.2** uses the **MPI** library, so I used `mpic++` to compile them. To run this exercises, the commands you can use are:
- `make`: compile the source file
- `make exe NP=<n>` or `mpiexec -np <n> ./main`: run the executable using `n` parallel processes, the default value for `make exe` is `n`$=4$
- `make all NP=<n>`: compile the source and run the executable using `n` parallel processes, the default is `n`$=4$
- `make clean`: delete all objects and executable files
- `make remove`: delete all the output files

