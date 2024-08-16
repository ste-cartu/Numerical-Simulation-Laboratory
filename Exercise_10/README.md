# Exercise 10

Here you can fnd $4$ directories, because I have implemented the solution to both the assignments of the exercises.
1. Directory **10.1** contains the parallelization of exercise **09.1** on the circumference and the square frameworks performing migrations.
2. Directory **10.2** contains the application of the parallelized genetic algorithm to the italian 'capoluoghi di provincia', where if the `migration_number` input parameter is set to $0$ migrations are not performed, or else if it is not zero it represents the number of generations between two migrations.
3. Directory **10.3** contains the application of the parallel tempering technique to the genetic algorithm.
4. Directory **10.4** contains a non parallel genetic algorithm that uses the number of temperatures (and so of parallel processes) of exercise **10.3** as the population size.

Exercises **10.1**, **10.2** and **10.3** uses the **MPI** library, so I used `mpic++` to compile them. To run this exercises, the commands you can use are:
- `make`: compile the source file
- `make exe NP=<n>` or `mpiexec -np <n> ./main`: run the executable using `n` parallel processes, the default value for `make exe` is `n`$=4$
- `make all NP=<n>`: compile the source and run the executable using `n` parallel processes, the default is `n`$=4$
- `make clean`: delete all objects and executable files
- `make remove`: delete all the output files

