#include <fstream>
#include <cmath>
#include <mpi.h>

#include "../../Libraries/geneticalg.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/fmtlib.hpp"


int main(int argc, char* argv[]) {

    int size, rank;

    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if(size < 2 or size > 11){
        fmt::print("ERROR! This program accept from 2 to 11 parallel processes!\n\n");
        MPI_Finalize();
        return 1;
    }

    Random rnd("../../Libraries/RandomGen/", rank);

    // initialization of the problem
    TSP salesman(&rnd, rank);
    salesman.Init();

    int pop_init[1] = {0};
    if(rank == 0) {         // rank 0 initialize the cities
        salesman.InitPopulation();
        pop_init[0] = 1;
    }

    // defining some variables from initialization
    string type = salesman.GetType();
    int n_cities = salesman.GetLen();
    int pop_size = salesman.GetDim();
    int n_gens = salesman.GetNGens();
    int n_migr = salesman.GetNMigr();

    // broadcast the message: now all other ranks can read the cities from file
    MPI_Bcast(pop_init, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    string cities_file = "cities_" + type + ".tsv";
    if(rank != 0 and pop_init[0] == 1) {salesman.FromFile(cities_file);}
    // now all ranks have the same cities

    // output file
    string optim = "Output/opt_" + type + "_rank=" + to_string(rank) + "_migr=" + to_string(n_migr) + ".csv";
    ofstream out(optim);
    out << "generation,best_loss,";
    for(int j=0 ; j<n_cities-1 ; j++) {out << "city_" << j+1 << ",";}
    out << "city_" << n_cities << endl;

    // store all the losses
    string losses = "Output/loss_" + type + "_rank=" + to_string(rank) + "_migr=" + to_string(n_migr) + ".csv";
    ofstream los(losses);
    los << "generation,";
    for(int j=0 ; j<pop_size-1 ; j++) {los << "path_" << j+1 << ",";}
    los << "path_" << pop_size << endl;

    // matrix to store the best individuals for migration
    umat best(n_cities, size);

    // evolutive cycle
    for(int i=0 ; i<n_gens ; i++) {
        Progress_Bar(i, n_gens-1);
        if(n_migr != 0 and i%n_migr == 0 and i != 0) {
            salesman.Order();

            // all processes wait for each other, then gather all the best paths to process 0, shuffle them and send them back
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(salesman[0].GetIndex().memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, best.memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

            if(rank == 0) {best = shuffle(best, 1);}

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(best.memptr(), best.size(), MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            salesman[0].SetIndex(best.col(rank));
        }

        salesman.Selection();

        // print optimization to file
        out << i+1 << "," << salesman[0].GetLoss() << ",";
        for(int j=0 ; j<n_cities-1 ; j++) {out << salesman[0][j] << ",";}
        out << salesman[0][n_cities-1] << endl;

        // print losses to file
        los << i+1 << ",";
        for(int j=0 ; j<pop_size-1 ; j++) {los << salesman[j].GetLoss() << ",";}
        los << salesman[pop_size-1].GetLoss() << endl;

        salesman.Mutations(); 
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0) {fmt::print("\n");}
 
    salesman.GetPop().Order();

    // finally, put in process 0 the best of all paths
    double min_local[1] = {salesman[0].GetLoss()};
    double min_global[1];

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(min_local, min_global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(min_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(salesman[0].GetLoss() == min_global[0]) {fmt::print("\nBest process: {}   Loss: {}\n\n", rank, min_global[0]);}

    MPI_Finalize();
    return 0;
}


