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
    TSP salesman_circle(&rnd, rank);
    salesman_circle.Init();
    TSP salesman_square(&rnd, rank);
    salesman_square.Init();

    salesman_circle.SetType("circle");
    salesman_square.SetType("square");

    int n_cities = salesman_circle.GetLen();
    int pop_size = salesman_circle.GetDim();
    //int dist_ord = salesman_circle.GetNorm();
    int n_gens = salesman_circle.GetNGens();
    int n_migr = salesman_circle.GetNMigr();

    int circle_init[1] = {0}, square_init[1] = {0};
    if(rank == 0) {         // rank 0 initialize the cities on the circle
        salesman_circle.Circle();
        circle_init[0] = 1;
    }
    if(rank == 1) {         // rank 1 initialize the cities on the square
        salesman_square.Square();
        square_init[0] = 1;
    }

    // broadcast the message: now all other ranks can read the cities from file
    MPI_Bcast(circle_init, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(square_init, 1, MPI_INTEGER, 1, MPI_COMM_WORLD);

    if(rank != 0 and circle_init[0] == 1) {salesman_circle.FromFile("cities_circle.tsv");}
    if(rank != 1 and square_init[0] == 1) {salesman_square.FromFile("cities_square.tsv");}
    // now all ranks have the same cities

    // output files
    string opt_circle = "Output/opt_circle_rank=" + to_string(rank) + ".csv";
    ofstream out_circle(opt_circle);
    out_circle << "generation,best_loss,";
    string opt_square = "Output/opt_square_rank=" + to_string(rank) + ".csv";
    ofstream out_square(opt_square);
    out_square << "generation,best_loss,";
    for(int j=0 ; j<n_cities-1 ; j++) {
        out_circle << "city_" << j+1 << ",";
        out_square << "city_" << j+1 << ",";
    }
    out_circle << "city_" << n_cities << endl;
    out_square << "city_" << n_cities << endl;

    // store all the losses
    string loss_circle = "Output/loss_circle_rank=" + to_string(rank) + ".csv";
    ofstream l_cir(loss_circle);
    string loss_square = "Output/loss_square_rank=" + to_string(rank) + ".csv";
    ofstream l_squ(loss_square);
    l_cir << "generation,";
    l_squ << "generation,";
    for(int i=0 ; i<pop_size-1 ; i++) {
        l_cir << "path_" << i+1 << ",";
        l_squ << "path_" << i+1 << ",";
    }
    l_cir << "path_" << pop_size << endl;
    l_squ << "path_" << pop_size << endl;

    // matrices to store the best individuals for migration
    umat best_circle(n_cities, size);
    umat best_square(n_cities, size);

    // evolutive cycle
    for(int i=0 ; i<n_gens ; i++) {
        Progress_Bar(i, n_gens-1);
        if(n_migr != 0 and i%n_migr == 0 and i != 0) {
            salesman_square.Order();
            salesman_circle.Order();

            // all processes wait for each other, then gather all the best paths to process 0, shuffle them and send them back
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(salesman_circle[0].GetIndex().memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, best_circle.memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Gather(salesman_square[0].GetIndex().memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, best_square.memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

            if(rank == 0) {
                best_circle = shuffle(best_circle, 1);
                best_square = shuffle(best_square, 1);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(best_circle.memptr(), best_circle.size(), MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(best_square.memptr(), best_square.size(), MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

            salesman_circle[0].SetIndex(best_circle.col(rank));
            salesman_square[0].SetIndex(best_square.col(rank));
        }

        salesman_square.Selection();
        salesman_circle.Selection();

        // print optimization to file
        out_circle << i+1 << "," << salesman_circle[0].GetLoss() << ",";
        out_square << i+1 << "," << salesman_square[0].GetLoss() << ",";
        for(int j=0 ; j<n_cities-1 ; j++) {
            out_circle << salesman_circle[0][j] << ",";
            out_square << salesman_square[0][j] << ",";
        }
        out_circle << salesman_circle[0][n_cities-1] << endl;
        out_square << salesman_square[0][n_cities-1] << endl;

        // print losses to file
        l_cir << i+1 << ",";
        l_squ << i+1 << ",";
        for(int j=0 ; j<pop_size-1 ; j++) {
            l_cir << salesman_circle[j].GetLoss() << ",";
            l_squ << salesman_square[j].GetLoss() << ",";
        }
        l_cir << salesman_circle[pop_size-1].GetLoss() << endl;
        l_squ << salesman_square[pop_size-1].GetLoss() << endl;  

        salesman_circle.Mutations();
        salesman_square.Mutations();
    }
    if(rank == 0) {fmt::print("\n");}

    salesman_circle.GetPop().Order();
    salesman_square.GetPop().Order();

    // finally, put in process 0 the best of all paths
    double min_cir_local[1] = {salesman_circle[0].GetLoss()};
    double min_squ_local[1] = {salesman_square[0].GetLoss()};
    double min_cir_global[1];
    double min_squ_global[1];

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(min_cir_local, min_cir_global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(min_squ_local, min_squ_global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(min_cir_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(min_squ_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(salesman_circle[0].GetLoss() == min_cir_global[0]) {fmt::print("Best circle process: {}   Loss: {}\n", rank, min_cir_global[0]);}
    if(salesman_square[0].GetLoss() == min_squ_global[0]) {fmt::print("Best square process: {}   Loss: {}\n", rank, min_squ_global[0]);}

    MPI_Finalize();
    return 0;

}


