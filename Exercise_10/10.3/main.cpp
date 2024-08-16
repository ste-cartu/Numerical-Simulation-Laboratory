#include <fstream>
#include <cmath>
#include <mpi.h>
#include <vector>
#include <string>

#include "../../Libraries/geneticalg.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/fmtlib.hpp"


int main(int argc, char* argv[]) {

    int size, rank;

    // initialization of the parallel processes
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(size < 2 or size > 11){
        fmt::print("ERROR! This program accept from 2 to 11 parallel processes! e in più sono stitico\n\n");
        MPI_Finalize();
        return 1;
    }

    Random rnd("../../Libraries/RandomGen/", rank);

    // output file to store the best losses
    ofstream out_best("best_processes.txt");
    out_best.close();

    // loop over three problem types: circle, square and italy
    vector<string> types = {"circle", "square", "italy"};
    for(auto type = types.begin() ; type != types.end() ; type++) {
        
        // initialization of the problem
        TSP salesman(&rnd, rank);
        
        string filename = "input_" + *type + ".txt";
        salesman.Init(filename);

        int n_cities = salesman.GetLen();
        int pop_size = salesman.GetDim();
        int n_gens = salesman.GetNGens();
        string cities_file = salesman.GetCitiesFile();
        double temp_i = salesman.GetTempHi();
        double temp_f = salesman.GetTempLo();

        // initialize temperatures for each process
        double delta_temp = fabs(temp_i - temp_f) / (size-1);
        double temp = rank*delta_temp + min(temp_i, temp_f);
        salesman.SetTemp(temp);

        // rank 0 initializes the cities
        bool init = false;
        if(rank == 0) {
            salesman.InitPopulation();
            init = true;
        }

        // broadcast the message: now all other ranks can read the cities from file
        MPI_Bcast(&init, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
        if(rank != 0 and init == true) {salesman.FromFile(cities_file);}
        // now all processes have the same cities

        // output files
        string opt = "Output/opt_" + *type + "_rank=" + to_string(rank) + ".csv";
        ofstream out(opt);
        out << "generation,best_loss,";
        for(int j=0 ; j<n_cities-1 ; j++) {out << "city_" << j+1 << ",";}
        out << "city_" << n_cities << endl;

        // store all the losses
        string loss = "Output/loss_" + *type + "_rank=" + to_string(rank) + ".csv";
        ofstream l(loss);
        l << "generation,";
        for(int i=0 ; i<pop_size-1 ; i++) {l << "path_" << i+1 << ",";}
        l << "path_" << pop_size << endl;

        // evolutive cycle
        for(int i=0 ; i<n_gens ; i++) {
            Progress_Bar(i, n_gens-1);

            salesman.Selection();
            // print optimization to file
            out << i+1 << "," << salesman[0].GetLoss() << ",";
            for(int j=0 ; j<n_cities-1 ; j++) {out << salesman[0][j] << ",";}
            out << salesman[0][n_cities-1] << endl;

            // print losses to file
            l << i+1 << ",";
            for(int j=0 ; j<pop_size-1 ; j++) {l << salesman[j].GetLoss() << ",";}
            l << salesman[pop_size-1].GetLoss() << endl;

            salesman.Mutations();

            // declaring variables for exchanging between processes
            uvec index_send = salesman[0].GetIndex();
            uvec index_recv(n_cities);
            int tag_index = 0;

            double loss_send = salesman[0].GetLoss();
            double loss_recv = 0;
            int tag_loss = 1;

            double temp_send = salesman.GetTemp();
            double temp_recv = 0;
            int tag_temp = 2;

            // proposing exchanges between processes with different temperatures
            if(i % 2 == 0) {
                if(rank % 2 == 0 and rank != 0) {
                    MPI_Send(index_send.memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, rank-1, tag_index, MPI_COMM_WORLD);
                    MPI_Send(&loss_send, 1, MPI_DOUBLE, rank-1, tag_loss, MPI_COMM_WORLD);
                    MPI_Send(&temp_send, 1, MPI_DOUBLE, rank-1, tag_temp, MPI_COMM_WORLD);
                }
                else if(rank % 2 != 0 and rank != size-1) {
                    MPI_Status status;
                    MPI_Recv(index_recv.memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, rank+1, tag_index, MPI_COMM_WORLD, &status);
                    MPI_Recv(&loss_recv, 1, MPI_DOUBLE, rank+1, tag_loss, MPI_COMM_WORLD, &status);
                    MPI_Recv(&temp_recv, 1, MPI_DOUBLE, rank+1, tag_temp, MPI_COMM_WORLD, &status);
                }
            }
            else {
                if(rank % 2 != 0) {
                    MPI_Send(index_send.memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, rank-1, tag_index, MPI_COMM_WORLD);
                    MPI_Send(&loss_send, 1, MPI_DOUBLE, rank-1, tag_loss, MPI_COMM_WORLD);
                    MPI_Send(&temp_send, 1, MPI_DOUBLE, rank-1, tag_temp, MPI_COMM_WORLD);

                }
                else if(rank % 2 == 0 and rank != size-1) {
                    MPI_Status status;
                    MPI_Recv(index_recv.memptr(), n_cities, MPI_UNSIGNED_LONG_LONG, rank+1, tag_index, MPI_COMM_WORLD, &status);
                    MPI_Recv(&loss_recv, 1, MPI_DOUBLE, rank+1, tag_loss, MPI_COMM_WORLD, &status);
                    MPI_Recv(&temp_recv, 1, MPI_DOUBLE, rank+1, tag_temp, MPI_COMM_WORLD, &status);
                }
            }

            // accepting or rejecting the proposed exchange
            double acc = min(1., exp(- (1./temp_send - 1./temp_recv) * (loss_recv - loss_send)));
            if(rnd.Rannyu() < acc) {
                salesman[0].SetIndex(index_recv);
                salesman[0].SetLoss(loss_recv);
            }

        }
        if(rank == 0) {fmt::print("\n");}

        // print best loss
        out_best.open("best_processes.txt", ios::app);
        //ofstream out_best("best_processes.txt", ios::app);
        if(rank == 0) {
            string best_str = fmt::format("Best {} loss: {:.5f}\n", *type, salesman[0].GetLoss());
            fmt::print("{}", best_str);
            out_best << best_str;
        }
        out_best.close();

    }

    MPI_Finalize();
    return 0;

}


