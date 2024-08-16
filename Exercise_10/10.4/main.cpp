#include <fstream>
#include <cmath>

#include "../../Libraries/geneticalg.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/fmtlib.hpp"


int main() {

    Random rnd("../../Libraries/RandomGen/");

    // output file to store the best losses
    ofstream out_best("best_processes.txt");
    out_best.close();

    // loop over three problem types: circle, square and italy
    vector<string> types = {"circle", "square", "italy"};
    for(auto type = types.begin() ; type != types.end() ; type++) {
        
        // initialization of the problem
        TSP salesman(&rnd);
        
        string filename = "input_" + *type + ".txt";
        salesman.Init(filename);

        int n_cities = salesman.GetLen();
        int pop_size = salesman.GetDim();
        int n_gens = salesman.GetNGens();

        // initialize the cities
        salesman.InitPopulation();

        // output files
        string opt = "Output/opt_" + *type + ".csv";
        ofstream out(opt);
        out << "generation,best_loss,";
        for(int j=0 ; j<n_cities-1 ; j++) {out << "city_" << j+1 << ",";}
        out << "city_" << n_cities << endl;

        // store all the losses
        string loss = "Output/loss_" + *type + ".csv";
        ofstream l(loss);
        l << "generation,";
        for(int i=0 ; i<pop_size-1 ; i++) {l << "path_" << i+1 << ",";}
        l << "path_" << pop_size << endl;

        // evolutive cycle
        fmt::print("\nOptimizing {}\n", *type);
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
        }
        
        fmt::print("\n");

        // print best loss
        out_best.open("best_processes.txt", ios::app);
        string best_str = fmt::format("Best {} loss: {:.5f}\n", *type, salesman[0].GetLoss());
        fmt::print("{}\n", best_str);
        out_best << best_str;
        out_best.close();

    }

    return 0;

}
