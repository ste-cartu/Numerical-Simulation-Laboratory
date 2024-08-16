#include <fstream>
#include <cmath>

#include "../../Libraries/geneticalg.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/fmtlib.hpp"


int main() {

    fmt::print("\n");

    Random rnd("../../Libraries/RandomGen/", 0);

    // initialization of the problem
    TSP salesman_circle(&rnd);
    salesman_circle.Init();
    TSP salesman_square(&rnd);
    salesman_square.Init();
    int n_cities = salesman_circle.GetLen();
    int pop_size = salesman_circle.GetDim();
    int n_gens = salesman_circle.GetNGens();

    salesman_circle.SetType("circle");
    salesman_square.SetType("square");

    salesman_circle.Circle();
    salesman_square.Square();
    ofstream out_circle;
    ofstream out_square;

    // output file
    out_circle.open("opt_circle.csv");
    out_circle << "generation,best_loss,";
    out_square.open("opt_square.csv");
    out_square << "generation,best_loss,";
    for(int j=0 ; j<n_cities-1 ; j++) {
        out_circle << "city_" << j+1 << ",";
        out_square << "city_" << j+1 << ",";
    }
    out_circle << "city_" << n_cities << endl;
    out_square << "city_" << n_cities << endl;

    // store all the losses
    ofstream l_cir("loss_circle.csv");
    ofstream l_squ("loss_square.csv");
    l_cir << "generation,";
    l_squ << "generation,";
    for(int i=0 ; i<pop_size-1 ; i++) {
        l_cir << "path_" << i+1 << ",";
        l_squ << "path_" << i+1 << ",";
    }
    l_cir << "path_" << pop_size << endl;
    l_squ << "path_" << pop_size << endl;

    // evolutive cycle
    for(int i=0 ; i<n_gens ; i++) {
        Progress_Bar(i, n_gens-1);
        salesman_square.Selection();
        salesman_circle.Selection();

        out_circle << i+1 << "," << salesman_circle[0].GetLoss() << ",";
        out_square << i+1 << "," << salesman_square[0].GetLoss() << ",";
        for(int j=0 ; j<n_cities-1 ; j++) {
            out_circle << salesman_circle[0][j] << ",";
            out_square << salesman_square[0][j] << ",";
        }
        out_circle << salesman_circle[0][n_cities-1] << endl;
        out_square << salesman_square[0][n_cities-1] << endl;

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
 
    salesman_circle.Order();
    salesman_square.Order();

    fmt::print("\n\nSmallest distance on the circle: {:.4f}", salesman_circle[0].GetLoss());
    fmt::print("\nSmallest distance in the square: {:.4f}", salesman_square[0].GetLoss());
    fmt::print("\n\n");

    return 0;

}
