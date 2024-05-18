#include <fstream>
#include <cmath>

#include "../../Libraries/fmtlib.h"
#include "../../Libraries/RandomGen/random.hpp"
#include "../../Libraries/classes.hpp"

using namespace std;


int main(int argc, char** argv){

    fmt::print("\n");

    // initial check
    cout << endl;
    if (argc != 3) {
        fmt::print("ERROR! Program usage: {} <n_extractions> <n_blocks>\n\n", argv[0]);
        return -1;
    }

    unsigned int extr = std::stod(argv[1]);     // number of extractions
    unsigned int blocks = stod(argv[2]);        // number of blocks

    if (extr%blocks != 0) {
        fmt::print("ERROR! <n_extractions> must be a mutiple of <n_blocks>\n\n");
        return -1;
    }
    // unsigned int dim = extr/blocks;

    // declaring integrand 
    cosine f(M_PI/2., M_PI/2., 0);
    

    /*––––––––––––––––––––––––––––––– UNIFORM SAMPLING –––––––––––––––––––––––––––––––*/

    // probability distribution for the sampling and inverse of its cumulative
    double (*unif_sampl)(double) = Identity_Distrib;
    double (*unif_inv)(double) = Identity_Inverse;

    BA_Integral int_unif(f, unif_sampl, unif_inv, extr, blocks);
    ofstream out("unif_sampl.csv");
    out << "blocks,extractions,integral,error" << endl;
    int_unif.Progressive(out);
    out.close();


    /*––––––––––––––––––––––––––––––– IMPORTANCE SAMPLING –––––––––––––––––––––––––––––––*/

    double (*impo_sampl)(double) = Importance_Distrib;
    double (*impo_inv)(double) = Importance_Inverse;

    BA_Integral int_impo(f, impo_sampl, impo_inv, extr, blocks);
    out.open("impo_sampl.csv");
    out << "blocks,extractions,integral,error" << endl;
    int_impo.Progressive(out);
    out.close();


    fmt::print("\n");
    
    return 0;
}