//#include <iostream>
#include <fstream>

#include "../../Libraries/fmtlib.h"
#include "../../Libraries/library.hpp"
#include "../../Libraries/RandomGen/random.hpp"
#include "../../Libraries/classes.hpp"



int main(int argc, char** argv){

    // initial check
    fmt::print("\n");
    if (argc != 5) {
        fmt::print("ERROR! Program usage: {} <needle_length> <line_distance> <n_throws> <n_blocks>\n\n", argv[0]);
        return -1;
    }

    // setting parameters
    double l = std::stod(argv[1]);                  // length of the needle
    double d = std::stod(argv[2]);                  // distance of the lines

    if (d <= l) {
        fmt::print("ERROR! <line_distance> must be greater than <needle_length>\n\n");
        return -1;
    }

    unsigned int throws = std::stod(argv[3]);       // number of throws
    unsigned int blocks = std::stod(argv[4]);       // number of blocks

    if (throws%blocks != 0) {
        fmt::print("ERROR! <n_throws> must be a mutiple of <n_blocks>\n\n");
        return -1;
    }

    // initializing random numbers generator
    Random rnd("../../Libraries/RandomGen/");

    std::ofstream out("pi.csv");
    out << "blocks,extractions,value of pi,error" << std::endl;

    BA_Buffon buff(l, d, throws, blocks);
    buff.Progressive(out);
    out.close();
    



    return 0;
}