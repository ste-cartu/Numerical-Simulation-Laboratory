#include <iostream>
#include <fstream>

#include "../../Libraries/fmtlib.hpp"
#include "../../Libraries/library.hpp"
#include "../../Libraries/RandomGen/random.hpp"



int main(int argc, char** argv) {

    // initial check
    fmt::print("\n");
    if (argc != 2) {
        fmt::print("ERROR! Program usage: {} <n_extractions>\n\n", argv[0]);
        return -1;
    }

    int len[5] = {0,1,2,10,100};                // lengths of the sums
    unsigned int extr = std::stod(argv[1]);     // number of extractions
    Random rnd("../../Libraries/RandomGen/");   // initializing random numbers generator

    double lambda = 1;          // decay rate for the exponential distribution
    double mu = 0, gamma = 1;   // mean and width for the Cauchy-Lorentz distribution

    double sum_unif = 0, sum_expo = 0, sum_calo = 0;
    double r_unif, r_expo, r_calo;

    // output files
    std::ofstream out_unif("uniform.csv");
    std::ofstream out_expo("exponential.csv");
    std::ofstream out_calo("cauchy-lorentz.csv");
    out_unif << "sums length,extraction,value\n";
    out_expo << "sums length,extraction,value\n";
    out_calo << "sums length,extraction,value\n";

    for (int i=0 ; i<extr ; i++) {

        for (int l=1 ; l<5 ; l++) {

            for (int j=0 ; j<len[l]-len[l-1] ; j++) {
                r_unif = rnd.Rannyu();
                r_expo = Exp_Distrib(r_unif, lambda);
                r_calo = CauLor_Distrib(r_unif, mu, gamma);

                sum_unif += r_unif;
                sum_expo += r_expo;
                sum_calo += r_calo;
            }

            out_unif << len[l] << "," << i+1 << "," << sum_unif/len[l] << std::endl;
            out_expo << len[l] << "," << i+1 << "," << sum_expo/len[l] << std::endl;
            out_calo << len[l] << "," << i+1 << "," << sum_calo/len[l] << std::endl;            

        }

        sum_unif = 0;
        sum_expo = 0;
        sum_calo = 0;

    }


    return 0;


}