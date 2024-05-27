#include <fstream>

#include "../../Libraries/fmtlib.hpp"
#include "../../Libraries/RandomGen/random.hpp"
#include "../../Libraries/blockaverage.hpp"



int main(){

    fmt::print("\n");

    unsigned int extr = 1e6;            // number of extractions
    unsigned int blocks = 1e2;          // number of blocks

    if (extr%blocks != 0) {
        fmt::print("ERROR! <n_extractions> must be a mutiple of <n_blocks>\n\n");
        return -1;
    }

    //unsigned int dim = extr/blocks;     // length of each block
    double S0 = 100;                    // asset price at t=0
    double T = 1;                       // delivery time
    double K = 100;                     // strike price
    double r = 0.1;                     // risk-free interest rate
    double sigma = 0.25;                // volatility


    /*––––––––––––––––––––––––––––––– DIRECT SAMPLING –––––––––––––––––––––––––––––––*/

    // call option
    double (*call_profit)(double, double) = Call_Profit;
    BA_Option call(extr, blocks, call_profit, S0, K, r, sigma);
    std::ofstream out_call("direct_call.csv");
    out_call << "blocks,extractions,profit,error" << std::endl;
    call.Progressive(out_call);
    out_call.close();

    // put option
    double (*put_profit)(double, double) = Put_Profit;
    BA_Option put(extr, blocks, put_profit, S0, K, r, sigma);
    std::ofstream out_put("direct_put.csv");
    out_put << "blocks,extractions,profit,error" << std::endl;
    put.Progressive(out_put);
    out_put.close();


    /*––––––––––––––––––––––––––––––– INDIRECT SAMPLING –––––––––––––––––––––––––––––––*/

    double T0 = 0;
    unsigned int nstep = 100;
    
    // call option
    call.SetSampling(T0, T, nstep);
    out_call.open("indirect_call.csv");
    out_call << "blocks,extractions,profit,error" << std::endl;
    call.Progressive(out_call);
    out_call.close();

    // put option
    put.SetSampling(T0, T, nstep);
    out_put.open("indirect_put.csv");
    out_put << "blocks,extractions,profit,error" << std::endl;
    put.Progressive(out_put);
    out_put.close();




    fmt::print("\n");
    return 0;
}