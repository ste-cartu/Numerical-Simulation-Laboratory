#include <fstream>

#include "../../Libraries/functions.hpp"
#include "../../Libraries/metropolis.hpp"
#include "../../Libraries/blockaverage.hpp"

int main() {

    fmt::print("\n");
    /*––––––––––––––––––––––––––––––––– FIRST EVALUATION –––––––––––––––––––––––––––––––––*/

    // blocking averagre variables
    unsigned int extr = 1e6;
    unsigned int blocks = 1e2;

    // initializing the integrand with parameters sigma and mu
    double mu = 0., sigma = 1.;
    wave_func2 psi(mu, sigma);
    hamiltonian H(mu, sigma);

    // Metropolis sampler
    double step = 1.;
    int unif = 0;
    arma::vec x = arma::zeros<arma::vec>(1);
    Metropolis metro_1(step, unif, psi, x);

    // setting acceptance at 50%
    double target = 0.5, prec = 0.01, nstep = 1e3;
    std::string filename = fmt::format("acceptance_m={:.1f}_s={:.1f}.txt", mu, sigma);
    std::ofstream acc(filename);
    step = metro_1.SetAcceptance(target, prec, step, nstep, &acc);
    acc.close();

    // blocking average
    BA_SimAnn energy_1(extr, blocks, metro_1, H);
    filename = fmt::format("energy_m={:.1f}_s={:.1f}.csv", mu, sigma);
    std::ofstream out(filename);
    out << "blocks,extractions,energy,error" << std::endl;
    energy_1.Progressive(out);
    out.close();
    std::string newname = fmt::format("sampling_m={:.1f}_s={:.1f}.csv", mu, sigma);
    rename("sampling.csv", newname.c_str());
    fmt::print("\n");


    /*––––––––––––––––––––––––––––––––– OTHER PARAMETERS VALUES –––––––––––––––––––––––––––––––––*/
    mu = 1., sigma = 0.5;
    psi.SetParameters(mu, sigma);
    H.SetParameters(mu, sigma);

    // building a new Metropolis sampler
    Metropolis metro_2(step, unif, psi, x);

    // setting acceptance at 50%
    step = 1.;
    filename = fmt::format("acceptance_m={:.1f}_s={:.1f}.txt", mu, sigma);
    acc.open(filename);
    step = metro_2.SetAcceptance(target, prec, step, nstep, &acc);
    acc.close();

    // blocking average
    BA_SimAnn energy_2(extr, blocks, metro_2, H);
    filename = fmt::format("energy_m={:.1f}_s={:.1f}.csv", mu, sigma);
    out.open(filename);
    out << "blocks,extractions,energy,error" << std::endl;
    energy_2.Progressive(out);
    out.close();
    newname = fmt::format("sampling_m={:.1f}_s={:.1f}.csv", mu, sigma).c_str();
    rename("sampling.csv", newname.c_str());


    fmt::print("\n\n");

    return 0;
}