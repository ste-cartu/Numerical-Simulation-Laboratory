#include <fstream>

#include "../../Libraries/functions.hpp"
#include "../../Libraries/metropolis.hpp"
#include "../../Libraries/blockaverage.hpp"

int main() {

    fmt::print("\n");

    // blocking averagre variables
    unsigned int extr = 1e6;
    unsigned int blocks = 1e2;

    // initializing the integrand with parameters sigma and mu
    double mu = 1., sigma = 1.;
    wave_func2 psi(mu, sigma);
    hamiltonian H(mu, sigma);

    // Metropolis variables
    double step = 1., acc;
    int unif = 0, nstep_acc = 1e3, nstep_eq = 1e4;
    double target_acc = 0.5, prec = 0.01;
    arma::vec x = arma::zeros<arma::vec>(1);

    // initializing variables for simulated annealing
    double T, beta;                 // temperature going from 1 to 0.001
    double T0 = 1.;
    double mu_old = mu, sigma_old = sigma;
    double mu_best, sigma_best;
    double ene = 0., ene_old = 0., ene_best = ene;
    double err = 0., err_best = 0.;
    double alpha, r;
    double delta_mu = 0.5, delta_sigma = 0.5;
    double d_mu, d_sigma;
    int counter = 1;

    std::ofstream opt, sim;
    opt.open("optimization.csv");
    opt << "step,temperature,mu,sigma,acceptance,energy,error" << std::endl;
    opt.close();
    std::string output;

    // random numbers generator
    std::string const rnd_path = "../../Libraries/RandomGen/";
    Random rnd(rnd_path);

    int N = 10000;
    fmt::print("Optimization of the parameters\n");
    for(int i=0 ; i<N ; i++) {
        beta = i+1;
        T = 1./beta;
        Progress_Bar(i, N-1);

        // updating parameters
        d_mu = d_mu > mu_old ? mu_old : delta_mu;
        d_sigma = d_sigma > sigma_old ? sigma_old : delta_sigma;
        mu = std::fabs(mu_old + rnd.Gauss(0.,d_mu*T/T0));
        sigma = std::fabs(sigma_old + rnd.Gauss(0.,d_sigma*T/T0));

        psi.SetParameters(mu, sigma);
        H.SetParameters(mu, sigma);

        // Metropolis sampler
        Metropolis metro(step, unif, psi, x);
        metro.Equilibrate(nstep_eq);
        step = metro.SetAcceptance(target_acc, prec, step, nstep_acc);
        acc = metro.GetAcceptance();

        // blocking average
        BA_SimAnn energy(extr, blocks, metro, H);
        sim.open("energy.csv");
        sim << "blocks,extractions,energy,error" << std::endl;
        energy.Progressive(sim, false);
        sim.close();
        ene = energy.GetValue();
        err = energy.GetError();

        // printing the progress
        opt.open("optimization.csv", std::ios::app);
        output = fmt::format("{},{:.5f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f}\n", counter, T, mu, sigma, acc, ene, err);
        opt << output;
        opt.close();
        counter++;

        // acception/rejection
        alpha = std::min(1., std::exp(-beta * (ene - ene_old)));
        r = rnd.Rannyu();
        if(r < alpha) {
            mu_old = mu;
            sigma_old = sigma;
            if(ene < ene_best) {
                mu_best = mu;
                sigma_best = sigma;
                ene_best = ene;
                err_best = err;
            }
        }
        else {
            mu = mu_old;
            sigma = sigma_old;
        }

        metro.~Metropolis();
        energy.~BA_SimAnn();
    }

    fmt::print("\n\nEstimation of the best energy\n");
    // calculating energy with the best parameters
    psi.SetParameters(mu_best, sigma_best);
    H.SetParameters(mu_best, sigma_best);

    Metropolis metro(step, unif, psi, x);
    metro.Equilibrate(nstep_eq);
    metro.SetAcceptance(target_acc, prec, step, nstep_acc);

    std::ofstream sampl("sampling.csv");
    sampl << "blk_step,position" << std::endl;
    BA_SimAnn energy(extr, blocks, metro, H, true);
    std::string filename = fmt::format("energy_best_mu={:.3f}_sigma={:.3f}.csv", mu_best, sigma_best);
    sim.open(filename);
    sim << "blocks,extractions,energy,error" << std::endl;
    energy.Progressive(sim);
    sim.close();

    fmt::print("\n\nBest energy: {} ± {}", ene_best, err_best);
    fmt::print("\nBest parameters: 𝝁 = {}, 𝝈 = {}", mu_best, sigma_best);
    fmt::print("\n\n");

    return 0;
}