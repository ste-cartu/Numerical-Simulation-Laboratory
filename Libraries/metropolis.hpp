# pragma once

#include <armadillo>

#include "RandomGen/random.hpp"
#include "functions.hpp"
#include "fmtlib.hpp"


/*––––––––––––––––––––––––––––––––––– M(RT)^2 ALGORITHM –––––––––––––––––––––––––––––––––––*/


class Metropolis {

    public:
        Metropolis(double step, int type, function& prob, arma::vec x = arma::zeros<arma::vec>(3)) : 
        rnd_("../../Libraries/RandomGen/"), prob_(prob) {
            step_ = step;
            type_ = type;
            dim_ = x.n_elem;

            x_.resize(dim_);
            x_ = x;
            xold_ = x;
        }
        ~Metropolis() {;}

        void Reset(double step, int type, arma::vec x = arma::zeros<arma::vec>(3));         // resets the sampler
        //double SetAcceptance(double target, double prec, double step, int nstep);
        double SetAcceptance(double target, double prec, double step, int nstep, std::ofstream* out = nullptr);     // sets step_ value to have target acceptance rate acc_
        void Propose();                     // proposes a move     
        bool Accept();                      // accepts or rejects the proposed move
        double Distance();                  // returns the distance of the sampled point from the origin
        void Equilibrate(int nstep);        // lets the sampler evolve nstep times to equilibrate the sampling
        void SetSampling(int type) {type_ = type;}
        int GetSampling() {return type_;}
        arma::vec GetPoint() {return x_;}
        double GetAcceptance() {return acc_;}

    protected:
        Random rnd_;                // random numbers generator
        double step_, acc_;         // step length and acceptance rate
        int dim_, type_;            // dimensionality of the problem and type of transition probability
        arma::vec x_, xold_;        // current and old sampled points
        function& prob_;            // probability density to sample

};