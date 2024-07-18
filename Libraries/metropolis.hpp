# pragma once

#include <armadillo>

#include "RandomGen/random.hpp"
#include "functions.hpp"
#include "fmtlib.hpp"


/*––––––––––––––––––––––––––––––– METROPOLIS ALGORITHM –––––––––––––––––––––––––––––––*/


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
        void Reset(double step, int type, arma::vec x = arma::zeros<arma::vec>(3));
        //double SetAcceptance(double target, double prec, double step, int nstep);
        double SetAcceptance(double target, double prec, double step, int nstep, std::ofstream* out = nullptr);
        void Propose();
        bool Accept();
        double Distance();
        void SetSampling(int type) {type_ = type;}
        int GetSampling() {return type_;}
        arma::vec GetPoint() {return x_;}
        double GetAcceptance() {return acc_;}
        void Equilibrate(int nstep);

    protected:
        Random rnd_;
        double step_, acc_;
        int dim_, type_;
        arma::vec x_, xold_;
        function& prob_;

};