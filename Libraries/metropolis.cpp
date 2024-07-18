#include "metropolis.hpp"
#include <fstream>


/*––––––––––––––––––––––––––––––– METROPOLIS ALGORITHM –––––––––––––––––––––––––––––––*/


void Metropolis :: Reset(double step, int type, arma::vec x) {
    type_ = type;
    step_ = step;
    x_ = x;
    xold_ = x;
}


/* double Metropolis :: SetAcceptance(double target, double prec, double step, int nstep) {
    
    int accepted, counter = 1;
    double corr = step;
    std::string output;

    do {
        accepted = 0;
        step_ = step;
        x_ = arma::zeros<arma::vec>(dim_);
        xold_ = x_;

        for(int j=0 ; j<nstep ; j++) {
            this->Propose();
            if(this->Accept()) accepted++;
        }
        acc_ = (double)accepted/nstep;
        counter++;

        if(corr > prec) corr /= 2.;
        else corr = prec;

        if(acc_ < target - prec) step -= corr;
        else if(acc_ > target + prec) step += corr;

    } while((acc_ < target - prec) or (acc_ > target + prec));

    return step_;
} */


double Metropolis :: SetAcceptance(double target, double prec, double step, int nstep, std::ofstream* out) {
    
    int accepted, counter = 1;
    double corr = step;
    std::string output;

    do {
        accepted = 0;
        step_ = step;
        x_ = arma::zeros<arma::vec>(dim_);
        xold_ = x_;

        for(int j=0 ; j<nstep ; j++) {
            this->Propose();
            if(this->Accept()) accepted++;
        }
        acc_ = (double)accepted/nstep;
        if(out) {
            output = fmt::format("{}) Step: {:.5f}\tAcceptance: {:.3f}\n", counter, step, acc_);
            (*out) << output;
        }
        counter++;

        if(corr > prec) corr /= 2.;
        else corr = prec;

        if(acc_ < target - prec) step -= corr;
        else if(acc_ > target + prec) step += corr;

    } while((acc_ < target - prec) or (acc_ > target + prec));

    return step_;
}


void Metropolis :: Propose() {
    if(type_ == 0)
        for(int i=0 ; i<dim_ ; i++) x_[i] += rnd_.Rannyu(-step_, step_);
    else if(type_ == 1)
        for(int i=0 ; i<dim_ ; i++) x_[i] += rnd_.Gauss(0, step_);
    else {
        fmt::print("ERROR! Unknown type of probability distribution for the proposed step!\n\n");
        exit(0);
    }
}


bool Metropolis :: Accept() {
    double alpha = std::min(1., prob_(x_)/prob_(xold_));
    double rand = rnd_.Rannyu();

    if(rand <= alpha) {
        xold_ = x_;
        return true;
    }
    else {
        x_ = xold_;
        return false;
    }
}


double Metropolis :: Distance() {return arma::norm(x_,2);}


void Metropolis :: Equilibrate(int nstep) {
    for(int i=0 ; i<nstep ; i++) {
        this->Propose();
        this->Accept();
    }
}