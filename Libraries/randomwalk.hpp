#pragma once

#include "blockaverage.hpp"
#include <armadillo>


/*––––––––––––––––––––––––––––––– RANDOM WALK –––––––––––––––––––––––––––––––*/


class Walker {

    public:
        Walker(double step=1, int dim=3, double x=0, double y=0, double z=0) : rnd_("../../Libraries/RandomGen/") {
            step_ = step;
            dim_ = dim;
            pos_.resize(dim_);
            oldpos_.resize(dim_);

            pos_[0] = x;
            pos_[1] = y;
            pos_[2] = z;

            oldpos_[0] = x;
            oldpos_[1] = y;
            oldpos_[2] = z;
        }
        ~Walker() {;}

        void Discrete_Step();                           // performs a random walk step in a discrete space (i.e. a lattice)
        void Discrete_Step(double r);                   // performs a random walk step in a discrete space (i.e. a lattice)
        void Continue_Step();                           // performs a random walk step in a continue space
        void Continue_Step(double the, double phi);     // performs a random walk step in a continue space
        double Distance2(double x0 = 0, double y0 = 0, double z0 = 0);      // returns the square value of the random walk distance from the origin
        double X() {return pos_[0];}
        double Y() {return pos_[1];}
        double Z() {return pos_[2];}

    private:
        Random rnd_;            // random numbers generator
        double step_;           // length of the random walk steps
        int dim_;               // dimensionality of the random walk
        arma::vec pos_;         // current position of the random walk
        arma::vec oldpos_;      // old position of the random walk

};