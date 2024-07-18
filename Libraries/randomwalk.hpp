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
        void Discrete_Step();
        void Discrete_Step(double r);
        void Continue_Step();
        void Continue_Step(double the, double phi);
        double Distance2(double x0 = 0, double y0 = 0, double z0 = 0);
        double X() {return pos_[0];}
        double Y() {return pos_[1];}
        double Z() {return pos_[2];}

    private:
        double step_;
        int dim_;
        arma::vec pos_;
        arma::vec oldpos_;
        Random rnd_;

};