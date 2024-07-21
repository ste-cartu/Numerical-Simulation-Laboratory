#include "randomwalk.hpp"


/*––––––––––––––––––––––––––––––– RANDOM WALK –––––––––––––––––––––––––––––––*/


void Walker :: Discrete_Step() {
    double r = rnd_.Rannyu(-dim_, dim_);
    for(int i=0 ; i<dim_ ; i++) {
        if(r >= -dim_ + i and r < -dim_ + i + 1) pos_[i] -= step_;
        else if(r >= i and r < i + 1) pos_[i] += step_;
    }
}


void Walker :: Discrete_Step(double r) {
    if (r < 0 || r >= 1) {
        fmt::print("ERROR! Random number must be in range [0,1)!\n\n");
        exit(0);
    }
    r = 2*dim_*r - dim_;
    for(int i=0 ; i<dim_ ; i++) {
        if(r >= -dim_ + i and r < -dim_ + i + 1) pos_[i] -= step_;
        else if(r >= i and r < i + 1) pos_[i] += step_;
    }
}


void Walker :: Continue_Step() {
    if(dim_ != 3) {
        fmt::print("ERROR! Spherical coordinate step must be in 3D!\n\n");
        exit(0);
    }

    double theta = std::acos(1 - 2*rnd_.Rannyu());
    double phi = rnd_.Rannyu(0, 2*M_PI);

    pos_[0] += step_ * std::sin(theta) * std::cos(phi);
    pos_[1] += step_ * std::sin(theta) * std::sin(phi);
    pos_[2] += step_ * std::cos(theta);
}


void Walker :: Continue_Step(double the, double phi) {
    if(dim_ != 3) {
        fmt::print("ERROR! Spherical coordinate step must be in 3D!\n\n");
        exit(0);
    }
    if (the < 0 || the >= 1 || phi < 0 || phi >= 1 ) {
        fmt::print("ERROR! Random number must be in range [0,1)!\n\n");
        exit(0);
    }
    double theta = std::acos(1 - 2*the);
    phi *= 2*M_PI;

    pos_[0] += step_ * std::sin(theta) * std::cos(phi);
    pos_[1] += step_ * std::sin(theta) * std::sin(phi);
    pos_[2] += step_ * std::cos(theta);
}


double Walker :: Distance2(double x0, double y0, double z0) {
    return pow(arma::norm(pos_,2), 2);
}
