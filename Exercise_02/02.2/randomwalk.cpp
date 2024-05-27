#include "randomwalk.hpp"


/*––––––––––––––––––––––––––––––– RANDOM WALK –––––––––––––––––––––––––––––––*/

void Walker :: Discrete_Step() {
    double r = rnd_.Rannyu(-3, 3);
    if (r >= -3 && r < -2) {x_ -= step_;}
    else if (r >= -2 && r < -1) {y_ -= step_;}
    else if (r >= -1 && r < 0) {z_ -= step_;}
    else if (r >= 0 && r < 1) {x_ += step_;}
    else if (r >= 1 && r < 2) {y_ += step_;}
    else if (r >= 2 && r < 3) {z_ += step_;}
}


void Walker :: Discrete_Step(double r) {
    if (r < 0 || r >= 1) {
        fmt::print("ERROR! Random number must be in range [0,1)!\n\n");
        exit(0);
    }
    r = 6*r - 3;

    if (r >= -3 && r < -2) {x_ -= step_;}
    else if (r >= -2 && r < -1) {y_ -= step_;}
    else if (r >= -1 && r < 0) {z_ -= step_;}
    else if (r >= 0 && r < 1) {x_ += step_;}
    else if (r >= 1 && r < 2) {y_ += step_;}
    else if (r >= 2 && r < 3) {z_ += step_;}
}


void Walker :: Continue_Step() {
    double theta = std::acos(1 - 2*rnd_.Rannyu());
    double phi = rnd_.Rannyu(0, 2*M_PI);

    x_ += step_ * std::sin(theta) * std::cos(phi);
    y_ += step_ * std::sin(theta) * std::sin(phi);
    z_ += step_ * std::cos(theta);
}


void Walker :: Continue_Step(double the, double phi) {
    double theta = std::acos(1 - 2*the);
    phi *= 2*M_PI;

    x_ += step_ * std::sin(theta) * std::cos(phi);
    y_ += step_ * std::sin(theta) * std::sin(phi);
    z_ += step_ * std::cos(theta);
}


double Walker :: Distance2(double x0, double y0, double z0) {
    return std::pow((x_ - x0), 2) + std::pow((y_ - y0), 2) + std::pow((z_ - z0), 2);
}