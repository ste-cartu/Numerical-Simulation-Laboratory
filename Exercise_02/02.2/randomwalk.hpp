#pragma once

#include "../../Libraries/blockaverage.hpp"


/*––––––––––––––––––––––––––––––– RANDOM WALK –––––––––––––––––––––––––––––––*/


class Walker {

    public:
        Walker(double step=1, double x=0, double y=0, double z=0) : rnd_("../../Libraries/RandomGen/") {
            step_ = step;
            x_ = x;
            y_ = y;
            z_ = z;
        }
        ~Walker() {;}
        void Discrete_Step();
        void Discrete_Step(double r);
        void Continue_Step();
        void Continue_Step(double the, double phi);
        double Distance2(double x0 = 0, double y0 = 0, double z0 = 0);
        double X() {return x_;}
        double Y() {return y_;}
        double Z() {return z_;}

    private:
        double step_;
        double x_ = 0, y_ = 0, z_ = 0;
        Random rnd_;

};