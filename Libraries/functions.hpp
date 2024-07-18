#pragma once

#include <cmath>
#include <armadillo>

#include "library.hpp"
#include "fmtlib.hpp"



/*––––––––––––––––––––––––––––––– FUNCTIONS –––––––––––––––––––––––––––––––*/

class function {

    public:
        function() {;}
        virtual double Eval(arma::vec) const = 0;
        double operator()(arma::vec x) {return Eval(x);}
        virtual ~function() {;}
        
};


// y = a*cos(b*x) + c
class cosine : public function {

    public:
        cosine() {a_=1; b_=1.; c_=0.;}
        cosine(double a, double b, double c) {a_=a; b_=b; c_=c;}
        virtual ~cosine() {;}
        double Eval(arma::vec x) const override {return a_*cos(b_*x[0]) + c_;}
        void SetA(double a) {a_=a;}
        double GetA() const {return a_;}
        void SetB(double b) {b_=b;}
        double GetB() const {return b_;}
        void SetC(double c) {c_=c;};
        double GetC() const {return c_;}

    private:
        double a_, b_, c_;

};


// ground state wave function of the hydrogen atom with Bohr radius (a0) = 1
class hydrogen_100 : public function {

    public:
        hydrogen_100() {;}
        virtual ~hydrogen_100() {;}
        double Eval(arma::vec r) const override {
            return (1./M_PI) * std::exp(-2. * arma::norm(r,2));
        }
        
};


// first excited state wave function of the hydrogen atom with Bohr radius (a0) = 1
class hydrogen_210 : public function {

    public:
        hydrogen_210() {;}
        virtual ~hydrogen_210() {;}
        double Eval(arma::vec v) const override {
            double r = arma::norm(v,2);
            double theta = 2.*std::atan(std::sqrt(v[0]*v[0] + v[1]*v[1]) / (v[3] + r));
            return (r*r)/(32.*M_PI) * std::exp(-r) * std::pow(std::cos(theta), 2);
        }

};


class variational : virtual public function {

    public:
        variational() {;}
        virtual ~variational() {;}
        void SetParameters(double mu, double sigma) {mu_ = mu, sigma_ = sigma;}
    
    protected:
        double mu_, sigma_;
        double min_(arma::vec x) const {return std::pow((x[0] - mu_), 2) / (2.*sigma_*sigma_);}
        double plu_(arma::vec x) const {return std::pow((x[0] + mu_), 2) / (2.*sigma_*sigma_);}

};


class wave_func : public variational {

    public:
        wave_func(double mu, double sigma) {mu_ = mu, sigma_ = sigma;}
        virtual ~wave_func() {;}
        double Eval(arma::vec x) const override {return std::exp(-min_(x)) + std::exp(-plu_(x));}

};


class wave_func2 : public variational {

    public:
        wave_func2(double mu, double sigma) {mu_ = mu, sigma_ = sigma;}
        virtual ~wave_func2() {;}
        double Eval(arma::vec x) const override {return std::pow(std::exp(-min_(x)) + std::exp(-plu_(x)), 2);}

};


class hamiltonian : public variational {

    public:
        hamiltonian(double mu, double sigma) {mu_ = mu, sigma_ = sigma;}
        virtual ~hamiltonian() {;}
        double Eval(arma::vec x) const override {
            //double psi = std::exp(-min_(x)) + std::exp(-plu_(x));
            wave_func psi(mu_, sigma_);

            double kin = - (1./(2.*sigma_*sigma_*psi(x))) * ( (2.*min_(x) - 1.)*std::exp(-min_(x)) + (2.*plu_(x) - 1.)*std::exp(-plu_(x)) );
            double pot = std::pow(x[0], 4) - (5./2.)*std::pow(x[0], 2);

            return kin + pot;
        }

};