#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <armadillo>

#include "RandomGen/random.hpp"
#include "library.hpp"
#include "fmtlib.hpp"
#include "functions.hpp"
#include "metropolis.hpp"


/*––––––––––––––––––––––––––––––– BLOCKING AVERAGE –––––––––––––––––––––––––––––––*/

// Blocking Average method (BA)
class BlockingAverage{

    public:
        BlockingAverage(unsigned int extr, unsigned int blocks) : rnd_("../../Libraries/RandomGen/") {
            if (extr%blocks != 0) {
                std::cout << "ERROR! <n_extractions> must be a mutiple of <n_blocks>" << std::endl << std::endl;
                exit(0);
            }
            extr_ = extr;
            blocks_ = blocks;
            dim_ = extr_/blocks_;
        }
        ~BlockingAverage() {;}
        void Progressive(std::ofstream&, bool progress = true);
        virtual double Increase() = 0;
        double GetValue() {return average_;}
        double GetError() {return error_;}

    protected:
        Random rnd_;
        unsigned int extr_, blocks_, dim_;
        double average_, error_;

};


// BA: expected value of a uniform distribution 
// used in: 01.1
class BA_Mean : public BlockingAverage {

    public:
        BA_Mean(unsigned int extr, unsigned int blocks) : BlockingAverage(extr, blocks) {;}
        virtual ~BA_Mean() {;}
        double Increase() override;

};


// BA: variance of a uniform distribution 
// used in: 01.2
class BA_Variance : public BlockingAverage {

    public:
        BA_Variance(unsigned int extr, unsigned int blocks) : BlockingAverage(extr, blocks) {;}
        virtual ~BA_Variance() {;}
        double Increase() override;

};


// BA: calculation of the value of pi with Buffon experiment simulation
// used in: 01.3
class BA_Buffon : public BlockingAverage {

    public:
        BA_Buffon(double length, double distance, unsigned int extr, unsigned int blocks) : BlockingAverage(extr, blocks) {
            length_ = length;
            distance_ = distance;
        }
        virtual ~BA_Buffon() {;}
        double Increase() override;

    private:
        double length_, distance_;

};


// BA: Monte Carlo integration with importance sampling
// used in: 02.1
class BA_Integral : public BlockingAverage {

    public:
        BA_Integral(function& f, double (*sampling)(double), double (*inverse)(double), unsigned int extr, unsigned int blocks) : 
        BlockingAverage(extr, blocks), f_(f), Sampl_(sampling), Inv_(inverse) {;}
        virtual ~BA_Integral() {;}
        double Increase() override;

    private:
        function& f_;
        double (*Sampl_)(double);
        double (*Inv_)(double);

};


// BA: price of an european option
// used in: 03.1
class BA_Option : public BlockingAverage {

    public:
        BA_Option(unsigned int extr, unsigned int blocks, double (*profit)(double, double), double S0, double K, double mu, double sigma) : 
        BlockingAverage(extr, blocks), Profit_(profit) {
            initial_ = S0;
            strike_ = K;
            drift_ = mu;
            volatility_ = sigma;
        }
        virtual ~BA_Option() {;}
        void SetSampling(double tinit, double tfin, unsigned int nstep);
        double Increase() override;

    private:
        double (*Profit_)(double, double);
        double initial_, strike_;
        double drift_, volatility_;
        double start_ = 0., stop_ = 1.;
        unsigned int nstep_ = 1;
        
};


// BA: distance of a set of points sampled using M(RT)^2 algorithm
// used in 05.1
class BA_Metro : public BlockingAverage {

    public:
        BA_Metro(unsigned int extr, unsigned int blocks, Metropolis& metro, bool save = true) : 
        BlockingAverage(extr, blocks), metro_(metro) {;}
        virtual ~BA_Metro() {;}
        //BA_Metro(unsigned int extr, unsigned int blocks, std::shared_ptr<Metropolis>& metro, bool save) : 
        double Increase() override;
        void Equilib(unsigned int nstep) {metro_.Equilibrate(nstep);}
        void Reset(int type, bool answer);
        void Reset(int type);
        void SetSave(bool answer) {save_ = answer;}

    
    protected:
        //std::shared_ptr<Metropolis> metro_;
        Metropolis& metro_;
        bool save_;

};

// BA: simulated annealing using M(RT)^2 algorithm
// used in 08.1, 08.2
class BA_SimAnn : public BlockingAverage{

    public:
        BA_SimAnn(unsigned int extr, unsigned int blocks, Metropolis& metro, function& hamiltonian, bool save = false) :
        BlockingAverage(extr, blocks), metro_(metro), H_(hamiltonian) {save_ = save;}
        virtual ~BA_SimAnn() {;}
        double Increase() override;
        void SetSave(bool answer) {save_ = answer;}

    protected:
        Metropolis& metro_;
        function& H_;
        bool save_;
        
};








