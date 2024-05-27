#ifndef classes_hpp
#define classes_hpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "RandomGen/random.hpp"
#include "library.hpp"
#include "fmtlib.h"




/*––––––––––––––––––––––––––––––– FUNCTIONS –––––––––––––––––––––––––––––––*/

class function {

    public:
        function() {m_a = 0;}
        virtual double Eval(double) const = 0;
        double operator()(double x) {return Eval(x);}
        virtual ~function() {;}

    protected:
        double m_a;
    
};


//y = a*cos(b*x) + c
class cosine : public function {

    public:
        cosine() {a_=1; b_=1.; c_=0.;}
        cosine(double a, double b, double c) {a_=a; b_=b; c_=c;}
        ~cosine() {;}
        double Eval (double x) const override {return a_*cos(b_*x) + c_;}
        void SetA (double a) {a_=a;}
        double GetA () const {return a_;}
        void SetB (double b) {b_=b;}
        double GetB () const {return b_;}
        void SetC (double c) {c_=c;};
        double GetC () const {return c_;}

    private:
        double a_, b_, c_;

};




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
        void Progressive(std::ofstream&);
        virtual double Increase() = 0;

    protected:
        Random rnd_;
        unsigned int extr_, blocks_, dim_;

};


// BA: expected value of a uniform distribution 
// used in: 01.1
class BA_Mean : public BlockingAverage {

    public:
        BA_Mean(unsigned int extr, unsigned int blocks) : BlockingAverage(extr, blocks) {;}
        ~BA_Mean() {;}
        double Increase() override;

};


// BA: variance of a uniform distribution 
// used in: 01.2
class BA_Variance : public BlockingAverage {

    public:
        BA_Variance(unsigned int extr, unsigned int blocks) : BlockingAverage(extr, blocks) {;}
        ~BA_Variance() {;}
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
        ~BA_Buffon() {;}
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
        ~BA_Integral() {;}
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
        ~BA_Option() {;}
        void SetSampling(double tinit, double tfin, unsigned int nstep);
        double Increase() override;

    private:

        double (*Profit_)(double, double);
        double initial_, strike_;
        double drift_, volatility_;
        double start_ = 0., stop_ = 1.;
        unsigned int nstep_ = 1;
        
};













#endif