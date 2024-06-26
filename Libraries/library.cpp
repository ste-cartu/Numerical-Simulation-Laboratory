#include "library.hpp"


double Error(double x, double x2, int n){
    if (n == 0) {return 0;}
    else {return std::sqrt((x2 - x*x) / n);}
}


double Exp_Distrib(double x, double lambda) {return -(1./lambda)*std::log(1 - x);}


double CauLor_Distrib(double x, double mu, double gamma) {return gamma*std::tan(M_PI*(x - 0.5)) - mu;}


double Identity_Distrib(double x) {return 2*x;}


double Identity_Inverse(double x) {return std::sqrt(x);}


double Importance_Distrib(double x) {return 2*(1 - x);}


double Importance_Inverse(double x) {return std::sqrt(1-x) + 1;}


double Call_Profit(double S, double K) {return std::max(0., S - K);}


double Put_Profit(double S, double K) {return std::max(0., K - S);}


void Progress_Bar(int progress, int total, int bar_width) {

    float percentage = static_cast<float>(progress) / total;
    int pos = static_cast<int>(bar_width * percentage);

    std::string bar;
    for (int i = 0; i < bar_width; ++i) {
        if (i <= pos) bar += "◼︎";
        else bar += " ";
    }

    fmt::print("|{}| {:3d} %\r", bar, int(percentage * 100.0));
    std::fflush(stdout);
}