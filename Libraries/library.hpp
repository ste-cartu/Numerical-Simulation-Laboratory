#pragma once

#include <cmath>
#include <algorithm>


// error: mean standard deviation: sqrt[(<x^2> - <x>^2) / (n-1)]
// used in: 01.1
double Error(double x, double x2, int n);


// obtaining exponential ditribution from a random variable (x) uniformly distributed
// used in: 01.2
double Exp_Distrib(double x, double lambda);


// obtaining Cauchy-Lorentz ditribution from a random variable (x) uniformly distributed
// used in: 01.2
double CauLor_Distrib(double x, double mu, double gamma);


// identity function
// used in: 02.1
double Identity_Distrib(double x);


// inverse identity
// used in: 02.1
double Identity_Inverse(double x);


// probability distribution for importance sampling
// used in: 02.1
double Importance_Distrib(double x);


// inverse of the cumulative
// used in: 02.1
double Importance_Inverse(double x);


// call option profit
// used in 03.1
double Call_Profit(double S, double K);


// put option profit
// used in 03.1
double Put_Profit(double S, double K);