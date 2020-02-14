// header file for misc_math functions
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#ifndef MISC_MATH_H
#define MISC_MATH_H

#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <random>

double logMultPMF(const std::vector <double>& k, const std::vector <double>& a);
// double logDirichMultPMF(const std::vector <double>& k, const std::vector <double>& a);
double ppoMendelian(const std::vector <int>& p1, const std::vector <int>& p2, const std::vector <int>& o);
int sampleC(const std::vector <int>& items, const std::vector <double>& probs, std::mt19937& rNum);
double varianceEstim(double n, double ss, double mean);

#endif
