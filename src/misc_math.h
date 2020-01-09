// header file for misc_math functions

#ifndef MISC_MATH_H
#define MISC_MATH_H

#include <Rcpp.h>
#include <math.h>
#include <vector>

double logDirichMultPMF(const std::vector <double>& k, const std::vector <double>& a);
double ppoMendelian(const std::vector <int>& p1, const std::vector <int>& p2, const std::vector <int>& o);

#endif
