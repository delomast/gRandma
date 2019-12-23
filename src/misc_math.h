// header file for misc_math functions

#ifndef MISC_MATH_H
#define MISC_MATH_H

#include <Rcpp.h>
#include <math.h>
#include <vector>

double logDirichMultPMF(std::vector <double> k, std::vector <double> a);
std::vector <std::vector <int> > listAllPairs(std::vector <int> allInds);

#endif
