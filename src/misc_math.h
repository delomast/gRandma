// header file for misc_math functions

#ifndef MISC_MATH_H
#define MISC_MATH_H

#include <Rcpp.h>
#include <math.h>
#include <vector>

double logDirichMultPMF(std::vector <double> k, std::vector <double> a);
std::vector <std::vector <int> > listAllPairs(std::vector <int> allInds);
std::vector <std::vector <int> > rcppMatrixToVectorInt(Rcpp::NumericMatrix inputMatrix);
std::vector <std::vector <double> > rcppListToVectorDouble(Rcpp::List inputList);
std::vector <std::vector <std::vector <int> > > convertGenotypeKey(Rcpp::List inputList);
double ppoMendelian(std::vector <int> p1, std::vector <int> p2, std::vector <int> o);

#endif
