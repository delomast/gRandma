// header file for utils functions

#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <vector>
#include "misc_math.h"

std::vector <std::vector <int> > listAllPairs(std::vector <int> allInds);
std::vector <std::vector <int> > rcppMatrixToVectorInt(Rcpp::NumericMatrix inputMatrix);
std::vector <std::vector <double> > rcppListToVectorDouble(Rcpp::List inputList);
std::vector <std::vector <std::vector <int> > > convertGenotypeKey(Rcpp::List inputList);
std::vector < std::vector <double> > rcppMatrixToVectorDouble(Rcpp::NumericMatrix inputMatrix);
std::vector <std::vector <std::vector <std::vector <double> > > > createPPOvector(std::vector <std::vector <std::vector <int> > > genotypeKeyC);
bool noAllelesInCommonGP(std::vector <int> gp1, std::vector <int> gp2, std::vector <int> d);

#endif
