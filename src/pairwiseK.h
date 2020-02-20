// header file for misc_math functions
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#ifndef PAIRWISEK_H
#define PAIRWISEK_H

#include <Rcpp.h>
#include <math.h>
#include <vector>
#include "misc_math.h"

void pairwiseK(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC, 
	                      const  std::vector <std::vector <double> >& lGenos_base, 
	                      const std::vector <std::vector <std::vector <double> > >& baselineParamsC,
	                      const std::vector <std::vector <std::vector <double> > >& unsampledPopParamsC,
	                      const int pop,
	                      const std::vector <double>& k_prob,
	                      std::vector <std::vector <std::vector <double> > >& lGenos_k);
void create_pairswise_OBSvector(const std::vector <std::vector <std::vector <double> > >& lGenos_k, 
                          const std::vector <std::vector <std::vector <double> > >& genotypeErrorRatesC, 
                          std::vector <std::vector <std::vector <double> > >& lGenos_k_OBS
                          );

#endif
