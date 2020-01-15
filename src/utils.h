// header file for utils functions

#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <vector>
#include "misc_math.h"

void listAllPairs(const std::vector <int>& allInds,
                                              std::vector <std::vector <int> >& pairs);
std::vector <std::vector <int> > rcppMatrixToVectorInt(Rcpp::NumericMatrix inputMatrix);
std::vector <std::vector <double> > rcppListToVectorDouble(Rcpp::List inputList);
std::vector <std::vector <std::vector <int> > > convertGenotypeKey(Rcpp::List inputList);
std::vector < std::vector <double> > rcppMatrixToVectorDouble(Rcpp::NumericMatrix inputMatrix);
void createPPOvector(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC,
                std::vector <std::vector <std::vector <std::vector <double> > > >& lGenos_ppo);
bool noAllelesInCommonGP(const std::vector <int>& gp1, const std::vector <int>& gp2, 
                         const std::vector <int>& d);
void createSSGPvector(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC, 
                      const  std::vector <std::vector <double> >& lGenos_base, 
                      const std::vector <std::vector <std::vector <double> > >& unsampledPopParamsC, 
                      const int pop, 
                      std::vector <std::vector <std::vector <std::vector <double> > > >& lGenos_ssGP);
void createOBSvector(const std::vector <std::vector <std::vector <std::vector <double> > > >& lGenos_ssGP, 
                          const std::vector< std::vector < std::vector <double> > >& genotypeErrorRatesC, 
                          std::vector <std::vector <std::vector <std::vector <double> > > >& lGenos_ssGP_OBS,
                          const std::vector <std::vector <double> >& lGenos_base,
                          const std::vector <std::vector <double> >& lGenos_unsamp,
                          std::vector <std::vector <std::vector <std::vector <double> > > >& lGenos_Unrelated_OBS);
void create_CORR_OBSvector(const std::vector< std::vector < std::vector <double> > >& genotypeErrorRatesC, 
                           const std::vector <std::vector <double> >& lGenos_randomDescendant,
                          const std::vector <std::vector <double> >& lGenos_base,
                          std::vector <std::vector <std::vector <std::vector <double> > > >& CORR_lGenos_OBS
                          );

#endif
