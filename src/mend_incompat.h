// header file for mend_incompat functions

#ifndef MEND_INCOMPAT_H
#define MEND_INCOMPAT_H

#include <Rcpp.h>
#include <vector>
#include "utils.h"

void calcProbMIperLocus(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC,
                        const std::vector< std::vector < std::vector <double> > >& genotypeErrorRatesC,
                        const std::vector <std::vector <std::vector <std::vector <double> > > >& lGenos_ssGP,
                        std::vector <double>& pMI);
void calcProbSumMI(const std::vector <double>& pMI, std::vector <double>& pTotalMI);

#endif
