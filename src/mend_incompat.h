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
void calcProbMIperLocus_obs(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC,
                            const std::vector <std::vector <std::vector <std::vector <double> > > >& lGenos_trio_obs,
                            std::vector <double>& pMI);
void calcProbSumMI(const std::vector <double>& pMI, std::vector <double>& pTotalMI);
void calcProbSumMI_returnAll(const std::vector <double>& pMI, std::vector <double>& pTotalMI,
                             std::vector <std::vector <double> >& all_pTotalMI);
void fwd_ssGP(const int miLimit, const int maxMissingGenos, const int nLoci, 
              const std::vector <double>& prob_missing_geno, 
              const std::vector <double>& pMI_noMiss, 
              std::vector <std::vector <std::vector <std::vector <std::vector <double> > > > >& all_states);
void fwd_sP(const int miLimit, const int maxMissingGenos, const int nLoci, 
            const std::vector <double>& prob_missing_geno, 
            const std::vector <double>& pMI_noMiss, 
            std::vector <std::vector <std::vector <std::vector <double> > > >& all_states);
void calcProbMIperLocus_sP(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC,
                        const std::vector <std::vector < std::vector <double> > >& genotypeErrorRatesC,
                        const std::vector <std::vector <std::vector <double> > >& lGenos_sP,
                        std::vector <double>& pMI);
void calcProbMIperLocus_Unrel(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC,
                        const std::vector< std::vector < std::vector <double> > >& genotypeErrorRatesC,
                        const std::vector <std::vector <double> >& lGenos_one,
                        const std::vector <std::vector <double> >& lGenos_two,
                        std::vector <double>& pMI);
void strat_lGenosBuilder_sP(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC,
                        const std::vector< std::vector < std::vector <double> > >& genotypeErrorRatesC,
                        const std::vector <std::vector <double> >& lGenos_base,
                        const std::vector <std::vector <double> >& lGenos_randomDescendant,
                        std::vector <std::vector <std::vector <double> > >& lGenos_OBS_sample_MI,
								std::vector <std::vector <std::vector <double> > >& lGenos_OBS_sample_no_MI);
void calcProbMIperLocus_Pair(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC,
                        const std::vector< std::vector < std::vector <double> > >& genotypeErrorRatesC,
                        const std::vector <std::vector <std::vector <double> > >& lGenos_tRel,
                        std::vector <double>& pMI);
void strat_lGenosBuilder_Pairwise(const std::vector <std::vector <std::vector <int> > >& genotypeKeyC,
                        const std::vector< std::vector < std::vector <double> > >& genotypeErrorRatesC,
                        const std::vector <std::vector <std::vector <double> > >& lGenos_tRel_OBS,
                        std::vector <std::vector <std::vector <double> > >& lGenos_OBS_sample_MI,
								std::vector <std::vector <std::vector <double> > >& lGenos_OBS_sample_no_MI);

#endif
