// functions for calculating the threshold number of Mendelian incompatibilities to use

#include <Rcpp.h>
#include <vector>
#include "utils.h"
#include "mend_incompat.h"

using namespace std;

// at each locus, calculate the probability of a MI | true grandparents
void calcProbMIperLocus(const vector <vector <vector <int> > >& genotypeKeyC,
                        const vector< vector < vector <double> > >& genotypeErrorRatesC,
                        const vector <vector <vector <vector <double> > > >& lGenos_ssGP,
                        vector <double>& pMI){
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){
		// summing accross all possible obs and true, calculate
		// P(obs, true, MI | true grandparents) = P(MI|o)P(o|t)P(t|true grandparents)
		for(int gp1 = 0, max3 = genotypeKeyC[i].size(); gp1 < max3; gp1++){ // for all possible OBS genotypes
			for(int gp2 = 0; gp2 < max3; gp2++){ 
				for(int d = 0; d < max3; d++){ 
					//check if observed MI
					if(!noAllelesInCommonGP(
						genotypeKeyC[i][gp1], // gp1 genotype as vector of alleles
						genotypeKeyC[i][gp2], // gp2 genotype as vector of alleles
						genotypeKeyC[i][d] // d genotype as vector of alleles
						)
   				) continue;
					for(int gp1T = 0; gp1T < max3; gp1T++){ // for all possible TRUE genotypes
						for(int gp2T = 0; gp2T < max3; gp2T++){
							for(int dT = 0; dT < max3; dT++){
								// make sure this is a valid genotype for a trio - must have at least one allele inherited
								if(noAllelesInCommonGP(
									genotypeKeyC[i][gp1T], // gp1 genotype as vector of alleles
									genotypeKeyC[i][gp2T], // gp2 genotype as vector of alleles
									genotypeKeyC[i][dT] // d genotype as vector of alleles
									)
								) continue;

								// P(MI|o)P(o|t)P(t|true grandparents)
								// P(MI|o) = 1 b/c checked if observed MI above
								pMI[i] += genotypeErrorRatesC[i][gp1T][gp1] * genotypeErrorRatesC[i][gp2T][gp2] * 
									genotypeErrorRatesC[i][dT][d] * lGenos_ssGP[i][gp1T][gp2T][dT];
							}
						}
					}
				}
			}
		}
	}
	
}

// calculate the probability of sum(MI) = x | true grandparents
void calcProbSumMI(const vector <double>& pMI, vector <double>& pTotalMI){
	pTotalMI[0] = 1; // start at 0
	for(int i = 0, max = pMI.size(); i < max; i++){ // for all loci
		vector <double> tempVec (i+2, 0);
		for(int j = 0; j < (i+1); j++){ // for all possible states to be in
			// for all possible states to move too
			// stay in state j
			tempVec[j] += pTotalMI[j] * (1 - pMI[i]);
			// move to state j + 1
			tempVec[j + 1] += pTotalMI[j] * pMI[i];
		}
		for(int j = 0; j < (i+2); j++) pTotalMI[j] = tempVec[j];
	}
}
