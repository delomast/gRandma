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

// calculate the probability of sum(MI) = x
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

//////////////////////
//////////////////////
// begin single parent MI functions
//////////////////////
//////////////////////

// at each locus, calculate the probability of a MI | true parent
void calcProbMIperLocus_sP(const vector <vector <vector <int> > >& genotypeKeyC,
                        const vector< vector < vector <double> > >& genotypeErrorRatesC,
                        const vector <vector <vector <double> > >& lGenos_sP,
                        vector <double>& pMI){
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){
		// summing accross all possible obs and true, calculate
		// P(obs, true, MI | true parent) = P(MI|o)P(o|t)P(t|true parent)
		for(int p1 = 0, max3 = genotypeKeyC[i].size(); p1 < max3; p1++){ // for all possible OBS genotypes
			for(int d = 0; d < max3; d++){ 
				//check if observed MI
				if(!noAllelesInCommonSP(
					genotypeKeyC[i][p1], // p1 genotype as vector of alleles
					genotypeKeyC[i][d] // d genotype as vector of alleles
					)
				) continue;
				for(int p1T = 0; p1T < max3; p1T++){ // for all possible TRUE genotypes
					for(int dT = 0; dT < max3; dT++){
						// make sure this is a valid genotype for a trio - must have at least one allele inherited
						if(noAllelesInCommonSP(
							genotypeKeyC[i][p1T], // p1 genotype as vector of alleles
							genotypeKeyC[i][dT] // d genotype as vector of alleles
							)
						) continue;
						// P(MI|o)P(o|t)P(t|true parent)
						// P(MI|o) = 1 b/c checked if observed MI above
						pMI[i] += genotypeErrorRatesC[i][p1T][p1] * 
							genotypeErrorRatesC[i][dT][d] * lGenos_sP[i][p1T][dT];
					}
				}
			}
		}
	}
}


// at each locus, calculate the probability of an observed MI
//   given log-likelihoods for TRUE genotypes
//   for two unrelated individuals
void calcProbMIperLocus_Unrel(const vector <vector <vector <int> > >& genotypeKeyC,
                        const vector< vector < vector <double> > >& genotypeErrorRatesC,
                        const vector <vector <double> >& lGenos_one,
                        const vector <vector <double> >& lGenos_two,
                        vector <double>& pMI){
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){
		// summing accross all possible obs and true, calculate
		// P(obs, true, MI | true parent) = P(MI|o)P(o|t)P(t|true parent)
		for(int p1 = 0, max3 = genotypeKeyC[i].size(); p1 < max3; p1++){ // for all possible OBS genotypes
			for(int d = 0; d < max3; d++){ 
				//check if observed MI
				if(!noAllelesInCommonSP(
					genotypeKeyC[i][p1], // p1 genotype as vector of alleles
					genotypeKeyC[i][d] // d genotype as vector of alleles
					)
				) continue;
				for(int p1T = 0; p1T < max3; p1T++){ // for all possible TRUE genotypes
					for(int dT = 0; dT < max3; dT++){
						// P(MI|o)P(o|t)P(t|genotype frequencies)
						// P(MI|o) = 1 b/c checked if observed MI above
						pMI[i] += genotypeErrorRatesC[i][p1T][p1] * 
							genotypeErrorRatesC[i][dT][d] * exp(lGenos_one[i][p1T] + lGenos_two[i][dT]);
					}
				}
			}
		}
	}
}


void strat_lGenosBuilder_sP(const vector <vector <vector <int> > >& genotypeKeyC,
                        const vector< vector < vector <double> > >& genotypeErrorRatesC,
                        const vector <vector <double> >& lGenos_base,
                        const vector <vector <double> >& lGenos_randomDescendant,
                        vector <vector <vector <double> > >& lGenos_OBS_sample_MI,
								vector <vector <vector <double> > >& lGenos_OBS_sample_no_MI){
	// initialize with 0's
	int nLoci = genotypeKeyC.size();
	for(int i = 0; i < nLoci; i++){ //for each locus
		vector <vector <double> > tempLocus;
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			vector <double> tempP1 (max2,0.0);
			tempLocus.push_back(tempP1);
		}
		lGenos_OBS_sample_MI.push_back(tempLocus);
		lGenos_OBS_sample_no_MI.push_back(tempLocus);
	}
	
	for(int i = 0; i < nLoci; i++){
		double MI_sum = 0;
		double no_MI_sum = 0;
		for(int p1 = 0, nGeno = genotypeKeyC[i].size(); p1 < nGeno; p1++){
			for(int d = 0; d < nGeno; d++){
				double tempOBSp = 0;
				for(int p1T = 0; p1T < nGeno; p1T++){ // for all possible TRUE genotypes
					for(int dT = 0; dT < nGeno; dT++){
						tempOBSp += exp(lGenos_base[i][p1T] + lGenos_randomDescendant[i][dT]) *
							genotypeErrorRatesC[i][p1T][p1] * genotypeErrorRatesC[i][dT][d];
					}
				}
				// determine if MI
				if(noAllelesInCommonSP(
				genotypeKeyC[i][p1], // p1 genotype as vector of alleles
				genotypeKeyC[i][d] // d genotype as vector of alleles
				)){
					lGenos_OBS_sample_MI[i][p1][d] = tempOBSp;
					MI_sum += tempOBSp;
				} else {
					lGenos_OBS_sample_no_MI[i][p1][d] = tempOBSp;
					no_MI_sum += tempOBSp;
				}
			}
		}
		
		// normalize within each
		for(int p1 = 0, nGeno = genotypeKeyC[i].size(); p1 < nGeno; p1++){
			for(int d = 0; d < nGeno; d++){
				lGenos_OBS_sample_MI[i][p1][d] /= MI_sum;
				lGenos_OBS_sample_no_MI[i][p1][d] /= no_MI_sum;
			}
		}
	}
}


// at each locus, calculate the probability of an observed MI
//   given likelihoods (NOT log) for TRUE genotypes
//   for two individuals (relatinoship unspecified, joint genotype probabilties given in input)
void calcProbMIperLocus_Pair(const vector <vector <vector <int> > >& genotypeKeyC,
                        const vector< vector < vector <double> > >& genotypeErrorRatesC,
                        const vector <vector <vector <double> > >& lGenos_tRel,
                        vector <double>& pMI){
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){
		// summing accross all possible obs and true, calculate
		// P(obs, true, MI | true parent) = P(MI|o)P(o|t)P(t|true parent)
		for(int p1 = 0, max3 = genotypeKeyC[i].size(); p1 < max3; p1++){ // for all possible OBS genotypes
			for(int d = 0; d < max3; d++){ 
				//check if observed MI
				if(!noAllelesInCommonSP(
					genotypeKeyC[i][p1], // p1 genotype as vector of alleles
					genotypeKeyC[i][d] // d genotype as vector of alleles
					)
				) continue;
				for(int p1T = 0; p1T < max3; p1T++){ // for all possible TRUE genotypes
					for(int dT = 0; dT < max3; dT++){
						// P(MI|o)P(o|t)P(t|genotype frequencies)
						// P(MI|o) = 1 b/c checked if observed MI above
						pMI[i] += genotypeErrorRatesC[i][p1T][p1] * 
							genotypeErrorRatesC[i][dT][d] * lGenos_tRel[i][p1T][dT];
					}
				}
			}
		}
	}
}


void strat_lGenosBuilder_Pairwise(const vector <vector <vector <int> > >& genotypeKeyC,
                        const vector< vector < vector <double> > >& genotypeErrorRatesC,
                        const vector <vector <vector <double> > >& lGenos_tRel_OBS,
                        vector <vector <vector <double> > >& lGenos_OBS_sample_MI,
								vector <vector <vector <double> > >& lGenos_OBS_sample_no_MI){
	// initialize with 0's
	int nLoci = genotypeKeyC.size();
	for(int i = 0; i < nLoci; i++){ //for each locus
		vector <vector <double> > tempLocus;
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			vector <double> tempP1 (max2,0.0);
			tempLocus.push_back(tempP1);
		}
		lGenos_OBS_sample_MI.push_back(tempLocus);
		lGenos_OBS_sample_no_MI.push_back(tempLocus);
	}
	
	for(int i = 0; i < nLoci; i++){
		double MI_sum = 0;
		double no_MI_sum = 0;
		for(int p1 = 0, nGeno = genotypeKeyC[i].size(); p1 < nGeno; p1++){
			for(int d = 0; d < nGeno; d++){
				// determine if MI
				if(noAllelesInCommonSP(
				genotypeKeyC[i][p1], // p1 genotype as vector of alleles
				genotypeKeyC[i][d] // d genotype as vector of alleles
				)){
					lGenos_OBS_sample_MI[i][p1][d] = exp(lGenos_tRel_OBS[i][p1][d]);
					MI_sum += exp(lGenos_tRel_OBS[i][p1][d]);
				} else {
					lGenos_OBS_sample_no_MI[i][p1][d] = exp(lGenos_tRel_OBS[i][p1][d]);
					no_MI_sum += exp(lGenos_tRel_OBS[i][p1][d]);
				}
			}
		}
		
		// normalize within each
		for(int p1 = 0, nGeno = genotypeKeyC[i].size(); p1 < nGeno; p1++){
			for(int d = 0; d < nGeno; d++){
				lGenos_OBS_sample_MI[i][p1][d] /= MI_sum;
				lGenos_OBS_sample_no_MI[i][p1][d] /= no_MI_sum;
			}
		}
	}
}