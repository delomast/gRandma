// functions to calculate pairwise likelihood of genotypes given relationship, allele freqs, and genotyping
//   error rates
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <math.h>
#include <vector>
#include "misc_math.h"
#include "pairwiseK.h"

using namespace std;

//
//
// NOTE to self: this function is probably not correct. Need to go through it carefully.
//
//
// given TRUE genotypes
void pairwiseK(const vector <vector <vector <int> > >& genotypeKeyC, 
	                      const  vector <vector <double> >& lGenos_base, 
	                      const vector <vector <vector <double> > >& baselineParamsC,
	                      const vector <vector <vector <double> > >& unsampledPopParamsC,
	                      const int pop,
	                      const vector <double>& k_prob,
	                      vector <vector <vector <double> > >& lGenos_k){

	// initialize with 0's
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ //for each locus
		vector <vector <double> > tempLocus;
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			vector <double> tempP1 (max2,0.0);
			tempLocus.push_back(tempP1);
		}
		lGenos_k.push_back(tempLocus);
	}
	
	// now calculate likelihoods: P(G1) * P(G2 | G1, k)
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ // for each locus
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ // for each p1 genotype
			for(int d = 0; d < max2; d++){ // for each d genotype
				
				double tempLH = exp(lGenos_base[i][p1]);
				vector <double> k (unsampledPopParamsC[pop][i].size(), 0);
				k[genotypeKeyC[i][d][0]] = 1;
				
				// determine if 0,1,2 alleles are identical in STATE, then calculate appropriately
				if((genotypeKeyC[i][d][0] == genotypeKeyC[i][p1][0] && 
     				genotypeKeyC[i][d][1] == genotypeKeyC[i][p1][1]) ||
					(genotypeKeyC[i][d][0] == genotypeKeyC[i][p1][1] && 
     				genotypeKeyC[i][d][1] == genotypeKeyC[i][p1][0])){
					// ibs = 2;
					double pIBD_1, pIBD_0;
					if(genotypeKeyC[i][d][0] == genotypeKeyC[i][d][1]){ // d is homozygous
						pIBD_1 = exp(logMultPMF(k, unsampledPopParamsC[pop][i]));
						pIBD_0 = exp(logMultPMF(k, unsampledPopParamsC[pop][i]) + 
							logMultPMF(k, baselineParamsC[pop][i]));
					} else {
						vector <double> k2 (unsampledPopParamsC[pop][i].size(), 0);
						k2[genotypeKeyC[i][d][1]] = 1;
						pIBD_1 = 0.5 * (exp(logMultPMF(k2, unsampledPopParamsC[pop][i])) +
							exp(logMultPMF(k, unsampledPopParamsC[pop][i])));
						pIBD_0 = exp(logMultPMF(k, unsampledPopParamsC[pop][i]) + 
								logMultPMF(k2, baselineParamsC[pop][i])) +
							exp(logMultPMF(k2, unsampledPopParamsC[pop][i]) + 
								logMultPMF(k, baselineParamsC[pop][i]));
					}
					tempLH *= k_prob[2] + (k_prob[1] * pIBD_1) + (k_prob[0] * pIBD_0);
					
				} else if (genotypeKeyC[i][d][0] == genotypeKeyC[i][p1][0] ||
					genotypeKeyC[i][d][1] == genotypeKeyC[i][p1][1] ||
					genotypeKeyC[i][d][0] == genotypeKeyC[i][p1][1] ||
					genotypeKeyC[i][d][1] == genotypeKeyC[i][p1][0]){
					// ibs = 1;
					if(genotypeKeyC[i][d][0] == genotypeKeyC[i][d][1]){ // d is homozygous
						tempLH *= k_prob[1] * 0.5 * exp(logMultPMF(k, unsampledPopParamsC[pop][i])) +
							k_prob[0] * exp(logMultPMF(k, unsampledPopParamsC[pop][i]) + 
							logMultPMF(k, baselineParamsC[pop][i]));
							
					} else {
						vector <double> k2 (unsampledPopParamsC[pop][i].size(), 0);
						k2[genotypeKeyC[i][d][1]] = 1;
						double t;
						// first determine which allele would by ibd (only one is ibs, so one one can be idb)
						if(genotypeKeyC[i][d][0] == genotypeKeyC[i][p1][0] || genotypeKeyC[i][d][0] == genotypeKeyC[i][p1][1]){
							t = exp(logMultPMF(k2, unsampledPopParamsC[pop][i]));
						} else {
							t = exp(logMultPMF(k, unsampledPopParamsC[pop][i]));
						}
						tempLH *= (k_prob[1] * t) + (k_prob[0] * (exp(logMultPMF(k, unsampledPopParamsC[pop][i]) + 
								logMultPMF(k2, baselineParamsC[pop][i])) +
							exp(logMultPMF(k2, unsampledPopParamsC[pop][i]) + 
								logMultPMF(k, baselineParamsC[pop][i]))));
					}
					
				} else {
					// ibs = 0;
					if(genotypeKeyC[i][d][0] == genotypeKeyC[i][d][1]){ // d is homozygous
						tempLH *= k_prob[0] * exp(logMultPMF(k, unsampledPopParamsC[pop][i]) + 
							logMultPMF(k, baselineParamsC[pop][i]));
					} else {
						vector <double> k2 (unsampledPopParamsC[pop][i].size(), 0);
						k2[genotypeKeyC[i][d][1]] = 1;
						
						tempLH *= k_prob[0] * (exp(logMultPMF(k, unsampledPopParamsC[pop][i]) + 
								logMultPMF(k2, baselineParamsC[pop][i])) +
							exp(logMultPMF(k2, unsampledPopParamsC[pop][i]) + 
								logMultPMF(k, baselineParamsC[pop][i])));
					}

				}
				lGenos_k[i][p1][d] = tempLH;
			} // end for each d
		} // end for p1
	} // end for locus
}



// calculate genotype LOG - likelihoods for pairwise relationship
//	given OBSERVED genotypes - but only if NO MISSING genotypes
void create_pairswise_OBSvector(const vector <vector <vector <double> > >& lGenos_k, 
                          const vector< vector < vector <double> > >& genotypeErrorRatesC, 
                          vector <vector <vector <double> > >& lGenos_k_OBS
                          ){
	// initialize with 0's - note that 0's will cause problems if not replaced b/c these are log-likelihoods...
	for(int i = 0, max = genotypeErrorRatesC.size(); i < max; i++){ //for each locus
		vector <vector <double> > tempLocus;
		for(int p1 = 0, max2 = genotypeErrorRatesC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			vector <double> tempP1 (max2,0.0);
			tempLocus.push_back(tempP1);
		}
		lGenos_k_OBS.push_back(tempLocus);
	}

	// calculate for all loci and observed genotype combinations
	for(int i = 0, max = genotypeErrorRatesC.size(); i < max; i++){
		for(int obs_p1 = 0, max3 = genotypeErrorRatesC[i].size(); obs_p1 < max3; obs_p1++){
			for(int obs_d = 0; obs_d < max3; obs_d++){
				double p_likelihood = 0.0;
				// marginalize over all possible true genotype combinations
				for(int p1 = 0; p1 < max3; p1++){
					for(int d = 0; d < max3; d++){
						p_likelihood += lGenos_k[i][p1][d] *
							genotypeErrorRatesC[i][p1][obs_p1] * 
							genotypeErrorRatesC[i][d][obs_d];

					}
				}
				lGenos_k_OBS[i][obs_p1][obs_d] = log(p_likelihood);
			}
		}
	}
}
