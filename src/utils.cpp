// miscellaneous helper functions

#include <Rcpp.h>
#include <vector>
#include <math.h>
#include "misc_math.h"
#include "utils.h"

using namespace std;

// list all possible pairs of individuals from a vector of individual ID's (as integers)
void listAllPairs(const vector <int>& allInds, vector <vector <int> >& pairs){
	vector <int> tempVec (2, -9);
	for(int i = 0, max = allInds.size() - 1; i < max; i++){
		for(int j = i+1, max2 = allInds.size(); j < max2; j++){
			tempVec[0] = allInds[i];
			tempVec[1] = allInds[j];
			pairs.push_back(tempVec);
		}
	}
}

// convert a NumericMatrix of ints to to nested vectors, index 1 is row index 2 is column
vector <vector <int> > rcppMatrixToVectorInt(Rcpp::NumericMatrix inputMatrix){
	vector <vector <int> > output;
	for(int i = 0, max = inputMatrix.nrow(); i < max; i++){
		vector <int> tempVec;
		for(int j = 0, max2 = inputMatrix.ncol(); j < max2; j++){
			tempVec.push_back(inputMatrix(i,j));
		}
		output.push_back(tempVec);
	}
	return output;
}

// convert a List of NumericVectors to nested vectors index 1 is list position, index 2 is vector position
vector <vector <double> > rcppListToVectorDouble(Rcpp::List inputList){
	vector <vector <double> > output;
	for(int i = 0, max = inputList.length(); i < max; i++){
		Rcpp::NumericVector current;
		current = inputList[i];
		vector <double> tempVec;
		for(int j = 0, max2 = current.length(); j < max2; j++){
			tempVec.push_back(current[j]);
		}
		output.push_back(tempVec);
	}
	return output;
}

// convert the GenotypeKey list to nested vectors
vector <vector <vector <int> > > convertGenotypeKey(Rcpp::List inputList){
	vector <vector <vector <int> > > output;
	for(int i = 0, max = inputList.length(); i < max; i++){ // for each locus
		Rcpp::NumericMatrix current = inputList[i];
		vector <vector <int> > tempLocus;
		for(int j = 0, max2 = current.nrow(); j < max2; j++){ // for each genotype
			vector <int> tempGeno (2, -9);
			tempGeno[0] = current(j,1);
			tempGeno[1] = current(j,2);
			tempLocus.push_back(tempGeno);
		}
		output.push_back(tempLocus);
	}
	return output;
}

// convert a NumericMatrix to a nexted vector with first index being the row and 
//   second index being the column
vector < vector <double> > rcppMatrixToVectorDouble(Rcpp::NumericMatrix inputMatrix){
	vector < vector <double> > output;
	for(int i = 0, max = inputMatrix.nrow(); i < max; i++){
		vector <double> tempVec;
		for(int j = 0, max2 = inputMatrix.ncol(); j < max2; j++){
			tempVec.push_back(inputMatrix(i,j));
		}
		output.push_back(tempVec);
	}
	return output;
}

// takes genotypeKeyC, calculates Parent parent offspring likelihoods for all genotype combinations
//   and then retuns it as a nested vector with :	 	Index 1 is locus, Index 2 is parent1 genotype
//   Index 3 is parent 2 genotype, Index 4 is offspring genotype, Value is likelihood
void createPPOvector(const vector <vector <vector <int> > >& genotypeKeyC,
                vector <vector <vector <vector <double> > > >& lGenos_ppo){
	
	// initialize with 0's
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ //for each locus
		vector <vector <vector <double> > > tempLocus;
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			vector <vector <double> > tempP1;
			for(int p2 = 0; p2 < max2; p2++){ //for each p2 genotype
				vector <double> tempP2 (max2,0.0);
				tempP1.push_back(tempP2);
			}
			tempLocus.push_back(tempP1);
		}
		lGenos_ppo.push_back(tempLocus);
	}

	// now calculate likelihoods
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ //for each locus
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			for(int p2 = p1; p2 < max2; p2++){ //for each p2 genotype
				for(int o = 0; o < max2; o++){ //for each offspring genotype
					lGenos_ppo[i][p1][p2][o] = ppoMendelian(genotypeKeyC[i][p1], genotypeKeyC[i][p2], genotypeKeyC[i][o]);
					if(p1 != p2) lGenos_ppo[i][p2][p1][o] = lGenos_ppo[i][p1][p2][o];
				}
			}
		}
	}
}

// dertermine whether there are any alleles in common between the four alleles in
//   two potential grandparents (gp1 and gp2) and a potential descendant (d). This
//   is only for one locus. Returns true if no allles are in common (ie a MI),
//   and false if there are one or more alleles in common.
// @param gp1 genotype as vector of alleles
// @param gp2 genotype as vector of alleles
// @param d genotype as vector of alleles

bool noAllelesInCommonGP(const vector <int>& gp1, const vector <int>& gp2, const vector <int>& d){
	if(d[0] == gp1[0] || d[0] == gp1[1] || d[1] == gp1[0] || d[1] == gp1[1] ||
    d[0] == gp2[0] || d[0] == gp2[1] || d[1] == gp2[0] || d[1] == gp2[1]){
		return false;
	} else {
		return true;
	}
}

// calculate genotype likelihoods from grandparent(maternal) - grandparent(maternal) - grandchild relationship
// given TRUE genotypes
void createSSGPvector(const vector <vector <vector <int> > >& genotypeKeyC, 
                      const  vector <vector <double> >& lGenos_base, 
                      const vector <vector <vector <double> > >& unsampledPopParamsC,
                      const int pop,
                      vector <vector <vector <vector <double> > > >& lGenos_ssGP){
		
	// initialize with 0's
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ //for each locus
		vector <vector <vector <double> > > tempLocus;
		for(int gp1 = 0, max2 = genotypeKeyC[i].size(); gp1 < max2; gp1++){ //for each gp1 genotype
			vector <vector <double> > tempP1;
			for(int gp2 = 0; gp2 < max2; gp2++){ //for each gp2 genotype
				vector <double> tempP2 (max2,0.0);
				tempP1.push_back(tempP2);
			}
			tempLocus.push_back(tempP1);
		}
		lGenos_ssGP.push_back(tempLocus);
	}
	
	// now calculate likelihoods: P(gp1|baselineParams)P(gp2|baselineParams)P(y|gp1, gp2, unsampledPopParams)
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ //for each locus
		for(int gp1 = 0, max2 = genotypeKeyC[i].size(); gp1 < max2; gp1++){ //for each gp1 genotype
			for(int gp2 = gp1; gp2 < max2; gp2++){ //for each gp2 genotype
				for(int d = 0; d < max2; d++){ //for each offspring genotype
					double tempLH = exp(lGenos_base[i][gp1] + lGenos_base[i][gp2]);
					double pInherit0, pInherit1;
					// calculate prob of inheriting each allele from gps
					pInherit0 = ((genotypeKeyC[i][d][0] == genotypeKeyC[i][gp1][0]) + 
						(genotypeKeyC[i][d][0] == genotypeKeyC[i][gp1][1]) + 
						(genotypeKeyC[i][d][0] == genotypeKeyC[i][gp2][0]) + 
						(genotypeKeyC[i][d][0] == genotypeKeyC[i][gp2][1])) / 4.0;
					pInherit1 = ((genotypeKeyC[i][d][1] == genotypeKeyC[i][gp1][0]) + 
						(genotypeKeyC[i][d][1] == genotypeKeyC[i][gp1][1]) + 
						(genotypeKeyC[i][d][1] == genotypeKeyC[i][gp2][0]) + 
						(genotypeKeyC[i][d][1] == genotypeKeyC[i][gp2][1])) / 4.0;
					if(pInherit0 == 0 && pInherit1 == 0) continue; // if MI, then just leave as 0
					
					vector <double> k (unsampledPopParamsC[pop][i].size(), 0);
					k[genotypeKeyC[i][d][0]] = 1;
					if(genotypeKeyC[i][d][0] == genotypeKeyC[i][d][1]){
						// P(gp1) * P(gp2) * prob allele from unsampled pop * prob allele from one of GPs
						tempLH *= exp(logDirichMultPMF(k, unsampledPopParamsC[pop][i])) * pInherit0;
					} else {
						vector <double> k2 (unsampledPopParamsC[pop][i].size(), 0);
						k2[genotypeKeyC[i][d][1]] = 1;
						// P(gp1) * P(gp2) * (prob allele1 from unsampled pop * prob allele2 from GPs +
						// prob allele2 from unsampled pop * prob allele1 from GPs)
						tempLH *= (exp(logDirichMultPMF(k, unsampledPopParamsC[pop][i])) * pInherit1) +
							(exp(logDirichMultPMF(k2, unsampledPopParamsC[pop][i])) * pInherit0);
					}
					
					lGenos_ssGP[i][gp1][gp2][d] = tempLH;
					if(gp1 != gp2) lGenos_ssGP[i][gp2][gp1][d] = tempLH;
				}
			}
		}
	}
}

// calculate genotype LOG - likelihoods for grandparent(maternal) - grandparent(maternal) - grandchild relationship
// and LOG - likelihoods for all three being unrelated
//	given OBSERVED genotypes - but only if NO MISSING genotypes in the trio
void createOBSvector(const vector <vector <vector <vector <double> > > >& lGenos_ssGP, 
                          const vector< vector < vector <double> > >& genotypeErrorRatesC, 
                          vector <vector <vector <vector <double> > > >& lGenos_ssGP_OBS,
                          const vector <vector <double> >& lGenos_base,
                          const vector <vector <double> >& lGenos_unsamp,
                          vector <vector <vector <vector <double> > > >& lGenos_Unrelated_OBS
                          ){
	// initialize with 0's - note that 0's will cause problems if not replaced b/c these are log-likelihoods...
	for(int i = 0, max = genotypeErrorRatesC.size(); i < max; i++){ //for each locus
		vector <vector <vector <double> > > tempLocus;
		for(int gp1 = 0, max2 = genotypeErrorRatesC[i].size(); gp1 < max2; gp1++){ //for each gp1 genotype
			vector <vector <double> > tempP1;
			for(int gp2 = 0; gp2 < max2; gp2++){ //for each gp2 genotype
				vector <double> tempP2 (max2,0.0);
				tempP1.push_back(tempP2);
			}
			tempLocus.push_back(tempP1);
		}
		lGenos_ssGP_OBS.push_back(tempLocus);
		lGenos_Unrelated_OBS.push_back(tempLocus);
	}
	
	// calculate for all loci and observed genotype combinations
	for(int i = 0, max = genotypeErrorRatesC.size(); i < max; i++){
		for(int obs_gp1 = 0, max3 = genotypeErrorRatesC[i].size(); obs_gp1 < max3; obs_gp1++){
			for(int obs_gp2 = obs_gp1; obs_gp2 < max3; obs_gp2++){
				for(int obs_d = 0; obs_d < max3; obs_d++){
					double gp_likelihood = 0.0;
					double u_likelihood = 0.0;
					// marginalize over all possible true genotype combinations
					for(int gp1 = 0; gp1 < max3; gp1++){
						for(int gp2 = 0; gp2 < max3; gp2++){
							for(int d = 0; d < max3; d++){
								
								gp_likelihood += lGenos_ssGP[i][gp1][gp2][d] *
									genotypeErrorRatesC[i][gp1][obs_gp1] * 
									genotypeErrorRatesC[i][gp2][obs_gp2] * 
									genotypeErrorRatesC[i][d][obs_d];
								
								u_likelihood += exp(lGenos_base[i][gp1] + lGenos_base[i][gp2] + lGenos_unsamp[i][d]) * 
									genotypeErrorRatesC[i][gp1][obs_gp1] * 
									genotypeErrorRatesC[i][gp2][obs_gp2] * 
									genotypeErrorRatesC[i][d][obs_d];
							}
						}
					}
					lGenos_ssGP_OBS[i][obs_gp1][obs_gp2][obs_d] = log(gp_likelihood);
					lGenos_Unrelated_OBS[i][obs_gp1][obs_gp2][obs_d] = log(u_likelihood);
					if(obs_gp1 != obs_gp2){
						lGenos_ssGP_OBS[i][obs_gp2][obs_gp1][obs_d] = log(gp_likelihood);
						lGenos_Unrelated_OBS[i][obs_gp2][obs_gp1][obs_d] = log(u_likelihood);
					}
				}
			}
		}
	}
}


// calculate genotype LOG - likelihoods for three unrelated individuals:
//   two from one pop and one a descendant of another pop
//	given OBSERVED genotypes - but only if NO MISSING genotypes in the trio
// This is used in "otherPopERRORssGP"
void create_CORR_OBSvector(const vector< vector < vector <double> > >& genotypeErrorRatesC, 
                           const vector <vector <double> >& lGenos_randomDescendant,
                          const vector <vector <double> >& lGenos_base,
                          vector <vector <vector <vector <double> > > >& CORR_lGenos_OBS // output to this
                          ){
	// initialize with 0's - note that 0's will cause problems if not replaced b/c these are log-likelihoods...
	for(int i = 0, max = genotypeErrorRatesC.size(); i < max; i++){ //for each locus
		vector <vector <vector <double> > > tempLocus;
		for(int gp1 = 0, max2 = genotypeErrorRatesC[i].size(); gp1 < max2; gp1++){ //for each gp1 genotype
			vector <vector <double> > tempP1;
			for(int gp2 = 0; gp2 < max2; gp2++){ //for each gp2 genotype
				vector <double> tempP2 (max2,0.0);
				tempP1.push_back(tempP2);
			}
			tempLocus.push_back(tempP1);
		}
		CORR_lGenos_OBS.push_back(tempLocus);
	}
	
	// calculate for all loci and observed genotype combinations
	for(int i = 0, max = genotypeErrorRatesC.size(); i < max; i++){
		for(int obs_gp1 = 0, max3 = genotypeErrorRatesC[i].size(); obs_gp1 < max3; obs_gp1++){
			for(int obs_gp2 = obs_gp1; obs_gp2 < max3; obs_gp2++){
				for(int obs_d = 0; obs_d < max3; obs_d++){
					double u_likelihood = 0.0;
					// marginalize over all possible true genotype combinations
					for(int gp1 = 0; gp1 < max3; gp1++){
						for(int gp2 = 0; gp2 < max3; gp2++){
							for(int d = 0; d < max3; d++){
								
								u_likelihood += exp(lGenos_base[i][gp1] + lGenos_base[i][gp2] + 
										lGenos_randomDescendant[i][d]) * 
									genotypeErrorRatesC[i][gp1][obs_gp1] * 
									genotypeErrorRatesC[i][gp2][obs_gp2] * 
									genotypeErrorRatesC[i][d][obs_d];
							}
						}
					}
					CORR_lGenos_OBS[i][obs_gp1][obs_gp2][obs_d] = log(u_likelihood);
					if(obs_gp1 != obs_gp2){
						CORR_lGenos_OBS[i][obs_gp2][obs_gp1][obs_d] = log(u_likelihood);
					}
				}
			}
		}
	}
}
