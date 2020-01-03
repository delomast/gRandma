#include <Rcpp.h>
#include <vector>
#include <math.h>
#include "misc_math.h"
#include "utils.h"

using namespace std;

//' calculating llr of single-sided grandparent pair vs unrelated
//' returning a matrix of descendent, grandparentPopulation, grandparent1, grandparent2, llr
//' @param baseline baseline individuals, col1 is pop, col2 is id, cols 3... are genotypes
//' @param mixture mixture individulas col1 is id, cols 2... are genotypes
//' @param crossRecords col1 is pop, col2 is id 1, col3 is id2
//' @param baselineParams Dirichlet parameters for allele frequencies
//' @param unsampledPopParams Dirichlet parameters for allele frequencies
//' @param genotypeKey list of matrix for each locus, col1 is genotype, col2 is allele 1, col3 is allele 2
//' @param genotypeErrorRates list of matrix for each locus, rows are actual genotype, columns are observed,
//'   values are probability
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame ssGP(Rcpp::NumericMatrix baseline, Rcpp::NumericMatrix mixture, 
                           Rcpp::NumericMatrix crossRecords, Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::List genotypeKey,
                           Rcpp::List genotypeErrorRates
                           ){

	// determine if cross records used or not
	bool useCR = true;
	if(crossRecords.nrow() == 0) useCR = false;
	
	//////////////////////
	// first, turn all inputs into non-Rcpp c++ classes
	// here, making vectors with rows each being a vector of ints
	vector <vector <int> > baselineC;
	baselineC = rcppMatrixToVectorInt(baseline);
	
	vector <vector <int> > mixtureC;
	mixtureC = rcppMatrixToVectorInt(mixture);
	
	vector <vector <int> > crossRecordsC;
	if(useCR) crossRecordsC = rcppMatrixToVectorInt(crossRecords);
	
	vector <vector <vector <double> > > baselineParamsC; //first index is pop, then locus, then allele
	for(int i = 0, max = baselineParams.length(); i < max; i++){
		baselineParamsC.push_back(rcppListToVectorDouble(baselineParams[i]));
	}

	vector <vector <vector <double> > > unsampledPopParamsC; //first index is pop, then locus, then allele
	for(int i = 0, max = unsampledPopParams.length(); i < max; i++){
		unsampledPopParamsC.push_back(rcppListToVectorDouble(unsampledPopParams[i]));
	}

	vector <vector <vector <int> > > genotypeKeyC; //first index is locus, then genotype, values are alleles
	genotypeKeyC = convertGenotypeKey(genotypeKey);
	
	vector< vector < vector <double> > > genotypeErrorRatesC; // first index is locus, then actual genotype, then observed
	for(int i = 0, max = genotypeErrorRates.length(); i < max; i++){
		genotypeErrorRatesC.push_back(rcppMatrixToVectorDouble(genotypeErrorRates[i]));
	}
	
	// end conversion of types
	////////////////////
	

	// calculate number of MI allowed to consider a pair of potential grandparents
	double miLimit = (baselineC[0].size() - 2) * .05; // this is just a placeholder for initial testing

	// calculate genotype likelihoods from parent - parent - offspring relationship
	// the same for all baseline pops b/c based on mendelian inheritance
	/*	For lGenos_ppo:
	 	Index 1 is locus
	 	Index 2 is parent1 genotype
		Index 3 is parent 2 genotype
		Index 4 is offspring genotype
		Value is likelihood
	 */
	vector <vector <vector <vector <double> > > > lGenos_ppo;
	lGenos_ppo = createPPOvector(genotypeKeyC);
	
	// initiate results storage
	vector <int> mixtureInd;
	vector <int> GrandparentPop;
	vector <int> grandparent1;
	vector <int> grandparent2;
	vector <double> llr;
	vector <int> mendIncompat;
		
	// for each baseline pop
	for(int pop = 0, maxP = baselineParamsC.size(); pop < maxP; pop++){

		// either get cross records, or make a list of all possible pairs
		vector <vector <int> > pairs; //all pairs to consider
		if(useCR){
			vector <int> tempVec (2,-9);
			for(int i = 0, max = crossRecordsC.size(); i < max; i++){
				if(crossRecordsC[i][0] == pop){
					tempVec[0] = crossRecordsC[i][1];
					tempVec[1] = crossRecordsC[i][2];
					pairs.push_back(tempVec);
				}
			}
		}
		if(pairs.size() == 0){ // either !useCR or no cross records for the current population
			vector <int> baseInds;
			for(int i=0, max = baselineC.size(); i < max; i++) {
				if(baselineC[i][0] == pop) baseInds.push_back(baselineC[i][1]);
			}
			pairs = listAllPairs(baseInds);
		}

		// calculate genotype likelihoods from random sampled genotype
		vector <vector <double> > lGenos_base; // one vector for each locus, second index is genotype
		vector <vector <double> > lGenos_unsamp;

		for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ // for each locus
			// baselineParamsC[pop][i] are the alpha values for the baseline pop
			// unsampledPopParamsC[pop][i] are the alpha values for the corresponding unsampled pop
			// genotypeKeyC[i] is the key for the locus i
			vector <double> tempBase;
			vector <double> tempUnsamp;
			for(int j = 0, max2 = genotypeKeyC[i].size(); j < max2; j++){ // for each genotype
				vector <double> k (baselineParamsC[pop][i].size(), 0);
				for(int l = 0, max3 = baselineParamsC[pop][i].size(); l < max3; l++){
					if(genotypeKeyC[i][j][0] == l) k[l]++;
					if(genotypeKeyC[i][j][1] == l) k[l]++;
				}
				tempBase.push_back(logDirichMultPMF(k, baselineParamsC[pop][i]));
				tempUnsamp.push_back(logDirichMultPMF(k, unsampledPopParamsC[pop][i]));
			}
			lGenos_base.push_back(tempBase);
			lGenos_unsamp.push_back(tempUnsamp);
		}

		// for each individual in the mixture
		for(int m = 0, maxM = mixtureC.size(); m < maxM; m++){
			cout<<m<<endl;
			for(int i = 0, max = pairs.size(); i < max; i++){ // for each pair
				if(i % 100 == 0) Rcpp::checkUserInterrupt();
				/* filter pairs based on MI
				 * for a pair of grandparents, number of MI is the number of loci for which the grandparent pair
				 * and potential descendant share no alleles in common
				 */
				int countMI = 0;
				bool skip = false;
				for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ // for each locus
					if(noAllelesInCommonGP(
							genotypeKeyC[j][baselineC[pairs[i][0]][j+2]], // gp1 genotype as vector of alleles
							genotypeKeyC[j][baselineC[pairs[i][1]][j+2]], // gp2 genotype as vector of alleles
							genotypeKeyC[j][mixtureC[m][j+1]] // o genotype as vector of alleles
						)
					) countMI++;
					if (countMI > miLimit){
						skip = true;
						break;
					}
				}
				if(skip) continue;
				// passed the MI filter,so calculate the llr
				double uLLH = 0.0; // log-likelihood unrelated
				double gpLLH = 0.0; // log-likelihood grandparent pair
				for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ //for each locus
					//observed genotypes
					int obs_gp1 = baselineC[pairs[i][0]][j+2];
					int obs_gp2 = baselineC[pairs[i][1]][j+2];
					int obs_d = mixtureC[m][j+1];
					
					// likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
					double u_likelihood = 0.0;
					double gp_likelihood = 0.0;
					
					// prob of observed given actual
					double pOA_gp1, pOA_gp2, pOA_d;
					
					// unrelated
					for(int gp1 = 0, max3 = genotypeKeyC[j].size(); gp1 < max3; gp1++){ // for each possible gp1 genotype
						for(int gp2 = 0; gp2 < max3; gp2++){ // for each possible gp2 genotype
							for(int d = 0; d < max3; d++){ // for each possible d genotype
								if(obs_gp1 == -9){
									pOA_gp1 = 1.0;
								} else {
									pOA_gp1 = genotypeErrorRatesC[j][gp1][obs_gp1];
								}
								if(obs_gp2 == -9){
									pOA_gp2 = 1.0;
								} else {
									pOA_gp2 = genotypeErrorRatesC[j][gp2][obs_gp2];
								}
								if(obs_d == -9){
									pOA_d = 1.0;
								} else {
									pOA_d = genotypeErrorRatesC[j][d][obs_d];
								}
								// P(gp1|baseline pop)P(gp2|baseline pop)P(d|unsampled pop)P(gp1|obs_gp1)P(gp2|obs_gp2)P(d|obs_d)
								u_likelihood += exp(lGenos_base[j][gp1] + lGenos_base[j][gp2] + lGenos_unsamp[j][d]) * 
									 pOA_gp1 * pOA_gp2 * pOA_d;
							}
						}
					}
					uLLH += log(u_likelihood);
					// grandparent pair
					for(int gp1 = 0, max3 = genotypeKeyC[j].size(); gp1 < max3; gp1++){ // for each possible gp1 genotype
						for(int gp2 = 0; gp2 < max3; gp2++){ // for each possible gp2 genotype
							for(int d = 0; d < max3; d++){ // for each possible d genotype
								if(obs_gp1 == -9){
									pOA_gp1 = 1.0;
								} else {
									pOA_gp1 = genotypeErrorRatesC[j][gp1][obs_gp1];
								}
								if(obs_gp2 == -9){
									pOA_gp2 = 1.0;
								} else {
									pOA_gp2 = genotypeErrorRatesC[j][gp2][obs_gp2];
								}
								if(obs_d == -9){
									pOA_d = 1.0;
								} else {
									pOA_d = genotypeErrorRatesC[j][d][obs_d];
								}
								double probSum = 0.0;
								for(int uf = 0; uf < max3; uf++){ //for each geno of unobserved wild parent
									for(int um = 0; um < max3; um++){ //for each geno of unobserved hatchery parent
										//summing over unsampled genotypes P(uf)P(um|m,f)P(d|uf,um)
										probSum += exp(lGenos_unsamp[j][uf]) * lGenos_ppo[j][gp1][gp2][um] * 
											lGenos_ppo[j][uf][um][d];
									}
								}
								
								// P(gp1|baseline pop)P(gp2|baseline pop)(probsum)P(gp1|obs_gp1)P(gp2|obs_gp2)P(d|obs_d)
								gp_likelihood += exp(lGenos_base[j][gp1] + lGenos_base[j][gp2]) * probSum *
									 pOA_gp1 * pOA_gp2 * pOA_d;
							}
						}
					}
					gpLLH += log(gp_likelihood);
				}
				
				// save results
				mixtureInd.push_back(m);
				GrandparentPop.push_back(pop);
				grandparent1.push_back(pairs[i][0]);
				grandparent2.push_back(pairs[i][1]);
				llr.push_back(gpLLH - uLLH);
				mendIncompat.push_back(countMI);

			} // end for each pair
			
			

			
		}

		
		
		
	}
	

	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Kid") = mixtureInd,
                                                   Rcpp::Named("Gpop") = GrandparentPop,
                                                   Rcpp::Named("Gma") = grandparent1,
                                                   Rcpp::Named("Gpa") = grandparent2,
                                                   Rcpp::Named("llr") = llr,
                                                   Rcpp::Named("MI") = mendIncompat
																	);
	
	return results;
}
