#include <Rcpp.h>
#include <vector>
#include <math.h>
#include "misc_math.h"
#include "utils.h"
#include "mend_incompat.h"

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
//' @param saveLLR the minimum LLR to include the output in the results (only used if filterLLR is true)
//' @param MIexcludeProb the maximum probability of exclusion for a true grandparent pair due to 
//'   Mendelian incompatibilities
//' @param filterLLR true to filter results based on saveLLR, false to not
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame ssGP(Rcpp::NumericMatrix baseline, Rcpp::NumericMatrix mixture, 
                           Rcpp::NumericMatrix crossRecords, Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::List genotypeKey,
                           Rcpp::List genotypeErrorRates, double saveLLR, double MIexcludeProb,
                           bool filterLLR
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
	
	// initiate results storage
	vector <int> mixtureInd;
	vector <int> GrandparentPop;
	vector <int> grandparent1;
	vector <int> grandparent2;
	vector <double> llr;
	vector <int> mendIncompat;
	
	// for each baseline pop
	for(int pop = 0, maxP = baselineParamsC.size(); pop < maxP; pop++){
		Rcpp::Rcout<<"Now beginning baseline population "<<pop + 1<<"\n";
		
		// calculate genotype log-likelihoods from random sampled genotype
		// first index is locus, second is genotype, value is LOG-likelihood
		vector <vector <double> > lGenos_base;
		vector <vector <double> > lGenos_unsamp;

		for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ // for each locus
			// baselineParamsC[pop][i] are the alpha values for the baseline pop
			// unsampledPopParamsC[pop][i] are the alpha values for the corresponding unsampled pop
			// genotypeKeyC[i] is the key for the locus i
			vector <double> tempBase;
			vector <double> tempUnsamp;
			for(int j = 0, max2 = genotypeKeyC[i].size(); j < max2; j++){ // for each genotype
				vector <double> k (baselineParamsC[pop][i].size(), 0);
				k[genotypeKeyC[i][j][0]]++;
				k[genotypeKeyC[i][j][1]]++;
				tempBase.push_back(logMultPMF(k, baselineParamsC[pop][i]));
				tempUnsamp.push_back(logMultPMF(k, unsampledPopParamsC[pop][i]));
			}
			lGenos_base.push_back(tempBase);
			lGenos_unsamp.push_back(tempUnsamp);
		}
		

		/*	
		 *	calculate genotype likelihoods for grandparent(maternal) - grandparent(maternal) - grandchild relationship
			given true genotypes
			different for each pop b/c it incorporates allele frequencies for both the baseline pop and the unsampled pop
		 	Index 1 is locus
		 	Index 2 is grandparent1 genotype
			Index 3 is grandparent 2 genotype
			Index 4 is grandchild genotype
			Value is likelihood (NOT log)
		 */
		vector <vector <vector <vector <double> > > > lGenos_ssGP;
		createSSGPvector(genotypeKeyC, lGenos_base, unsampledPopParamsC, pop, lGenos_ssGP);
		
		/*	
		 *	calculate genotype LOG-likelihoods for grandparent(maternal) - grandparent(maternal) - grandchild relationship
		 *	and for all three being unrelated
			given OBSERVED genotypes - but only if NO MISSING genotypes in the trio
		 		Note that this is not the most elegant way to do it, would be better to just calculate with missing as
		 		another "genotype" in this reference "table", but that would require recoding the missing genotype,
		 		(probably as 0), and then all the genotypes are now +1, which means recoding their access, etc, etc.
		 		Hopefully, not much missing data, so this lookup table will be used most of the time. And if there is
		 		very little missing data, this may be faster as it prevents calculating all possible missing combinations
		 		whether they are needed or not.
			different for each pop b/c it incorporates allele frequencies for both the baseline pop and the unsampled pop
		 	Index 1 is locus
		 	Index 2 is grandparent1 genotype
			Index 3 is grandparent 2 genotype
			Index 4 is grandchild genotype
			Value is LOG-likelihood
		 */
		vector <vector <vector <vector <double> > > > lGenos_ssGP_OBS;
		vector <vector <vector <vector <double> > > > lGenos_Unrelated_OBS;
		createOBSvector(lGenos_ssGP, genotypeErrorRatesC, lGenos_ssGP_OBS, 
                  lGenos_base, lGenos_unsamp, lGenos_Unrelated_OBS);

		/* calculate number of MI allowed to consider a pair of potential grandparents
		 * for each locus calculating P(obs, true, MI | true grandparents), then summing to get marginal P(MI | grandparents)
		 * then sum across loci to get P(sum(MI) > threshold | true grandparents)
		 * MIexcludeProb is the target maximum probability that a true trio is removed for having too many MI
		 * miLimit is the maximum number of MI allowed to satisfy MIexcludeProb
		 */
		int miLimit = genotypeKeyC.size();
		if(MIexcludeProb > 0){
			// first calculate the probability of MI at each locus | true grandparents
			vector <double> pMI (genotypeKeyC.size(), 0.0); // p(MI|true grandparents) for each locus
			calcProbMIperLocus(genotypeKeyC, genotypeErrorRatesC, lGenos_ssGP, pMI);
			
			// now calculate the probability of sum(MI) = x | true grandparents
			vector <double> pTotalMI (genotypeKeyC.size() + 1, -9); // probability of having x observed MI | true grandparents
			calcProbSumMI(pMI, pTotalMI);
			
			// now find miLimit such that P(sum(MI) > miLimit) < MIexcludeProb
			double runningSum = 0.0;
			for(int i = 0, max = pTotalMI.size(); i < max; i++){
				runningSum += pTotalMI[i];
				if(runningSum > (1 - MIexcludeProb)){
					miLimit = i;
					break;
				}
			}
			// for(int i=0;i<10;i++) Rcpp::Rcout<<pTotalMI[i]<<"\n"; //testing
			Rcpp::Rcout<<"The maximum number of Mendalian incompatibilities allowed"<<
				" is: "<<miLimit<<". The probability of exclusion for a true grandparent pair (given no missing genotypes) is estimated as: "<<
					1 - runningSum<<".\n";
		} else{
			Rcpp::Rcout<<"MIexcludeProb was not greater than zero, so likelihoods for all possible combinations"<<
				" will be calculated.\n";
		}

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
			listAllPairs(baseInds, pairs);
		}

		Rcpp::Rcout<<"Number of potential grandparent pairs: "<< pairs.size()<<"\n\n";
		// for each individual in the mixture
		for(int m = 0, maxM = mixtureC.size(); m < maxM; m++){
			// cout<<m<<endl; //testing
			for(int i = 0, max = pairs.size(); i < max; i++){ // for each pair
				if(i % 100 == 0) Rcpp::checkUserInterrupt();
				// if(i % 500 == 0) cout<<"pair number: "<<i<<endl;  //testing
				
				/* filter pairs based on MI
				 * for a pair of grandparents, number of MI is the number of loci for which the grandparent pair
				 * and potential descendant share no alleles in common
				 */
				
				int countMI = 0;
				bool skip = false;
				for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ // for each locus
					int gp1Geno = baselineC[pairs[i][0]][j+2];
					int gp2Geno = baselineC[pairs[i][1]][j+2];
					int dGeno = mixtureC[m][j+1];
					if(gp1Geno == -9 || gp2Geno == -9 || dGeno == -9) continue; // skip if any missing
					if(noAllelesInCommonGP(
							genotypeKeyC[j][gp1Geno], // gp1 genotype as vector of alleles
							genotypeKeyC[j][gp2Geno], // gp2 genotype as vector of alleles
							genotypeKeyC[j][dGeno] // o genotype as vector of alleles
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
					// observed genotypes
					int obs_gp1 = baselineC[pairs[i][0]][j+2];
					int obs_gp2 = baselineC[pairs[i][1]][j+2];
					int obs_d = mixtureC[m][j+1];
					
					// if desc missing, no information about the relationship
					if(obs_d == -9) continue;

					// if missing genotypes, information if only one missing OR if an unsampled pop
					//   b/c allele freqs will be different whether its a descendant or not
					if(obs_gp1 == -9 || obs_gp2 == -9){
						// could make lookup tables for these
						double u_likelihood = 0.0; // likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
						double gp_likelihood = 0.0;
						double pOA_gp1, pOA_gp2, pOA_d; // prob of observed given actual
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
									pOA_d = genotypeErrorRatesC[j][d][obs_d];
									
									// unrelated
									// P(gp1|baseline pop)P(gp2|baseline pop)P(d|unsampled pop)P(gp1|obs_gp1)P(gp2|obs_gp2)P(d|obs_d)
									u_likelihood += exp(lGenos_base[j][gp1] + lGenos_base[j][gp2] + lGenos_unsamp[j][d]) * 
										 pOA_gp1 * pOA_gp2 * pOA_d;
									// grandparent pair
									gp_likelihood += lGenos_ssGP[j][gp1][gp2][d] *
										pOA_gp1 * pOA_gp2 * pOA_d;
								}
							}
						}
						uLLH += log(u_likelihood);
						gpLLH += log(gp_likelihood);
					} else { // no missing genotypes
						uLLH += lGenos_Unrelated_OBS[j][obs_gp1][obs_gp2][obs_d];
						gpLLH += lGenos_ssGP_OBS[j][obs_gp1][obs_gp2][obs_d];
					}
				}
				
				// save results
				if((gpLLH - uLLH) >= saveLLR || !filterLLR){
					mixtureInd.push_back(mixtureC[m][0]);
					GrandparentPop.push_back(pop);
					grandparent1.push_back(pairs[i][0]);
					grandparent2.push_back(pairs[i][1]);
					llr.push_back(gpLLH - uLLH);
					mendIncompat.push_back(countMI);
				}

			} // end for each pair

		} // end for each mixture individual
	
	} // end for each baseline pop
	

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
