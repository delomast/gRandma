#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <random>
#include "misc_math.h"
#include "utils.h"
#include "mend_incompat.h"

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

//' estimating false positive error rates for single-sided grandparent pair vs unrelated
//' with individuls from one population being assigned to a different population
//' Importance sampling Monte Carlo used for estimating false positive
//' @param baselineParams Dirichlet parameters for allele frequencies
//' @param unsampledPopParams Dirichlet parameters for allele frequencies
//' @param missingParams Beta parameters for missing genotypes (failure to genotype rate)
//' @param genotypeKey list of matrix for each locus, col1 is genotype, col2 is allele 1, col3 is allele 2
//' @param genotypeErrorRates list of matrix for each locus, rows are actual genotype, columns are observed,
//'   values are probability
//' @param llrToTest Vector of llr's to test as threshold values
//' @param N number of samples to take
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame otherPopERRORssGP(Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::List missingParams,
                           Rcpp::List genotypeKey,
                           Rcpp::List genotypeErrorRates, Rcpp::NumericVector llrToTest,
                           int N,
                           int seed
){
	//////////////////////
	// first, turn all inputs into non-Rcpp c++ classes

	vector <vector <vector <double> > > baselineParamsC; //first index is pop, then locus, then allele
	for(int i = 0, max = baselineParams.length(); i < max; i++){
		baselineParamsC.push_back(rcppListToVectorDouble(baselineParams[i]));
	}

	vector <vector <vector <double> > > unsampledPopParamsC; //first index is pop, then locus, then allele
	for(int i = 0, max = unsampledPopParams.length(); i < max; i++){
		unsampledPopParamsC.push_back(rcppListToVectorDouble(unsampledPopParams[i]));
	}

	vector <vector <double> > missingParamsC; // first index is locus, then missing and not missing parameter
	missingParamsC = rcppListToVectorDouble(missingParams);

	vector <vector <vector <int> > > genotypeKeyC; //first index is locus, then genotype, values are alleles
	genotypeKeyC = convertGenotypeKey(genotypeKey);

	vector< vector < vector <double> > > genotypeErrorRatesC; // first index is locus, then actual genotype, then observed
	for(int i = 0, max = genotypeErrorRates.length(); i < max; i++){
		genotypeErrorRatesC.push_back(rcppMatrixToVectorDouble(genotypeErrorRates[i]));
	}

	vector <double> llrToTestC;
	for(int i = 0, max = llrToTest.length(); i < max; i++) llrToTestC.push_back(llrToTest[i]);
	
	// end conversion of types
	////////////////////

	int nLoci = genotypeKeyC.size();
	
	// initiate random number generator
	mt19937 rg (seed);

	// random number generators for missing genotypes for each locus
	// simple Bernoulli
	vector <bernoulli_distribution> randMissing;
	for(int i = 0; i < nLoci; i++){ 
		bernoulli_distribution tempDist (missingParamsC[i][0] / (missingParamsC[i][0] + missingParamsC[i][1]));
		randMissing.push_back(tempDist);
	}
	
	// since using true trio as importance sampling distribution, can calculate pos and neg at the same time
	vector <int> popResults; // pop that trios are simulated from
	vector <int> popAssign; // pop that offspring are assigned to
	vector <double> llrRecord;
	vector <double> falsePosAll;
	vector <double> SDfalsePosAll;

	// for each baseline population
	for(int pop = 0, max = baselineParamsC.size(); pop < max; pop++){
		Rcpp::Rcout<<"Now considering descendants from baseline population "<<pop + 1<<"\n";
		
		// calculate genotype log-likelihoods from random sampled genotype
		// first index is locus, second is genotype, value is LOG-likelihood
		vector <vector <double> > CORR_lGenos_base;
		vector <vector <double> > CORR_lGenos_unsamp;

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
				tempBase.push_back(logDirichMultPMF(k, baselineParamsC[pop][i]));
				tempUnsamp.push_back(logDirichMultPMF(k, unsampledPopParamsC[pop][i]));
			}
			CORR_lGenos_base.push_back(tempBase);
			CORR_lGenos_unsamp.push_back(tempUnsamp);
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
		// only need this for createOBSvector below, could separate createOBSvector into two functions
		// if speed become an issue and this seems to be slowing things down
		vector <vector <vector <vector <double> > > > CORR_lGenos_ssGP;
		createSSGPvector(genotypeKeyC, CORR_lGenos_base, unsampledPopParamsC, pop, CORR_lGenos_ssGP);
		
		/*	
		 *	calculate genotype LOG-likelihoods for grandparent(maternal) - grandparent(maternal) - grandchild relationship
		 *	and for all three being unrelated
			given OBSERVED genotypes - but only if NO MISSING genotypes in the trio
			different for each pop b/c it incorporates allele frequencies for both the baseline pop and the unsampled pop
		 	Index 1 is locus
		 	Index 2 is grandparent1 genotype
			Index 3 is grandparent 2 genotype
			Index 4 is grandchild genotype
			Value is LOG-likelihood
		 */
		vector <vector <vector <vector <double> > > > CORR_lGenos_ssGP_OBS;
		vector <vector <vector <vector <double> > > > CORR_lGenos_Unrelated_OBS;
		createOBSvector(CORR_lGenos_ssGP, genotypeErrorRatesC, CORR_lGenos_ssGP_OBS,
                  CORR_lGenos_base, CORR_lGenos_unsamp, CORR_lGenos_Unrelated_OBS);
		
		for(int pop2 = 0, max = baselineParamsC.size(); pop2 < max; pop2++){
			if(pop2 == pop) continue;
			Rcpp::Rcout<<"Calculating error with baseline population "<<pop2 + 1<<"\n";
			// calculate genotype log-likelihoods from random sampled genotype
			// first index is locus, second is genotype, value is LOG-likelihood
			vector <vector <double> > lGenos_base;
			vector <vector <double> > lGenos_unsamp;
	
			for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ // for each locus
				// baselineParamsC[pop2][i] are the alpha values for the baseline pop
				// unsampledPopParamsC[pop2][i] are the alpha values for the corresponding unsampled pop
				// genotypeKeyC[i] is the key for the locus i
				vector <double> tempBase;
				vector <double> tempUnsamp;
				for(int j = 0, max2 = genotypeKeyC[i].size(); j < max2; j++){ // for each genotype
					vector <double> k (baselineParamsC[pop2][i].size(), 0);
					k[genotypeKeyC[i][j][0]]++;
					k[genotypeKeyC[i][j][1]]++;
					tempBase.push_back(logDirichMultPMF(k, baselineParamsC[pop2][i]));
					tempUnsamp.push_back(logDirichMultPMF(k, unsampledPopParamsC[pop2][i]));
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
			createSSGPvector(genotypeKeyC, lGenos_base, unsampledPopParamsC, pop2, lGenos_ssGP);
			
			/*	
			 *	calculate genotype LOG-likelihoods for grandparent(maternal) - grandparent(maternal) - grandchild relationship
			 *	and for all three being unrelated
				given OBSERVED genotypes - but only if NO MISSING genotypes in the trio
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
			
			/*
			 * rearranging lGenos_ssGP_OBS to be easier to use for sampling trio genotypes
			 * So for each locus, combo has the index of each trio Genotype combination
			 * probability has the probability of each combo
			 * trioGenosLookup has a vector giving the three genotypes in each combo (gpa, gma, d)
			 */
			vector <vector <int> > combo; // index of the genotype combos
			vector <vector <double> > prob; // the probabilities of each genotype combo
			vector <vector <vector <int> > > trioGenosLookup; // the genotype combos
			for(int i = 0; i < nLoci; i++){
				vector <int> tempCombo;
				vector <double> tempProb;
				vector <vector <int> > tempGenos;
				int trioIndex = 0;
				for(int gpa = 0, max = lGenos_ssGP_OBS[i].size(); gpa < max; gpa++){
					for(int gma = 0; gma < max; gma++){
						for(int d = 0; d < max; d++){
							tempProb.push_back(exp(lGenos_ssGP_OBS[i][gpa][gma][d]));
							tempCombo.push_back(trioIndex);
							vector <int> tempVec (3,-9);
							tempVec[0] = gpa;
							tempVec[1] = gma;
							tempVec[2] = d;
							tempGenos.push_back(tempVec);
							trioIndex++;
						}
					}
				}
				combo.push_back(tempCombo);
				prob.push_back(tempProb);
				trioGenosLookup.push_back(tempGenos);
			}
			
			// inititate results storage for this pop - one entry for each llr
			vector <double> falsePos (llrToTestC.size(), 0);
			vector <double> SSfalsePos (llrToTestC.size(), 0); // sum of squares for variance calcs

			for(int n = 0; n < N; n++){
				if(n % 100 == 0) Rcpp::checkUserInterrupt();
				
				// simulate trio
				vector <int> gmaGenos (nLoci, 0);
				vector <int> gpaGenos (nLoci, 0);
				vector <int> dGenos (nLoci, 0);
				vector <int> sampTrio (3, -1);
				for(int i = 0; i < nLoci; i++){ // for each locus
					// sample observed genotypes
					sampTrio = trioGenosLookup[i][sampleC(combo[i], prob[i], rg)]; // sample a trio of observed genotypes
	
					// add missing genotypes as appropriate
					for(int j = 0; j < 3; j++) if(randMissing[i](rg)) sampTrio[j] = -9;
					gpaGenos[i] = sampTrio[0];
					gmaGenos[i] = sampTrio[1];
					dGenos[i] = sampTrio[2];
				}
				
				// calc LLR
				double uLLH = 0.0;
				double gpLLH = 0.0;
				double CORR_uLLH = 0.0; // for correction ratio for importance sampling
				for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ //for each locus
					// observed genotypes
					int obs_gp1 = gpaGenos[j];
					int obs_gp2 = gmaGenos[j];
					int obs_d = dGenos[j];
	
					// if missing genotypes
					if(obs_gp1 == -9 || obs_gp2 == -9 || obs_d == -9){
						double u_likelihood = 0.0; // likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
						double gp_likelihood = 0.0;
						double CORR_u_likelihood = 0.0;
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
									if(obs_d == -9){
										pOA_d = 1.0;
									} else {
										pOA_d = genotypeErrorRatesC[j][d][obs_d];
									}
									// unrelated
									// P(gp1|baseline pop)P(gp2|baseline pop)P(d|unsampled pop)P(gp1|obs_gp1)P(gp2|obs_gp2)P(d|obs_d)
									u_likelihood += exp(lGenos_base[j][gp1] + lGenos_base[j][gp2] + lGenos_unsamp[j][d]) * 
										 pOA_gp1 * pOA_gp2 * pOA_d;
									// grandparent pair
									gp_likelihood += lGenos_ssGP[j][gp1][gp2][d] *
										pOA_gp1 * pOA_gp2 * pOA_d;
									// correction for importance sampling
									CORR_u_likelihood += exp(CORR_lGenos_base[j][gp1] + CORR_lGenos_base[j][gp2] + CORR_lGenos_unsamp[j][d]) * 
										 pOA_gp1 * pOA_gp2 * pOA_d;
								}
							}
						}
						uLLH += log(u_likelihood);
						gpLLH += log(gp_likelihood);
						CORR_uLLH += log(CORR_u_likelihood);
					} else { // no missing genotypes
						uLLH += lGenos_Unrelated_OBS[j][obs_gp1][obs_gp2][obs_d];
						gpLLH += lGenos_ssGP_OBS[j][obs_gp1][obs_gp2][obs_d];
						CORR_uLLH += CORR_lGenos_Unrelated_OBS[j][obs_gp1][obs_gp2][obs_d];
					}
				}
				double llr = gpLLH - uLLH;
				// determine number under and over threshold
				for(int i = 0, max = llrToTestC.size(); i < max; i++){
					if(llr >= llrToTestC[i]){
						falsePos[i] += exp(CORR_uLLH - gpLLH); // correction for sampling distribution
						SSfalsePos[i] += pow(exp(-llr), 2);
					}
				}
			} // end for N
			for(int i = 0, max = llrToTestC.size(); i < max; i++){
				double fpMean = falsePos[i] / N; // need for variance calc

				popResults.push_back(pop);
				popAssign.push_back(pop2);
				llrRecord.push_back(llrToTestC[i]);
				falsePosAll.push_back(fpMean);
				SDfalsePosAll.push_back(sqrt(varianceEstim(N, SSfalsePos[i], fpMean)));
			}
		
		} // end for pop2
			
	} // end for pop

	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Pop_offspring") = popResults,
                                                   Rcpp::Named("Pop_baseline") = popAssign,
                                                   Rcpp::Named("llrThreshold") = llrRecord,
                                                   Rcpp::Named("falsePos") = falsePosAll,
                                                   Rcpp::Named("falsePosSD") = SDfalsePosAll
																	);
	
	return results;
}
