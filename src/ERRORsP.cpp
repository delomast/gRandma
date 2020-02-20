#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <random>
#include "misc_math.h"
#include "utils.h"
#include "mend_incompat.h"
#include "pairwiseK.h"

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

//' estimatign error rates for parent - offspring pair vs unrelated
//' Monte Carlo used for estimating false negative
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
Rcpp::DataFrame ERRORsP(Rcpp::List baselineParams,
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
	vector <int> popResults;
	vector <double> llrRecord;
	vector <double> falsePosAll;
	vector <double> falseNegAll;
	vector <double> SDfalsePosAll;
	vector <double> SDfalseNegAll;
	vector <double> falsePosAll_aunt;
	vector <double> SDfalsePosAll_aunt;
	vector <double> falsePosAll_halfaunt;
	vector <double> SDfalsePosAll_halfaunt;
	vector <double> falsePosAll_parcous;
	vector <double> SDfalsePosAll_parcous;
	
	// for each baseline population
	for(int pop = 0, max = baselineParamsC.size(); pop < max; pop++){
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
		 * calculate genotype likelihoods for aunt/uncle - niece/nephew
		 * Index 1 is locus, 2 is aunt/uncle geno, 3 is niece/nephew geno
		 * value is likelihood (NOT log)
		 * same for half aunt - niece
		 * same for cousin of true parent - descendant
		 */
		vector <vector <vector <double> > > lGenos_aunt;
		vector <vector <vector <double> > > lGenos_halfaunt;
		vector <vector <vector <double> > > lGenos_parcous;
		vector <double> k_prob(3,0);
		k_prob[0] = 0.5;
		k_prob[1] = 0.5;
		// k_prob[2] = 0;
		pairwiseK_calc2(genotypeKeyC, lGenos_base, baselineParamsC, unsampledPopParamsC, pop,
	                      k_prob, lGenos_aunt);
		k_prob[0] = 0.75;
		k_prob[1] = 0.25;
		// k_prob[2] = 0;
		pairwiseK_calc2(genotypeKeyC, lGenos_base, baselineParamsC, unsampledPopParamsC, pop,
	                      k_prob, lGenos_halfaunt);
		k_prob[0] = 0.875;
		k_prob[1] = 0.125;
		// k_prob[2] = 0;
		
		// testing
		// k_prob[0] = 1;
		// k_prob[1] = 0;
		// end testing
		
		pairwiseK_calc2(genotypeKeyC, lGenos_base, baselineParamsC, unsampledPopParamsC, pop,
	                      k_prob, lGenos_parcous);

		/*	
		 *	calculate genotype LOG-likelihoods for aunt/uncle - niece/nephew relationship
		 *	and for all three being unrelated
			given OBSERVED genotypes - but only if NO MISSING genotypes in the pair
		 	Index 1 is locus, 2 is aunt/uncle geno, 3 is niece/nephew geno
			Value is LOG-likelihood
		 * same for half aunt - niece
		 * same for cousin of true parent - descendant
		 */
		vector <vector <vector <double> > > lGenos_aunt_OBS;
		vector <vector <vector <double> > > lGenos_halfaunt_OBS;
		vector <vector <vector <double> > > lGenos_parcous_OBS;
		create_pairswise_OBSvector(lGenos_aunt, genotypeErrorRatesC, lGenos_aunt_OBS);
		create_pairswise_OBSvector(lGenos_halfaunt, genotypeErrorRatesC, lGenos_halfaunt_OBS);
		create_pairswise_OBSvector(lGenos_parcous, genotypeErrorRatesC, lGenos_parcous_OBS);
		
		/*	
		 *	calculate genotype likelihoods for parent - offspring relationship
			given true genotypes
			different for each pop b/c it incorporates allele frequencies for both the baseline pop and the unsampled pop
		 	Index 1 is locus
		 	Index 2 is parent genotype
			Index 3 is offspring genotype
			Value is likelihood (NOT log)
		 */
		vector <vector <vector <double> > > lGenos_sP;
		createSPvector(genotypeKeyC, lGenos_base, unsampledPopParamsC, pop, lGenos_sP);

		// testing
		// for(int t = 0; t < genotypeKeyC[0].size(); t++){
		// 	for(int t2 = 0; t2 < genotypeKeyC[0].size(); t2++){
		// 		// Rcpp::Rcout << lGenos_sP[0][t][t2] << " " << lGenos_parcous[0][t][t2] << "\n";
		// 		Rcpp::Rcout << exp(lGenos_base[0][t] + lGenos_base[0][t2]) << " " << lGenos_parcous[0][t][t2] << "\n";
		// 	}
		// }
		// Rcpp::stop("here");
		// end testing
		
		/*	
		 *	calculate genotype LOG-likelihoods for parent - offspring relationship
		 *	and for all three being unrelated
			given OBSERVED genotypes - but only if NO MISSING genotypes in the pair
			different for each pop b/c it incorporates allele frequencies for both the baseline pop and the unsampled pop
		 	Index 1 is locus
		 	Index 2 is parent1 genotype
			Index 3 is offspring genotype
			Value is LOG-likelihood
		 */
		vector <vector <vector <double> > > lGenos_sP_OBS;
		vector <vector <vector <double> > > lGenos_Unrelated_OBS;
		createSP_OBSvector(lGenos_sP, genotypeErrorRatesC, lGenos_sP_OBS, 
                  lGenos_base, lGenos_unsamp, lGenos_Unrelated_OBS);
		
		/*
		 * rearranging lGenos_ssGP_OBS to be easier to use for sampling pair genotypes
		 * So for each locus, combo has the index of each trio Genotype combination
		 * probability has the probability of each combo
		 * trioGenosLookup has a vector giving the two genotypes in each combo (p1, d)
		 */
		vector <vector <int> > combo; // index of the genotype combos
		vector <vector <double> > prob; // the probabilities of each genotype combo
		vector <vector <vector <int> > > pairGenosLookup; // the genotype combos
		for(int i = 0; i < nLoci; i++){
			vector <int> tempCombo;
			vector <double> tempProb;
			vector <vector <int> > tempGenos;
			int pairIndex = 0;
			for(int p1 = 0, max = lGenos_sP_OBS[i].size(); p1 < max; p1++){
				for(int d = 0; d < max; d++){
					tempProb.push_back(exp(lGenos_sP_OBS[i][p1][d]));
					tempCombo.push_back(pairIndex);
					vector <int> tempVec (2,-9);
					tempVec[0] = p1;
					tempVec[1] = d;
					tempGenos.push_back(tempVec);
					pairIndex++;
				}
			}
			combo.push_back(tempCombo);
			prob.push_back(tempProb);
			pairGenosLookup.push_back(tempGenos);
		}
		
		// inititate results storage for this pop - one entry for each llr
		vector <double> falsePos (llrToTestC.size(), 0);
		vector <double> falseNeg (llrToTestC.size(), 0);
		vector <double> SSfalsePos (llrToTestC.size(), 0); // sum of squares for variance calcs
		vector <double> falsePos_aunt (llrToTestC.size(), 0);
		vector <double> SSfalsePos_aunt (llrToTestC.size(), 0);
		vector <double> falsePos_halfaunt (llrToTestC.size(), 0);
		vector <double> SSfalsePos_halfaunt (llrToTestC.size(), 0);
		vector <double> falsePos_parcous (llrToTestC.size(), 0);
		vector <double> SSfalsePos_parcous (llrToTestC.size(), 0);

		// don't need sum squares for falseNeg right now b/c all are either 0 or 1 (no importance sampling)
		
		for(int n = 0; n < N; n++){
			if(n % 100 == 0) Rcpp::checkUserInterrupt();
			
			// simulate pair
			vector <int> p1Genos (nLoci, 0);
			vector <int> dGenos (nLoci, 0);
			vector <int> sampPair (2, -1);
			for(int i = 0; i < nLoci; i++){ // for each locus
				// sample observed genotypes
				sampPair = pairGenosLookup[i][sampleC(combo[i], prob[i], rg)]; // sample a trio of observed genotypes

				// add missing genotypes as appropriate
				for(int j = 0; j < 2; j++) if(randMissing[i](rg)) sampPair[j] = -9;
				p1Genos[i] = sampPair[0];
				dGenos[i] = sampPair[1];
			}
			
			// calc LLR
			double uLLH = 0;
			double pLLH = 0;
			double aunt_LLH = 0; // correction for sampling distribution
			double halfaunt_LLH = 0;
			double parcous_LLH = 0;
			for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ //for each locus
				// observed genotypes
				int obs_p1 = p1Genos[j];
				int obs_d = dGenos[j];
				
				// if missing genotypes
				if(obs_p1 == -9 || obs_d == -9){
					double u_likelihood = 0; // likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
					double p_likelihood = 0;
					double aunt_likelihood = 0;
					double halfaunt_likelihood = 0;
					double parcous_likelihood = 0;
					double pOA_p1, pOA_d; // prob of observed given actual
					for(int p1 = 0, max3 = genotypeKeyC[j].size(); p1 < max3; p1++){ // for each possible p1 genotype
							for(int d = 0; d < max3; d++){ // for each possible d genotype
								if(obs_p1 == -9){
									pOA_p1 = 1;
								} else {
									pOA_p1 = genotypeErrorRatesC[j][p1][obs_p1];
								}
								if(obs_d == -9){
									pOA_d = 1;
								} else {
									pOA_d = genotypeErrorRatesC[j][d][obs_d];
								}
								// unrelated
								// P(p1|baseline pop)P(d|unsampled pop)P(p1|obs_p1)P(d|obs_d)
								u_likelihood += exp(lGenos_base[j][p1] + lGenos_unsamp[j][d]) * 
									 pOA_p1 * pOA_d;
								// grandparent pair
								p_likelihood += lGenos_sP[j][p1][d] *
									pOA_p1 * pOA_d;
								// aunt
								aunt_likelihood += lGenos_aunt[j][p1][d] *
									pOA_p1 * pOA_d;
								// half aunt
								halfaunt_likelihood += lGenos_halfaunt[j][p1][d] *
									pOA_p1 * pOA_d;
								// cousing of true parent
								parcous_likelihood += lGenos_parcous[j][p1][d] *
									pOA_p1 * pOA_d;
							}
						
					}
					uLLH += log(u_likelihood);
					pLLH += log(p_likelihood);
					aunt_LLH += log(aunt_likelihood);
					halfaunt_LLH += log(halfaunt_likelihood);
					parcous_LLH += log(parcous_likelihood);
				} else { // no missing genotypes
					uLLH += lGenos_Unrelated_OBS[j][obs_p1][obs_d];
					pLLH += lGenos_sP_OBS[j][obs_p1][obs_d];
					aunt_LLH += lGenos_aunt_OBS[j][obs_p1][obs_d];
					halfaunt_LLH += lGenos_halfaunt_OBS[j][obs_p1][obs_d];
					parcous_LLH += lGenos_parcous_OBS[j][obs_p1][obs_d];
				}
			}
			double llr = pLLH - uLLH;
			// determine number under and over threshold
			for(int i = 0, max = llrToTestC.size(); i < max; i++){
				if(llr < llrToTestC[i]){
					falseNeg[i]++;
				} else {
					falsePos[i] += exp(-llr); // correction for sampling distribution
					SSfalsePos[i] += pow(exp(-llr), 2);
					falsePos_aunt[i] += exp(aunt_LLH - pLLH);
					SSfalsePos_aunt[i] += pow(exp(aunt_LLH - pLLH), 2);
					falsePos_halfaunt[i] += exp(halfaunt_LLH - pLLH);
					SSfalsePos_halfaunt[i] += pow(exp(halfaunt_LLH - pLLH), 2);
					falsePos_parcous[i] += exp(parcous_LLH - pLLH);
					SSfalsePos_parcous[i] += pow(exp(parcous_LLH - pLLH), 2);
				}
			}
		} // end for N
		for(int i = 0, max = llrToTestC.size(); i < max; i++){
			double fpMean = falsePos[i] / N; // need for variance calc
			double fnMean = falseNeg[i] / N;
			double fpMean_aunt = falsePos_aunt[i] / N;
			double fpMean_halfaunt = falsePos_halfaunt[i] / N;
			double fpMean_parcous = falsePos_parcous[i] / N;

			popResults.push_back(pop);
			llrRecord.push_back(llrToTestC[i]);
			falsePosAll.push_back(fpMean);
			falseNegAll.push_back(fnMean);
			SDfalsePosAll.push_back(sqrt(varianceEstim(N, SSfalsePos[i], fpMean)));
			SDfalseNegAll.push_back(sqrt(varianceEstim(N, falseNeg[i], fnMean))); // since all are 1 SS is just sum
			falsePosAll_aunt.push_back(fpMean_aunt);
			SDfalsePosAll_aunt.push_back(sqrt(varianceEstim(N, SSfalsePos_aunt[i], fpMean_aunt)));
			falsePosAll_halfaunt.push_back(fpMean_halfaunt);
			SDfalsePosAll_halfaunt.push_back(sqrt(varianceEstim(N, SSfalsePos_halfaunt[i], fpMean_halfaunt)));
			falsePosAll_parcous.push_back(fpMean_parcous);
			SDfalsePosAll_parcous.push_back(sqrt(varianceEstim(N, SSfalsePos_parcous[i], fpMean_parcous)));
		}
	
	} // end for pop
	
	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Pop") = popResults,
                                                   Rcpp::Named("llrThreshold") = llrRecord,
                                                   Rcpp::Named("falsePos") = falsePosAll,
                                                   Rcpp::Named("falsePosSD") = SDfalsePosAll,
                                                   Rcpp::Named("falseNeg") = falseNegAll,
                                                   Rcpp::Named("falseNegSD") = SDfalseNegAll,
                                                   Rcpp::Named("falsePosAunt") = falsePosAll_aunt,
                                                   Rcpp::Named("falsePosAuntSD") = SDfalsePosAll_aunt,
                                                   Rcpp::Named("falsePosHalfAunt") = falsePosAll_halfaunt,
                                                   Rcpp::Named("falsePosHalfAuntSD") = SDfalsePosAll_halfaunt,
                                                   Rcpp::Named("falsePosParCous") = falsePosAll_parcous,
                                                   Rcpp::Named("falsePosParCousSD") = SDfalsePosAll_parcous
																	);
	
	return results;
}
