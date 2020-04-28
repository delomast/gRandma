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
Rcpp::DataFrame falseNeg_ERRORsP(Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::List missingParams,
                           Rcpp::List genotypeKey,
                           Rcpp::List genotypeErrorRates, Rcpp::NumericVector llrToTest,
                           int N,
                           int seed, double MIexcludeProb, int maxMissingGenos
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
	
	// since using true trio as importance sampling distribution, can calculate pos and neg at the same time
	vector <int> popResults;
	vector <double> llrRecord;
	vector <double> falseNegAll;
	vector <double> SDfalseNegAll;
	
	if(MIexcludeProb == 0) Rcpp::Rcout<<"MIexcludeProb was not greater than zero, so only the LLR will be used to"<<
		" accept or reject assignments.\n";
	
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
		
		/* calculate number of MI allowed to consider a pair of potential grandparents
		 * for each locus calculating P(obs, true, MI | true parent), then summing to get marginal P(MI | parent)
		 * then sum across loci to get P(sum(MI) > threshold | true parent)
		 * MIexcludeProb is the target maximum probability that a true trio is removed for having too many MI
		 * miLimit is the maximum number of MI allowed to satisfy MIexcludeProb
		 */
		int miLimit = genotypeKeyC.size();
		if(MIexcludeProb > 0){
			// first calculate the probability of MI at each locus | true parent
			vector <double> pMI (genotypeKeyC.size(), 0.0); // p(MI|true parent) for each locus
			calcProbMIperLocus_sP(genotypeKeyC, genotypeErrorRatesC, lGenos_sP, pMI);
			
			// now calculate the probability of sum(MI) = x
			vector <double> pTotalMI (genotypeKeyC.size() + 1, -9); // probability of having x observed MI | true parent
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
				" is: "<<miLimit<<". The probability of exclusion for a true parent (given no missing genotypes) is estimated as: "<<
					1 - runningSum<<".\n";
		}
		
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
		
		// forward algorithm for missing genotypes
		vector <double> prob_missing_geno (nLoci, 0); // probability that a genotype at locus i is missing
		for(int i = 0; i < nLoci; i++) prob_missing_geno[i] = missingParamsC[i][0] / 
			(missingParamsC[i][0] + missingParamsC[i][1]);
		
		vector <double> probTotalMiss (nLoci + 1, 0); // probability of ending the chain in each state
		vector <vector <double> > probTotalMiss_all; // one entry for each step of the chain
		calcProbSumMI_returnAll(prob_missing_geno, probTotalMiss, probTotalMiss_all);
		
		// set up vectors of states and relative probabilities for sampling
		vector <int> numMissingState (maxMissingGenos + 1);
		vector <double> numMissingProb (maxMissingGenos + 1);
		double tSum = 0;
		for(int i = 0; i <= maxMissingGenos; i++){
			numMissingState[i] = i;
			numMissingProb[i] = probTotalMiss[i];
			tSum += probTotalMiss[i];
		}
		// normalize for kicks
		for(int i = 0; i <= maxMissingGenos; i++) numMissingProb[i] /= tSum;
		
		// some vectors to use with sampleC in the backwards algorithm
		vector <int> select_miss (2, 0);
		select_miss[0] = 1;
		vector <double> prob_miss_nomiss (2, 0);
		
		// inititate results storage for this pop - one entry for each llr
		vector <double> falseNeg (llrToTestC.size(), 0);

		for(int n = 0; n < N; n++){
			if(n % 100 == 0) Rcpp::checkUserInterrupt();
			
			// simulate pair
			vector <int> p1Genos (nLoci, 0);
			vector <int> dGenos (nLoci, 0);
			vector <int> sampPair (2, -1);
			for(int i = 0; i < nLoci; i++){ // for each locus
				// sample observed genotypes
				sampPair = pairGenosLookup[i][sampleC(combo[i], prob[i], rg)]; // sample a trio of observed genotypes

				p1Genos[i] = sampPair[0];
				dGenos[i] = sampPair[1];
			}
			
			// backwards algorithm to add missing genotypes
			// choose state : number of missing genotypes
			int state_p1 = sampleC(numMissingState, numMissingProb, rg);
			int state_d = sampleC(numMissingState, numMissingProb, rg);
			// choose missing genotypes
			for(int i = nLoci; i > 0; i--){
				if(state_p1 + state_d == 0) break;
				//gpa
				prob_miss_nomiss[0] = 0;
				if(state_p1 > 0) prob_miss_nomiss[0] = probTotalMiss_all[i-1][state_p1 - 1] * prob_missing_geno[i-1];
				prob_miss_nomiss[1] = probTotalMiss_all[i-1][state_p1] * (1 - prob_missing_geno[i-1]);
				if(sampleC(select_miss, prob_miss_nomiss, rg) == 1){
					p1Genos[i-1] = -9;
					state_p1--;
				}
				//d
				prob_miss_nomiss[0] = 0;
				if(state_d > 0) prob_miss_nomiss[0] = probTotalMiss_all[i-1][state_d - 1] * prob_missing_geno[i-1];
				prob_miss_nomiss[1] = probTotalMiss_all[i-1][state_d] * (1 - prob_missing_geno[i-1]);
				if(sampleC(select_miss, prob_miss_nomiss, rg) == 1){
					dGenos[i-1] = -9;
					state_d--;
				}
			}
			
			// filter by number of obsereved MI (with missing data), if desired
			if(MIexcludeProb > 0){
				int countMI = 0;
				bool skip = false;
				for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ // for each locus
					int p1Geno = p1Genos[j];
					int dGeno = dGenos[j];
					if(p1Geno == -9 || dGeno == -9) continue; // skip if any missing
					if(noAllelesInCommonSP(
							genotypeKeyC[j][p1Geno], // p1 genotype as vector of alleles
                      genotypeKeyC[j][dGeno] // o genotype as vector of alleles
					)
					) countMI++;
					if (countMI > miLimit){
						skip = true;
						break;
					}
				}
				if(skip){
					// if too many MIs, then this is a false negative for all LLRs
					for(int i = 0, max = llrToTestC.size(); i < max; i++) falseNeg[i]++;
					continue;
				}
			}
			
			// calc LLR
			double uLLH = 0;
			double pLLH = 0;

			for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ //for each locus
				// observed genotypes
				int obs_p1 = p1Genos[j];
				int obs_d = dGenos[j];
				
				if(obs_d == -9) continue; // no information
				
				// if parent is missing, there is information IF there is an "unsampled pop"
				// because allele freqs will be different whether it is an offspring or not
				if(obs_p1 == -9){
					// at some point, can make another lookup table for this
					double u_likelihood = 0; // likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
					double p_likelihood = 0;
					double pOA_d; // prob of observed given actual
					for(int p1 = 0, max3 = genotypeKeyC[j].size(); p1 < max3; p1++){ // for each possible p1 genotype
							for(int d = 0; d < max3; d++){ // for each possible d genotype
								pOA_d = genotypeErrorRatesC[j][d][obs_d];
								// unrelated
								// P(p1|baseline pop)P(d|unsampled pop)P(p1|obs_p1)P(d|obs_d)
								u_likelihood += exp(lGenos_base[j][p1] + lGenos_unsamp[j][d]) * pOA_d;
								// grandparent pair
								p_likelihood += lGenos_sP[j][p1][d] * pOA_d;
							}
						
					}
					uLLH += log(u_likelihood);
					pLLH += log(p_likelihood);
				} else { // no missing genotypes
					uLLH += lGenos_Unrelated_OBS[j][obs_p1][obs_d];
					pLLH += lGenos_sP_OBS[j][obs_p1][obs_d];
				}
			}
			double llr = pLLH - uLLH;
			// determine number under and over threshold
			for(int i = 0, max = llrToTestC.size(); i < max; i++){
				if(llr < llrToTestC[i]) falseNeg[i]++;
			}
		} // end for N
		for(int i = 0, max = llrToTestC.size(); i < max; i++){
			double fnMean = falseNeg[i] / N;

			popResults.push_back(pop);
			llrRecord.push_back(llrToTestC[i]);
			falseNegAll.push_back(fnMean);
			SDfalseNegAll.push_back(sqrt(varianceEstim(N, falseNeg[i], fnMean))); // since all are 1 SS is just sum
		}
	
	} // end for pop
	
	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Pop") = popResults,
                                                   Rcpp::Named("llrThreshold") = llrRecord,
                                                   Rcpp::Named("falseNeg") = falseNegAll,
                                                   Rcpp::Named("falseNegSD") = SDfalseNegAll
																	);
	
	return results;
}
