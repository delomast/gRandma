#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <random>
#include <string>
#include <limits>
#include "misc_math.h"
#include "utils.h"
#include "mend_incompat.h"
#include "pairwiseK.h"

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

//' estimating positive error rates for parent - offspring pairs
//' stratified Monte Carlo used for estimating false positive
//' @param baselineParams Dirichlet parameters for allele frequencies
//' @param unsampledPopParams Dirichlet parameters for allele frequencies
//' @param missingParams Beta parameters for missing genotypes (failure to genotype rate)
//' @param genotypeKey list of matrix for each locus, col1 is genotype, col2 is allele 1, col3 is allele 2
//' @param genotypeErrorRates list of matrix for each locus, rows are actual genotype, columns are observed,
//'   values are probability
//' @param llrToTest Vector of llr's to test as threshold values
//' @param N number of samples to take
//' @param trueRel an integer indicating the true relationship to simulate: 0 = unrelated, 1 = aunt,
//'   2 = half aunt, 3 = cousin of parent
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::List strat_ERRORsP(Rcpp::List baselineParams,
                         Rcpp::List unsampledPopParams, Rcpp::List missingParams,
                         Rcpp::List genotypeKey,
                         Rcpp::List genotypeErrorRates, Rcpp::NumericVector llrToTest,
                         Rcpp::NumericVector itersPerMI,
                         int seed, int trueRel, double MIexcludeProb, 
                         double maxMissingGenos){
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
	
	vector <double> itersPerMIC;
	for(int i = 0, max = itersPerMI.length(); i < max; i++) itersPerMIC.push_back(itersPerMI[i]);
	
	vector <double> llrToTestC;
	for(int i = 0, max = llrToTest.length(); i < max; i++) llrToTestC.push_back(llrToTest[i]);
	
	// end conversion of types
	////////////////////
	
	int nLoci = genotypeKeyC.size();
	
	// initiate random number generator
	mt19937 rg (seed);
	
	vector <int> popResults;
	vector <double> llrRecord;
	vector <double> falsePosAll;
	vector <double> SDfalsePosAll;
	
	// for output of full results
	Rcpp::List fullOut;
	
	// choose values depending of true relationship being simulated
	vector <double> k_prob(3,0);
	string fpColName ("falsePos");
	string fpSDColName;
	if(trueRel == 0){
		// Unrelated
		k_prob[0] = 1;
		k_prob[1] = 0;
		// k_prob[2] = 0;
		fpColName += "Unrel";
	} else if (trueRel == 1){
		// Aunt
		k_prob[0] = 0.5;
		k_prob[1] = 0.5;
		// k_prob[2] = 0;
		fpColName += "Aunt";
	} else if (trueRel == 2){
		// half aunt
		k_prob[0] = 0.75;
		k_prob[1] = 0.25;
		// k_prob[2] = 0;
		fpColName += "HalfAunt";
	} else if (trueRel == 3){
		// cousin of parent
		k_prob[0] = 0.875;
		k_prob[1] = 0.125;
		// k_prob[2] = 0;
		fpColName += "ParCous";
	} else {
		Rcpp::stop("Internal error: trueRel not recognized");
	}
	fpSDColName = fpColName + "SD";
	
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
		* calculate genotype likelihoods for relationship (example here is aunt/uncle - niece/nephew)
		* Index 1 is locus, 2 is aunt/uncle geno, 3 is niece/nephew geno
		* value is likelihood (NOT log)
		*/
		vector <vector <vector <double> > > lGenos_tRel;
		pairwiseK_calc2(genotypeKeyC, lGenos_base, baselineParamsC, unsampledPopParamsC, pop,
                  k_prob, lGenos_tRel);
		
		/*	
		*	calculate genotype LOG-likelihoods for (example aunt/uncle - niece/nephew) relationship
		*	and for all three being unrelated
		given OBSERVED genotypes - but only if NO MISSING genotypes in the pair
		Index 1 is locus, 2 is aunt/uncle geno, 3 is niece/nephew geno
		Value is LOG-likelihood
		*/
		vector <vector <vector <double> > > lGenos_tRel_OBS;
		create_pairwise_OBSvector(lGenos_tRel, genotypeErrorRatesC, lGenos_tRel_OBS);
		
		
		/*	Note that this could be done with pairwise K functions, but is different for historical reasons
		* 	ie I wrote it before writing the pairwise K functions!
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
		*	and for being unrelated
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
		
		int miLimit = genotypeKeyC.size();
		if(MIexcludeProb > 0){
			// first calculate the probability of MI at each locus | true parent
			vector <double> pMI_sP (genotypeKeyC.size(), 0.0); // p(MI|true parent) for each locus
			calcProbMIperLocus_sP(genotypeKeyC, genotypeErrorRatesC, lGenos_sP, pMI_sP);
			
			// now calculate the probability of sum(MI) = x
			vector <double> pTotalMI_sP (genotypeKeyC.size() + 1, -9); // probability of having x observed MI | true parent
			calcProbSumMI(pMI_sP, pTotalMI_sP);
			
			// now find miLimit such that P(sum(MI) > miLimit) < MIexcludeProb
			double runningSum = 0.0;
			for(int i = 0, max = pTotalMI_sP.size(); i < max; i++){
				runningSum += pTotalMI_sP[i];
				if(runningSum > (1 - MIexcludeProb)){
					miLimit = i;
					break;
				}
			}
			Rcpp::Rcout<<"The maximum number of Mendalian incompatibilities allowed"<<
				" is: "<<miLimit<<". The probability of exclusion for a true parent (given no missing genotypes) is estimated as: "<<
					1 - runningSum<<".\n";
		}
		
		if(static_cast<int>(itersPerMIC.size()) < (miLimit + 1)) Rcpp::stop("Error: the length of itersPerMI is not long enough for the number of allowed MIs in this population.");
		if((miLimit + 1) * pow(maxMissingGenos + 1, 2) * nLoci > 345000000) Rcpp::Rcout << "\nUsing a large amount of memory...\n\n";

		vector <vector <int> > combo; // index of the genotype combos
		vector <vector <double> > prob_MI; // the probabilities of each genotype combo
		vector <vector <double> > prob_no_MI; // the probabilities of each genotype combo
		vector <vector <double> > prob_all; // the probabilities of each genotype combo
		vector <vector <vector <int> > > pairGenosLookup; // the genotype combos
		
		{ // scope to conserve some memory
			/*
			 * calculate probability of pairs of OBSERVED genotypes given MI at a locus or given no MI at a locus 
			 * indices are: locus, baseline genotype, desc genotype
			 * value is probability NOT log
			 */
			vector <vector <vector <double> > > lGenos_OBS_sample_MI;
			vector <vector <vector <double> > > lGenos_OBS_sample_no_MI;
			strat_lGenosBuilder_Pairwise(genotypeKeyC, genotypeErrorRatesC, lGenos_tRel_OBS,
	                               lGenos_OBS_sample_MI,
	                               lGenos_OBS_sample_no_MI);
			/*
			 * now making lGenos_OBS_sample_(MI / no_MI) easier to sample from
			 * So for each locus, combo has the index of each trio Genotype combination
			 * probability has the probability of each combo
			 * pairGenosLookup has a vector giving the two genotypes in each combo (p1, d)
			 */
			for(int i = 0; i < nLoci; i++){
				vector <int> tempCombo;
				double MIsum = 0;
				double noMIsum = 0;
				vector <double> tempProb_MI;
				vector <double> tempProb_no_MI;
				vector <double> tempProb_all;
				vector <vector <int> > tempGenos;
				int pairIndex = 0;
				for(int p1 = 0, max = lGenos_OBS_sample_MI[i].size(); p1 < max; p1++){
					for(int d = 0; d < max; d++){
						tempProb_MI.push_back(lGenos_OBS_sample_MI[i][p1][d]);
						tempProb_no_MI.push_back(lGenos_OBS_sample_no_MI[i][p1][d]);
						// one is zero and the other isn't, prob of either is just sum, disjoint
						tempProb_all.push_back(lGenos_OBS_sample_MI[i][p1][d] + lGenos_OBS_sample_no_MI[i][p1][d]);
						MIsum += lGenos_OBS_sample_MI[i][p1][d];
						noMIsum += lGenos_OBS_sample_no_MI[i][p1][d];
						tempCombo.push_back(pairIndex);
						vector <int> tempVec (2,-9);
						tempVec[0] = p1;
						tempVec[1] = d;
						tempGenos.push_back(tempVec);
						pairIndex++;
					}
				}
				// normalize
				for(int j = 0, max = tempCombo.size(); j < max; j++){
					tempProb_MI[j] /= MIsum;
					tempProb_no_MI[j] /= noMIsum;
					tempProb_all[j] /= (MIsum + noMIsum);
				}
				combo.push_back(tempCombo);
				prob_MI.push_back(tempProb_MI);
				prob_no_MI.push_back(tempProb_no_MI);
				prob_all.push_back(tempProb_all);
				pairGenosLookup.push_back(tempGenos);
			}
		} // end scoping
		
		// forward algorithm to calculate probabilities of all 
		//   states (#MI, #miss poten p1, #miss d)
		
		vector <double> prob_missing_geno (nLoci, 0); // probability that a genotype at locus i is missing
		for(int i = 0; i < nLoci; i++) prob_missing_geno[i] = missingParamsC[i][0] / 
			(missingParamsC[i][0] + missingParamsC[i][1]);
		// calculate probability of observing an MI at each locus given that all genotypes are present
		vector <double> pMI_noMiss (genotypeKeyC.size(), 0.0); // p(MI) observed for each locus
		calcProbMIperLocus_Pair(genotypeKeyC, genotypeErrorRatesC, lGenos_tRel, pMI_noMiss);
		
		// save the states vector for each step
		vector <vector <vector <vector <double> > > > all_states;
		fwd_sP(miLimit, maxMissingGenos, nLoci, prob_missing_geno, pMI_noMiss, all_states);
		
		// probability of an individual having more than the allowed number of missing genotypes
		double probFail = 1;
		{
			vector <double> probTotalFail (nLoci + 1, 0);
			// can use this even though it's not MIs
			calcProbSumMI(prob_missing_geno, probTotalFail);
			for(int i = 0; i <= maxMissingGenos; i++) probFail -= probTotalFail[i];
		}
		// probability of one or more individuals having more than the allowed number of missing genotypes
		double probAnyFail = 1 - pow(1 - probFail, 2); // just complement of both succeeding
		
		// calculate probability of being in each strata, given maximum number of missing genotypes per individual
		vector <double> pNumMI_OBS (miLimit + 1, 0);
		for(int mi = 0; mi <= miLimit; mi++){ // #mi
			for(int jp1 = 0; jp1 <= maxMissingGenos ; jp1++){ // #miss p1
				for(int jd = 0; jd <= maxMissingGenos; jd++){ // #miss d
					pNumMI_OBS[mi] += all_states[nLoci][mi][jp1][jd];
				}
			}
		}
		// probability of each strata, given that none have too many missing genotypes
		for(int mi = 0; mi <= miLimit; mi++) pNumMI_OBS[mi] /= (1 - probAnyFail);
		
		// some results storage for full results reporting
		vector <vector <double> > falsePosPerStrata (miLimit + 1, 
              vector <double> (llrToTestC.size(), 0)); // indices: num of MI, llr threshold
		
		// inititate results storage for this pop - one entry for each llr
		vector <double> fp_SD_total (llrToTestC.size(), 0);
		
		// two vectors to use with sampleC below
		vector <double> probs_MI_NoMI (2, 0);
		vector <int> select_MI (2);
		select_MI[0] = 1;
		select_MI[1] = 0;
		
		for(int mi = 0; mi <= miLimit; mi++){
			Rcpp::checkUserInterrupt();
			// check if prob of seeing this many MIs is zero
			// if it's zero, throw a warning and skip this stratum
			if (pNumMI_OBS[mi] <= numeric_limits<double>::min() && itersPerMIC[mi] > 0){
				string errMess = "The probability of observing ";
				errMess += to_string(mi); 
				errMess += " Mendelian incompatibilities for this relationship is 0. "; 
				errMess += "Skipping all iterations in this stratum."; 
				Rcpp::warning(errMess);
				continue;
			}
			
			// set up to sample from states within this strata
			vector <int> stateIndex;
			vector <double> stateProb;
			vector <vector <int> > statesToSample;
			int iState = 0;
			double tempSum = 0;
			for(int jp1 = 0; jp1 <= maxMissingGenos ; jp1++){ // #miss p1
				for(int jd = 0; jd <= maxMissingGenos; jd++){ // #miss d
					vector <int> tempStateVec (3);
					tempStateVec[0] = mi;
					tempStateVec[1] = jp1;
					tempStateVec[2] = jd;
					stateIndex.push_back(iState);
					stateProb.push_back(all_states[nLoci][mi][jp1][jd]);
					statesToSample.push_back(tempStateVec);
					tempSum += all_states[nLoci][mi][jp1][jd];
					iState++;
				}
			}
			// normalize probabilities
			for(int s = 0, max = stateProb.size(); s < max; s++) stateProb[s] /= tempSum;
			
			// vector to use with sampleC below
			vector <vector <int> > possibleMove (5, vector <int> (3, 0)); // all possible moves between states
			// order is mi, miss p1, miss d
			vector <int> moveIndex (5);
			for(int i = 0; i < 5; i++) moveIndex[i] = i;
			possibleMove[1][0] = 1;
			int tempi = 2;
			for(int m1 = 0; m1 < 2; m1++){ // m1-3: 0 is not missing, 1 is missing
				for(int m3 = 0; m3 < 2; m3++){
					if(m1 + m3 == 0) continue; // already handled no missing above
					possibleMove[tempi][1] = m1;
					possibleMove[tempi][2] = m3;
					tempi++;
				}
			}

			for(int n = 0; n < itersPerMIC[mi]; n++){
				if(n % 100 == 0) Rcpp::checkUserInterrupt();
				
				// simulate pair - backwards algorithm
				vector <int> p1Genos (nLoci, 0);
				vector <int> dGenos (nLoci, 0);
				
				// first choose state (given number of mi)
				vector <int> curState (statesToSample[sampleC(stateIndex, stateProb, rg)]);
				// choose loci with MIs and missing genotypes
				for(int i = nLoci; i > 0; i--){
					vector <double> moveProb (5, 0); // probabilities of moves
					// no missing, no MI (probability you were in curState, and you stayed there)
					moveProb[0] = all_states[i-1][curState[0]][curState[1]][curState[2]] * 
						pow(1 - prob_missing_geno[i-1], 2) * (1 - pMI_noMiss[i-1]);
					// no missing, MI
					if(curState[0] > 0) moveProb[1] = all_states[i-1][curState[0] - 1][curState[1]][curState[2]] * 
						pow(1 - prob_missing_geno[i-1], 2) * pMI_noMiss[i-1];
					
					// all missing options
					tempi = 2;
					for(int m1 = 0; m1 < 2; m1++){ // m1-3: 0 is not missing, 1 is missing
						if(curState[1] - m1 < 0) continue;
						for(int m3 = 0; m3 < 2; m3++){
							if(curState[2] - m3 < 0) {
								tempi++;
								continue;
							}
							if(m1 + m3 == 0) continue; // already handled no missing above
							moveProb[tempi] = all_states[i-1][curState[0]][curState[1] - m1][curState[2] - m3] *
								pow(prob_missing_geno[i-1], m1 + m3) * 
								pow(1 - prob_missing_geno[i-1], 2 - (m1 + m3));
							tempi++;
						}
					}
					
					// select the move made at this locus
					vector <int> move (possibleMove[sampleC(moveIndex, moveProb, rg)]);
					
					// sample genotypes given move
					vector <int> chosenGenos (2);
					if(move[0] == 1){
						// if MI, sample from MI genos
						chosenGenos = pairGenosLookup[i-1][sampleC(combo[i-1], prob_MI[i-1], rg)];
					} else if (move[1] + move[2] == 0){
						// if all observed and no MI, sample from non MI genos
						chosenGenos = pairGenosLookup[i-1][sampleC(combo[i-1], prob_no_MI[i-1], rg)];
					} else {
						// if some are missing, sample from all genos, then add missing
						chosenGenos = pairGenosLookup[i-1][sampleC(combo[i-1], prob_all[i-1], rg)];
						if(move[1] == 1) chosenGenos[0] = -9;
						if(move[2] == 1) chosenGenos[1] = -9;
					}
					
					p1Genos[i-1] = chosenGenos[0];
					dGenos[i-1] = chosenGenos[1];

					// change current state for next locus
					for(int j = 0; j < 3; j++) curState[j] -= move[j];
					
				} // end simulating pair
				
				// no need to filter b/c only simulating pairs with valid numbers of observed MIs

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
								// pair
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
					if(llr >= llrToTestC[i]) falsePosPerStrata[mi][i]++;
				}
			} // end for n
		} // end for mi
		
		Rcpp::NumericVector col_pop;
		Rcpp::NumericVector col_llr;
		Rcpp::NumericVector col_mi;
		Rcpp::NumericVector col_p;
		Rcpp::NumericVector col_fp;
		Rcpp::NumericVector col_fpSD;
		for(int i = 0, ttmax = llrToTestC.size(); i < ttmax; i++){
			for(int j = 0, tmax = falsePosPerStrata.size(); j < tmax; j++){
				col_pop.push_back(pop);
				col_llr.push_back(llrToTestC[i]);
				col_mi.push_back(j);
				col_p.push_back(pNumMI_OBS[j]);
				if(itersPerMIC[j] == 0 || pNumMI_OBS[j] <= numeric_limits<double>::min()){
					col_fp.push_back(0);
					col_fpSD.push_back(0);
				} else {
					col_fp.push_back(falsePosPerStrata[j][i] / itersPerMIC[j]);
					double tVar = varianceEstim(itersPerMIC[j], falsePosPerStrata[j][i], 
                                 falsePosPerStrata[j][i] / itersPerMIC[j]);
					col_fpSD.push_back(sqrt(tVar));
					// calculate SD as the strata SD are calculated
					fp_SD_total[i] += pow(pNumMI_OBS[j], 2) * (tVar / itersPerMIC[j]);
				}
			}
		}
		fullOut.push_back(Rcpp::DataFrame::create(Rcpp::Named("Pop") = col_pop,
                                            Rcpp::Named("llrThreshold") = col_llr,
                                            Rcpp::Named("count_MI") = col_mi,
                                            Rcpp::Named("prob_count_MI") = col_p,
                                            Rcpp::Named(fpColName) = col_fp,
                                            Rcpp::Named(fpSDColName) = col_fpSD
		));
		for(int i = 0, max2 = llrToTestC.size(); i < max2; i++){
			// calculate the mean accross strata
			double strata_fpMean = 0;
			for(int j = 0, tmax = falsePosPerStrata.size(); j < tmax; j++) if(itersPerMIC[j] > 
                    0) strata_fpMean += pNumMI_OBS[j] * (falsePosPerStrata[j][i] / itersPerMIC[j]);
			popResults.push_back(pop);
			llrRecord.push_back(llrToTestC[i]);
			falsePosAll.push_back(strata_fpMean);
			SDfalsePosAll.push_back(sqrt(fp_SD_total[i]));
		}
	} // end for pop
	
	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Pop") = popResults,
                                                   Rcpp::Named("llrThreshold") = llrRecord,
                                                   Rcpp::Named(fpColName) = falsePosAll,
                                                   Rcpp::Named(fpSDColName) = SDfalsePosAll
	);
	fullOut.push_front(results);
	return fullOut;
}


// this older function does not use Markov for missing genos and starta are obs number of MI
//   BEFORE missing genotypes are "chosen" - not very intuitive
// 
//' estimating positive error rates for parent - offspring pairs
//' stratified Monte Carlo used for estimating false positive
//' @param baselineParams Dirichlet parameters for allele frequencies
//' @param unsampledPopParams Dirichlet parameters for allele frequencies
//' @param missingParams Beta parameters for missing genotypes (failure to genotype rate)
//' @param genotypeKey list of matrix for each locus, col1 is genotype, col2 is allele 1, col3 is allele 2
//' @param genotypeErrorRates list of matrix for each locus, rows are actual genotype, columns are observed,
//'   values are probability
//' @param llrToTest Vector of llr's to test as threshold values
//' @param N number of samples to take
//' @param trueRel an integer indicating the true relationship to simulate: 0 = unrelated, 1 = aunt,
//'   2 = half aunt, 3 = cousin of parent
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::List old_strat_ERRORsP(Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::List missingParams,
                           Rcpp::List genotypeKey,
                           Rcpp::List genotypeErrorRates, Rcpp::NumericVector llrToTest,
                           Rcpp::NumericVector itersPerMI,
                           int seed, int trueRel, double MIexcludeProb){
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

	vector <double> itersPerMIC;
	for(int i = 0, max = itersPerMI.length(); i < max; i++) itersPerMIC.push_back(itersPerMI[i]);
	
	vector <double> llrToTestC;
	for(int i = 0, max = llrToTest.length(); i < max; i++) llrToTestC.push_back(llrToTest[i]);
	
	// end conversion of types
	////////////////////

	if(itersPerMIC.size() != (genotypeKeyC.size() + 1)) Rcpp::stop("Error: the length of itersPerMI is not equal to the number of loci + 1.");

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
	

	vector <int> popResults;
	vector <double> llrRecord;
	vector <double> falsePosAll;
	vector <double> SDfalsePosAll;

	// for output of full results
	Rcpp::List fullOut;
	
	// choose values depending of true relationship being simulated
	vector <double> k_prob(3,0);
	string fpColName ("falsePos");
	string fpSDColName ("falsePos");
	if(trueRel == 0){
		// Unrelated
		k_prob[0] = 1;
		k_prob[1] = 0;
		// k_prob[2] = 0;
		fpColName += "Unrel";
		fpSDColName += "Unrel";
	} else if (trueRel == 1){
		// Aunt
		k_prob[0] = 0.5;
		k_prob[1] = 0.5;
		// k_prob[2] = 0;
		fpColName += "Aunt";
		fpSDColName += "Aunt";
	} else if (trueRel == 2){
		// half aunt
		k_prob[0] = 0.75;
		k_prob[1] = 0.25;
		// k_prob[2] = 0;
		fpColName += "HalfAunt";
		fpSDColName += "HalfAunt";
	} else if (trueRel == 3){
		// cousin of parent
		k_prob[0] = 0.875;
		k_prob[1] = 0.125;
		// k_prob[2] = 0;
		fpColName += "ParCous";
		fpSDColName += "ParCous";
	} else {
		Rcpp::stop("Internal error: trueRel not recognized");
	}
	fpSDColName += "SD";
	
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
		 * calculate genotype likelihoods for relationship (example here is aunt/uncle - niece/nephew)
		 * Index 1 is locus, 2 is aunt/uncle geno, 3 is niece/nephew geno
		 * value is likelihood (NOT log)
		 */
		vector <vector <vector <double> > > lGenos_tRel;
		pairwiseK_calc2(genotypeKeyC, lGenos_base, baselineParamsC, unsampledPopParamsC, pop,
	                      k_prob, lGenos_tRel);

		/*	
		 *	calculate genotype LOG-likelihoods for (example aunt/uncle - niece/nephew) relationship
		 *	and for all three being unrelated
			given OBSERVED genotypes - but only if NO MISSING genotypes in the pair
		 	Index 1 is locus, 2 is aunt/uncle geno, 3 is niece/nephew geno
			Value is LOG-likelihood
		 */
		vector <vector <vector <double> > > lGenos_tRel_OBS;
		create_pairwise_OBSvector(lGenos_tRel, genotypeErrorRatesC, lGenos_tRel_OBS);

		
		/*	Note that this could be done with pairwise K functions, but is different for historical reasons
		 * 	ie I wrote it before writing the pairwise K functions!
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
		 *	and for being unrelated
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
		
		// these are all p(MI | neither genotype is missing)

		// first calculate the probability of MI at each locus | relationship
		vector <double> pMI (genotypeKeyC.size(), 0.0); // p(MI) observed for each locus
		calcProbMIperLocus_Pair(genotypeKeyC, genotypeErrorRatesC, lGenos_tRel, pMI);
		
		// then calculate probability of sum(MI) = x for x in [0, num_loci]
		vector <double> pTotalMI (genotypeKeyC.size() + 1, 0); // probability of having x observed MI
		vector <vector <double> > all_pTotalMI; // pTotalMI for each step
		calcProbSumMI_returnAll(pMI, pTotalMI, all_pTotalMI);
		
		int miLimit = genotypeKeyC.size();
		if(MIexcludeProb > 0){
			// first calculate the probability of MI at each locus | true parent
			vector <double> pMI_sP (genotypeKeyC.size(), 0.0); // p(MI|true parent) for each locus
			calcProbMIperLocus_sP(genotypeKeyC, genotypeErrorRatesC, lGenos_sP, pMI_sP);
			
			// now calculate the probability of sum(MI) = x
			vector <double> pTotalMI_sP (genotypeKeyC.size() + 1, -9); // probability of having x observed MI | true parent
			calcProbSumMI(pMI_sP, pTotalMI_sP);
			
			// now find miLimit such that P(sum(MI) > miLimit) < MIexcludeProb
			double runningSum = 0.0;
			for(int i = 0, max = pTotalMI_sP.size(); i < max; i++){
				runningSum += pTotalMI_sP[i];
				if(runningSum > (1 - MIexcludeProb)){
					miLimit = i;
					break;
				}
			}
			Rcpp::Rcout<<"The maximum number of Mendalian incompatibilities allowed"<<
				" is: "<<miLimit<<". The probability of exclusion for a true parent (given no missing genotypes) is estimated as: "<<
					1 - runningSum<<".\n";
		}
		
		
		/*
		 * calculate probability of pairs of OBSERVED genotypes given MI at a locus or given no MI at a locus 
		 * indices are: locus, baseline genotype, desc genotype
		 * value is probability NOT log
		 */
		vector <vector <vector <double> > > lGenos_OBS_sample_MI;
		vector <vector <vector <double> > > lGenos_OBS_sample_no_MI;
		strat_lGenosBuilder_Pairwise(genotypeKeyC, genotypeErrorRatesC, lGenos_tRel_OBS,
                        lGenos_OBS_sample_MI,
								lGenos_OBS_sample_no_MI);
		
		/*
		 * now making lGenos_OBS_sample_(MI / no_MI) easier to sample from
		 * So for each locus, combo has the index of each trio Genotype combination
		 * probability has the probability of each combo
		 * pairGenosLookup has a vector giving the two genotypes in each combo (p1, d)
		 */
		vector <vector <int> > combo; // index of the genotype combos
		vector <vector <double> > prob_MI; // the probabilities of each genotype combo
		vector <vector <double> > prob_no_MI; // the probabilities of each genotype combo
		vector <vector <vector <int> > > pairGenosLookup; // the genotype combos
		for(int i = 0; i < nLoci; i++){
			vector <int> tempCombo;
			vector <double> tempProb_MI;
			vector <double> tempProb_no_MI;
			vector <vector <int> > tempGenos;
			int pairIndex = 0;
			for(int p1 = 0, max = lGenos_OBS_sample_MI[i].size(); p1 < max; p1++){
				for(int d = 0; d < max; d++){
					tempProb_MI.push_back(lGenos_OBS_sample_MI[i][p1][d]);
					tempProb_no_MI.push_back(lGenos_OBS_sample_no_MI[i][p1][d]);
					tempCombo.push_back(pairIndex);
					vector <int> tempVec (2,-9);
					tempVec[0] = p1;
					tempVec[1] = d;
					tempGenos.push_back(tempVec);
					pairIndex++;
				}
			}
			combo.push_back(tempCombo);
			prob_MI.push_back(tempProb_MI);
			prob_no_MI.push_back(tempProb_no_MI);
			pairGenosLookup.push_back(tempGenos);
		}
		
		// some results storage for full results reporting
		vector <vector <double> > falsePosPerStrata; // indices: num of MI, llr threshold
		for(int i=0, tmax = itersPerMIC.size(); i < tmax; i++){
			vector <double> tVec (llrToTestC.size(), 0);
			falsePosPerStrata.push_back(tVec);
		}

		// inititate results storage for this pop - one entry for each llr
		vector <double> fp_SD_total (llrToTestC.size(), 0);
		
		// two vectors to use with sampleC below
		vector <double> probs_MI_NoMI (2, 0);
		vector <int> select_MI (2);
		select_MI[0] = 1;
		select_MI[1] = 0;

		for(int mi = 0, maxMI = itersPerMIC.size(); mi < maxMI; mi++){
			Rcpp::checkUserInterrupt();
			// check if prob of seeing this many MIs is zero
			// if it's zero, throw a warning and skip this stratum
			if (pTotalMI[mi] <= numeric_limits<double>::min() && itersPerMIC[mi] > 0){
				string errMess = "The probability of observing ";
				errMess += to_string(mi); 
				errMess += " Mendelian incompatibilities for this relationship is 0. "; 
				errMess += "Skipping all iterations in this stratum."; 
				Rcpp::warning(errMess);
				continue;
			}
			
			for(int n = 0; n < itersPerMIC[mi]; n++){
				if(n % 100 == 0) Rcpp::checkUserInterrupt();
				
				// simulate pair
				vector <int> p1Genos (nLoci, 0);
				vector <int> dGenos (nLoci, 0);
				
				// backwards algorithm
				// first choose loci with MI
				// all_pTotalMI[i][j] this is, after locus i, prob of being in stage j 
				vector <bool> lociWithMI (nLoci, false);
				int state = mi; // the current state
				for(int i = nLoci; i > 0; i--){
					if(state == 0) break;
					// so, need probabilities of MI at this locus
					// the probability you were in (state - 1), then observed an MI
					probs_MI_NoMI[0] = 0;
					if(state > 0) probs_MI_NoMI[0] = all_pTotalMI[i-1][state - 1] * pMI[i-1];
					// probability you were previously in state, then did not observe an MI
					probs_MI_NoMI[1] = all_pTotalMI[i-1][state] * (1 - pMI[i-1]);
					
					if(sampleC(select_MI, probs_MI_NoMI, rg) == 1){
						lociWithMI[i-1] = true;
						state--;
					}
				}

				// now choose genotypes
				for(int i = 0; i < nLoci; i++){
					// determine if MI
					int tPair;
					if(lociWithMI[i]){ // with MI
						tPair = sampleC(combo[i], prob_MI[i], rg);
					} else { // without MI
						tPair = sampleC(combo[i], prob_no_MI[i], rg);
					}
					p1Genos[i] = pairGenosLookup[i][tPair][0];
					dGenos[i] = pairGenosLookup[i][tPair][1];
					// add missing genotypes as appropriate
					if(randMissing[i](rg)) p1Genos[i] = -9;
					if(randMissing[i](rg)) dGenos[i] = -9;
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
					if(skip) continue;
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
									// pair
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
					if(llr >= llrToTestC[i]) falsePosPerStrata[mi][i]++;
				}
			} // end for n
		} // end for mi
		
		Rcpp::NumericVector col_pop;
		Rcpp::NumericVector col_llr;
		Rcpp::NumericVector col_mi;
		Rcpp::NumericVector col_p;
		Rcpp::NumericVector col_fp;
		Rcpp::NumericVector col_fpSD;
		for(int i = 0, ttmax = llrToTestC.size(); i < ttmax; i++){
			for(int j = 0, tmax = falsePosPerStrata.size(); j < tmax; j++){
				col_pop.push_back(pop);
				col_llr.push_back(llrToTestC[i]);
				col_mi.push_back(j);
				col_p.push_back(pTotalMI[j]);
				if(itersPerMIC[j] == 0 || pTotalMI[j] <= numeric_limits<double>::min()){
					col_fp.push_back(0);
					col_fpSD.push_back(0);
				} else {
					col_fp.push_back(falsePosPerStrata[j][i] / itersPerMIC[j]);
					double tVar = varianceEstim(itersPerMIC[j], falsePosPerStrata[j][i], 
                                           falsePosPerStrata[j][i] / itersPerMIC[j]);
					col_fpSD.push_back(sqrt(tVar));
					// calculate SD as the strata SD are calculated
					fp_SD_total[i] += pow(pTotalMI[j], 2) * (tVar / itersPerMIC[j]);
				}
			}
		}
		fullOut.push_back(Rcpp::DataFrame::create(Rcpp::Named("Pop") = col_pop,
                                                Rcpp::Named("llrThreshold") = col_llr,
                                                Rcpp::Named("count_MI") = col_mi,
                                                Rcpp::Named("prob_count_MI") = col_p,
                                                Rcpp::Named(fpColName) = col_fp,
                                                Rcpp::Named(fpSDColName) = col_fpSD
																));
		for(int i = 0, max2 = llrToTestC.size(); i < max2; i++){
			// calculate the mean accross strata
			double strata_fpMean = 0;
			for(int j = 0, tmax = falsePosPerStrata.size(); j < tmax; j++) if(itersPerMIC[j] > 
						0) strata_fpMean += pTotalMI[j] * (falsePosPerStrata[j][i] / itersPerMIC[j]);
			popResults.push_back(pop);
			llrRecord.push_back(llrToTestC[i]);
			falsePosAll.push_back(strata_fpMean);
			SDfalsePosAll.push_back(sqrt(fp_SD_total[i]));
		}
	} // end for pop
	
	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Pop") = popResults,
                                                   Rcpp::Named("llrThreshold") = llrRecord,
                                                   Rcpp::Named(fpColName) = falsePosAll,
                                                   Rcpp::Named(fpSDColName) = SDfalsePosAll
																	);
	fullOut.push_front(results);
	return fullOut;
}
