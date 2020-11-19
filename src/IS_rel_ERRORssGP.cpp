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

//' estimatign error rates for single-sided grandparent pair for other true relationship types
//' Importance sampling Monte Carlo used for estimating false positive
//' @param baselineParams Dirichlet parameters for allele frequencies
//' @param unsampledPopParams Dirichlet parameters for allele frequencies
//' @param missingParams Beta parameters for missing genotypes (failure to genotype rate)
//' @param genotypeKey list of matrix for each locus, col1 is genotype, col2 is allele 1, col3 is allele 2
//' @param genotypeErrorRates list of matrix for each locus, rows are actual genotype, columns are observed,
//'   values are probability
//' @param llrToTest Vector of llr's to test as threshold values
//' @param N number of samples to take
//' @param trueRel integer corresponding to true relationship type
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame IS_rel_ERRORssGP(Rcpp::List baselineParams,
                          Rcpp::List unsampledPopParams, Rcpp::List missingParams,
                          Rcpp::List genotypeKey,
                          Rcpp::List genotypeErrorRates, Rcpp::NumericVector llrToTest,
                          int N, int trueRel,
                          int seed, double MIexcludeProb,
                          int maxMissingGenos
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
	vector <double> falsePosAll;
	vector <double> SDfalsePosAll;
	vector <double> falseNegAll;
	vector <double> SDfalseNegAll;
	
	// choose values depending on true relationships being simulated
	vector <double> k_prob_1(3,0); // relationship of poten gp1 TO TRUE gp1
	vector <double> k_prob_2(3,0); // relationship of poten gp2 TO TRUE gp2
	string fpColName ("falsePos");
	string fpSDColName;
	if(trueRel == 0){
		// Unrelated and Unrelated
		k_prob_1[0] = 1;
		k_prob_2[0] = 1;
		fpColName += "Unrel";
	} else if (trueRel == 1) {
		// True GP and great-aunt
		k_prob_1[2] = 1;
		k_prob_2[0] = 0.25; // full sibs
		k_prob_2[1] = 0.5;
		k_prob_2[2] = 0.25;
		fpColName += "True_GAunt";
	} else if (trueRel == 2) {
		k_prob_1[2] = 1;
		k_prob_2[0] = 1;
		fpColName += "True_Unrel";
	} else if (trueRel == 3) {
		k_prob_1[2] = 1;
		k_prob_2[0] = .5; // half-sibs
		k_prob_2[1] = .5;
		fpColName += "True_HGAunt";
	} else if (trueRel == 4) {
		k_prob_1[2] = 1;
		k_prob_2[0] = .75; // cousin
		k_prob_2[1] = .25;
		fpColName += "True_GpCous";
	} else if (trueRel == 5) {
		k_prob_1[0] = 1;
		k_prob_2[0] = 0.25; // full sibs
		k_prob_2[1] = 0.5;
		k_prob_2[2] = 0.25;
		fpColName += "GAunt_Unrel";
	} else if (trueRel == 6) {
		k_prob_1[0] = 1;
		k_prob_2[0] = .5; // half-sibs
		k_prob_2[1] = .5;
		fpColName += "HGAunt_Unrel";
	} else if (trueRel == 7) {
		k_prob_1[0] = 1;
		k_prob_2[0] = .75; // cousin
		k_prob_2[1] = .25;
		fpColName += "GpCous_Unrel";
	} else if (trueRel == 8) {
		k_prob_1[0] = 0.25; // full sibs
		k_prob_1[1] = 0.5;
		k_prob_1[2] = 0.25;
		k_prob_2[0] = 0.25; // full sibs
		k_prob_2[1] = 0.5;
		k_prob_2[2] = 0.25;
		fpColName += "GAunt";
	} else if (trueRel == 9) {
		k_prob_1[0] = 0.25; // full sibs
		k_prob_1[1] = 0.5;
		k_prob_1[2] = 0.25;
		k_prob_2[0] = .5; // half-sibs
		k_prob_2[1] = .5;
		fpColName += "GAunt_HGAunt";
	} else if (trueRel == 10) {
		k_prob_1[0] = 0.25; // full sibs
		k_prob_1[1] = 0.5;
		k_prob_1[2] = 0.25;
		k_prob_2[0] = .75; // cousin
		k_prob_2[1] = .25;
		fpColName += "Gaunt_GpCous";
	} else if (trueRel == 11) {
		k_prob_1[0] = .5; // half-sibs
		k_prob_1[1] = .5;
		k_prob_2[0] = .5; // half-sibs
		k_prob_2[1] = .5;
		fpColName += "HGAunt";
	} else if (trueRel == 12) {
		k_prob_1[0] = .5; // half-sibs
		k_prob_1[1] = .5;
		k_prob_2[0] = .75; // cousin
		k_prob_2[1] = .25;
		fpColName += "HGAunt_GpCous";
	} else if (trueRel == 13) {
		k_prob_1[0] = .75; // cousin
		k_prob_1[1] = .25;
		k_prob_2[0] = .75; // cousin
		k_prob_2[1] = .25;
		fpColName += "GpCous";
	} else {
		Rcpp::stop("Internal error: trueRel not recognized");
	}
	fpSDColName = fpColName + "SD";
	
	// testing
	if(k_prob_1[0] + k_prob_1[1] + k_prob_1[2] != 1) Rcpp::stop("Internal error k_prob_1.");
	if(k_prob_2[0] + k_prob_2[1] + k_prob_2[2] != 1) Rcpp::stop("Internal error k_prob_2.");
	

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
		
		// determine max number of MIs
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
			Rcpp::Rcout<<"The maximum number of Mendalian incompatibilities allowed"<<
				" is: "<<miLimit<<". The probability of exclusion for a true grandparent pair (given no missing genotypes) is estimated as: "<<
					1 - runningSum<<".\n";
		}
		
		/*
		 * calculate log likelihood of Observed and true trio genotypes given specified true relationships
		 * above
		 * These are rearranged and used to sample trio genotypes given Mi or no MI
		 * Indices are: locus, poten gp1 genotype, poten gp2 geno, grandchild geno
		 * Value of lGenos_trio_genos_OBS is LOG likelihood
		 * Value of lGenos_trio_genos is likelihood (NOT log)
		 */
		vector <vector <vector <vector <double> > > > lGenos_trio_genos;
		vector <vector <vector <vector <double> > > > lGenos_trio_genos_OBS;
		// scoped to save memory but no need to make separate function
		{
			vector <vector <vector <double> > > lGenos_tRel_gp1;
			pairwiseK_calc2(genotypeKeyC, lGenos_base, baselineParamsC, baselineParamsC, pop,
                   k_prob_1, lGenos_tRel_gp1);
			
			vector <vector <vector <double> > > lGenos_tRel_gp2;
			pairwiseK_calc2(genotypeKeyC, lGenos_base, baselineParamsC, baselineParamsC, pop,
                   k_prob_2, lGenos_tRel_gp2);
		
			for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ // for each locus
				// then calc probabilities of all potential trio genotypes
				// initialize with 0's
				vector <vector <vector <double> > > tempLoc_true_probs (genotypeKeyC[i].size(), 
                                                            vector <vector <double> > (genotypeKeyC[i].size(), 
                                                                                       vector <double> (genotypeKeyC[i].size(), 0)));
				lGenos_trio_genos.push_back(vector <vector <vector <double> > > (genotypeKeyC[i].size(),
                                                                         vector <vector <double> > (genotypeKeyC[i].size(),
                                                                                                    vector <double> (genotypeKeyC[i].size(), 0))));
				
				for(int gp1 = 0, mg = genotypeKeyC[i].size(); gp1 < mg; gp1++){
					for(int gp2 = 0; gp2 < mg; gp2++){ // b/c true relationships may be different, have to calc all combos
						for(int d = 0; d < mg; d++){
							// P(d)*P(gp1|d)*P(gp2|d)
							// for all possible true grandparent genotypes
							for(int tgp1 = 0; tgp1 < mg; tgp1++){
								for(int tgp2 = 0; tgp2 < mg; tgp2++){
									// prob of d and true gp genos given poten gp genos 
									tempLoc_true_probs[gp1][gp2][d] += exp(log(lGenos_ssGP[i][tgp1][tgp2][d]) + 
										log(lGenos_tRel_gp1[i][gp1][tgp1]) +
										log(lGenos_tRel_gp2[i][gp2][tgp2]) - lGenos_base[i][tgp1] - lGenos_base[i][tgp2]);
								}
							}
							lGenos_trio_genos[i][gp1][gp2][d] = tempLoc_true_probs[gp1][gp2][d];
						}
					}
				}
				
				// then for each OBS genotype, calc probability given probs of true genotypes and err rates
				// initialize with 0's
				lGenos_trio_genos_OBS.push_back(vector <vector <vector <double> > > (genotypeKeyC[i].size(),
                                                                       vector <vector <double> > (genotypeKeyC[i].size(),
                                                                                                  vector <double> (genotypeKeyC[i].size(), 0))));
				for(int gp1 = 0, mg = genotypeKeyC[i].size(); gp1 < mg; gp1++){
					for(int gp2 = 0; gp2 < mg; gp2++){ // b/c true relationships may be different, have to calc all combos
						for(int d = 0; d < mg; d++){
							// for each possible true genotype
							for(int tgp1 = 0; tgp1 < mg; tgp1++){
								for(int tgp2 = 0; tgp2 < mg; tgp2++){
									for(int td = 0; td < mg; td++){
										lGenos_trio_genos_OBS[i][gp1][gp2][d] += tempLoc_true_probs[tgp1][tgp2][td] * 
											genotypeErrorRatesC[i][tgp1][gp1] * genotypeErrorRatesC[i][tgp2][gp2] * 
											genotypeErrorRatesC[i][td][d];
									}
								}
							} // end for tgp1
							lGenos_trio_genos_OBS[i][gp1][gp2][d] = log(lGenos_trio_genos_OBS[i][gp1][gp2][d]); // change to log likelihood
						}
					}
				}
			} // end for each locus
		} // end scoping
		
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
		vector <double> falsePos (llrToTestC.size(), 0);
		vector <double> falseNeg (llrToTestC.size(), 0);
		vector <double> SSfalsePos (llrToTestC.size(), 0); // sum of squares for variance calcs
		// don't need sum squares for falseNeg right now b/c all are either 0 or 1 (no importance sampling)
		
		for(int n = 0; n < N; n++){
			if(n % 100 == 0) Rcpp::checkUserInterrupt();
			
			// simulate trio
			vector <int> gmaGenos (nLoci, 0);
			vector <int> gpaGenos (nLoci, 0);
			vector <int> dGenos (nLoci, 0);
			vector <int> sampTrio (3, -1);
			for(int i = 0; i < nLoci; i++){ // for each locus
				// sample observed genotypes
				sampTrio = trioGenosLookup[i][sampleC(combo[i], prob[i], rg)];
				gpaGenos[i] = sampTrio[0];
				gmaGenos[i] = sampTrio[1];
				dGenos[i] = sampTrio[2];
			}
			
			// backwards algorithm to add missing genotypes
			// choose state : number of missing genotypes
			int state_gpa = sampleC(numMissingState, numMissingProb, rg);
			int state_gma = sampleC(numMissingState, numMissingProb, rg);
			int state_d = sampleC(numMissingState, numMissingProb, rg);
			// choose missing genotypes
			for(int i = nLoci; i > 0; i--){
				if(state_gpa + state_gma + state_d == 0) break;
				//gpa
				prob_miss_nomiss[0] = 0;
				if(state_gpa > 0) prob_miss_nomiss[0] = probTotalMiss_all[i-1][state_gpa - 1] * prob_missing_geno[i-1];
				prob_miss_nomiss[1] = probTotalMiss_all[i-1][state_gpa] * (1 - prob_missing_geno[i-1]);
				if(sampleC(select_miss, prob_miss_nomiss, rg) == 1){
					gpaGenos[i-1] = -9;
					state_gpa--;
				}
				//gma
				prob_miss_nomiss[0] = 0;
				if(state_gma > 0) prob_miss_nomiss[0] = probTotalMiss_all[i-1][state_gma - 1] * prob_missing_geno[i-1];
				prob_miss_nomiss[1] = probTotalMiss_all[i-1][state_gma] * (1 - prob_missing_geno[i-1]);
				if(sampleC(select_miss, prob_miss_nomiss, rg) == 1){
					gmaGenos[i-1] = -9;
					state_gma--;
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
			
			// determine if have allowable number of MIs
			if(MIexcludeProb > 0){
				int countMI = 0;
				bool skip = false;
				for(int j = 0; j < nLoci; j++){ // for each locus
					int gp1Geno = gpaGenos[j];
					int gp2Geno = gmaGenos[j];
					int dGeno = dGenos[j];
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
			}
			
			// calc LLR
			double uLLH = 0.0;
			double gpLLH = 0.0;
			double corrLLH = 0.0;
			for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ //for each locus
				// observed genotypes
				int obs_gp1 = gpaGenos[j];
				int obs_gp2 = gmaGenos[j];
				int obs_d = dGenos[j];
				
				// if desc missing, no information about the relationship
				if(obs_d == -9) continue;
				
				// if missing genotypes, information if only one missing OR if an unsampled pop
				//   b/c allele freqs will be different whether its a descendant or not
				if(obs_gp1 == -9 || obs_gp2 == -9){
					// could make lookup tables for these
					double u_likelihood = 0.0; // likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
					double gp_likelihood = 0.0;
					double corr_likelihood = 0.0;
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
								// correction
								corr_likelihood += lGenos_trio_genos[j][gp1][gp2][d] *
									pOA_gp1 * pOA_gp2 * pOA_d;
							}
						}
					}
					uLLH += log(u_likelihood);
					gpLLH += log(gp_likelihood);
					corrLLH += log(corr_likelihood);
				} else { // no missing genotypes
					uLLH += lGenos_Unrelated_OBS[j][obs_gp1][obs_gp2][obs_d];
					gpLLH += lGenos_ssGP_OBS[j][obs_gp1][obs_gp2][obs_d];
					corrLLH += lGenos_trio_genos_OBS[j][obs_gp1][obs_gp2][obs_d];
				}
			}
			double llr = gpLLH - uLLH;
			double corr_ratio = exp(corrLLH - gpLLH);
			// determine number under and over threshold
			for(int i = 0, max = llrToTestC.size(); i < max; i++){
				if(llr < llrToTestC[i]){
					falseNeg[i]++;
				} else {
					falsePos[i] += corr_ratio; // correction for sampling distribution
					SSfalsePos[i] += pow(corr_ratio, 2);
				}
			}
		} // end for N
		for(int i = 0, max = llrToTestC.size(); i < max; i++){
			double fpMean = falsePos[i] / N; // need for variance calc
			double fnMean = falseNeg[i] / N;
			
			popResults.push_back(pop);
			llrRecord.push_back(llrToTestC[i]);
			falsePosAll.push_back(fpMean);
			falseNegAll.push_back(fnMean);
			SDfalsePosAll.push_back(sqrt(varianceEstim(N, SSfalsePos[i], fpMean)));
			SDfalseNegAll.push_back(sqrt(varianceEstim(N, falseNeg[i], fnMean))); // since all are 1 SS is just sum
		}
		
	} // end for pop
	
	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Pop") = popResults,
                                                   Rcpp::Named("llrThreshold") = llrRecord,
                                                   Rcpp::Named(fpColName) = falsePosAll,
                                                   Rcpp::Named(fpSDColName) = SDfalsePosAll,
                                                   Rcpp::Named("falseNeg") = falseNegAll,
                                                   Rcpp::Named("falseNegSD") = SDfalseNegAll
	);
	
	return results;
}
