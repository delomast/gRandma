#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <string>
#include <random>
#include <limits>
#include "misc_math.h"
#include "utils.h"
#include "mend_incompat.h"

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]




//' estimating false positive error rates for single-sided grandparent pair vs unrelated
//' with individuls from one population being assigned to a different population
//' stratified Monte Carlo used for estimating false positive
//' strata are number of observed mendelian incompatibilities in the trio
//' 
//' 
//' 
//' @param baselineParams Dirichlet parameters for allele frequencies
//' @param unsampledPopParams Dirichlet parameters for allele frequencies
//' @param missingParams Beta parameters for missing genotypes (failure to genotype rate)
//' @param genotypeKey list of matrix for each locus, col1 is genotype, col2 is allele 1, col3 is allele 2
//' @param genotypeErrorRates list of matrix for each locus, rows are actual genotype, columns are observed,
//'   values are probability
//' @param llrToTest Vector of llr's to test as threshold values
//' @param itersPerMI number of samples to take in each strata
//' @param seed random number generation seed
//' @param skipBaseline added unsampled pops to skip as baseline
//' @param MIexcludeProb maximum prob to lead to GP exclusion
//' @param maxMissingGenos maximum number of missing genotypes allowed in an individual
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::List strat_otherPopERRORssGP(Rcpp::List baselineParams,
                                     Rcpp::List unsampledPopParams, Rcpp::List missingParams,
                                     Rcpp::List genotypeKey,
                                     Rcpp::List genotypeErrorRates, Rcpp::NumericVector llrToTest,
                                     Rcpp::NumericVector itersPerMI,
                                     int seed, Rcpp::NumericVector skipBaseline, double MIexcludeProb,
                                     int maxMissingGenos){
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
	
	vector <double> itersPerMIC;
	for(int i = 0, max = itersPerMI.length(); i < max; i++) itersPerMIC.push_back(itersPerMI[i]);
	
	vector <int> skipBaselineC;
	for(int i = 0, max = skipBaseline.length(); i < max; i++) skipBaselineC.push_back(skipBaseline[i]);
	
	// end conversion of types
	////////////////////
	
	int nLoci = genotypeKeyC.size();
	if(maxMissingGenos > nLoci) Rcpp::stop("maxMissingGenos cannot be greater than the number of loci.");
	
	if(MIexcludeProb == 0) Rcpp::Rcout << "MIexcludeProb was not greater than zero, " <<
		"so only the LLR will be used to accept or reject assignments.\n";
	
	// initiate random number generator
	mt19937 rg (seed);
	
	vector <int> popResults; // pop that trios are simulated from
	vector <int> popAssign; // pop that offspring are assigned to
	vector <double> llrRecord;
	vector <double> falsePosAll;
	vector <double> SDfalsePosAll;
	
	// for output of full results
	Rcpp::List fullOut;
	
	// for each baseline population
	for(int pop = 0, max = baselineParamsC.size(); pop < max; pop++){
		Rcpp::Rcout<<"Now considering descendants from baseline population "<<pop + 1<<"\n";

		/*
		* calculate genotype likelihoods for a descendant of this population sampled at random
		* Index 1 is locus
		* Index 2 is genotype
		* Value is LOG-likelihood
		*/
		vector <vector <double> > lGenos_randomDescendant;
		for(int i = 0; i < nLoci; i++){
			vector <double> tempLH;
			for(int j = 0, max2 = genotypeKeyC[i].size(); j < max2; j++){
				vector <double> k1 (baselineParamsC[pop][i].size(), 0);
				vector <double> k2 (baselineParamsC[pop][i].size(), 0);
				k1[genotypeKeyC[i][j][0]]++;
				k2[genotypeKeyC[i][j][1]]++;
				if(genotypeKeyC[i][j][0] == genotypeKeyC[i][j][1]){ // alleles are the same
					// prob sample allele from both
					tempLH.push_back(logMultPMF(k1, baselineParamsC[pop][i]) + 
						logMultPMF(k2, unsampledPopParamsC[pop][i])
					);
				} else {
					// prob sample allele1 from one, allele2 from other + prob of vice versa
					tempLH.push_back(log(
							exp(logMultPMF(k1, baselineParamsC[pop][i]) + 
								logMultPMF(k2, unsampledPopParamsC[pop][i])) + 
								exp(logMultPMF(k2, baselineParamsC[pop][i]) + 
								logMultPMF(k1, unsampledPopParamsC[pop][i]))
					));
				}
			}
			lGenos_randomDescendant.push_back(tempLH);
		}

		for(int pop2 = 0, maxp2 = baselineParamsC.size(); pop2 < maxp2; pop2++){
			if(pop2 == pop) continue;
			bool skipBaseBool = false;
			for(int i = 0, maxSkip = skipBaselineC.size(); i < maxSkip; i++){
				if(pop2 == skipBaselineC[i]){
					skipBaseBool = true;
					break;
				}
			}
			if (skipBaseBool) continue;
			
			Rcpp::Rcout<<"Calculating error with baseline population "<<pop2 + 1<<"\n";
			// calculate genotype log-likelihoods from random sampled genotype
			// first index is locus, second is genotype, value is LOG-likelihood
			vector <vector <double> > lGenos_base;
			vector <vector <double> > lGenos_unsamp;
			
			for(int i = 0, max3 = genotypeKeyC.size(); i < max3; i++){ // for each locus
				// baselineParamsC[pop2][i] are the alpha values for the baseline pop
				// unsampledPopParamsC[pop2][i] are the alpha values for the corresponding unsampled pop
				// genotypeKeyC[i] is the key for the locus i
				vector <double> tempBase;
				vector <double> tempUnsamp;
				for(int j = 0, max2 = genotypeKeyC[i].size(); j < max2; j++){ // for each genotype
					vector <double> k (baselineParamsC[pop2][i].size(), 0);
					k[genotypeKeyC[i][j][0]]++;
					k[genotypeKeyC[i][j][1]]++;
					tempBase.push_back(logMultPMF(k, baselineParamsC[pop2][i]));
					tempUnsamp.push_back(logMultPMF(k, unsampledPopParamsC[pop2][i]));
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
			
			int miLimit = genotypeKeyC.size();
			if(MIexcludeProb > 0){
				// first calculate the probability of MI at each locus | true grandparents
				vector <double> pMI_ssGP (genotypeKeyC.size(), 0.0); // p(MI|true grandparents) for each locus
				calcProbMIperLocus(genotypeKeyC, genotypeErrorRatesC, lGenos_ssGP, pMI_ssGP);
				
				// now calculate the probability of sum(MI) = x | true grandparents
				vector <double> pTotalMI_ssGP (genotypeKeyC.size() + 1, -9); // probability of having x observed MI | true grandparents
				calcProbSumMI(pMI_ssGP, pTotalMI_ssGP);
				
				// now find miLimit such that P(sum(MI) > miLimit) < MIexcludeProb
				double runningSum = 0.0;
				for(int i = 0, max = pTotalMI_ssGP.size(); i < max; i++){
					runningSum += pTotalMI_ssGP[i];
					if(runningSum > (1 - MIexcludeProb)){
						miLimit = i;
						break;
					}
				}
				Rcpp::Rcout<<"The maximum number of Mendalian incompatibilities allowed"<<
					" is: "<<miLimit<<". The probability of exclusion for a true grandparent pair (given no missing genotypes) is estimated as: "<<
						1 - runningSum<<".\n";
			}
			
			if(static_cast<int>(itersPerMIC.size()) < (miLimit + 1)) Rcpp::stop("Error: the length of itersPerMI is not long enough for the number of allowed MIs in this population.");
			if((miLimit + 1) * pow(maxMissingGenos + 1, 3) * nLoci > 345000000) Rcpp::Rcout << "\nUsing a large amount of memory...\n\n";

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
			 * calculate probabilities of Observed trio genotypes given specified true relationships
			 * above
			 * These are rearranged and used to sample trio genotypes given Mi or no MI
			 * Indices are: locus, poten gp1 genotype, poten gp2 geno, grandchild geno
			 * Value is probability NOT log
			 */
			vector <vector <vector <vector <double> > > > prob_trio_genos_OBS;
			// scoped to save memory but no need to make separate function
			{
				for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ // for each locus
					// then calc probabilities of all potential trio genotypes (TRUE genos, not OBS)
					// initialize with 0's
					vector <vector <vector <double> > > tempLoc_true_probs (genotypeKeyC[i].size(), 
                       vector <vector <double> > (genotypeKeyC[i].size(), 
                               vector <double> (genotypeKeyC[i].size(), 0)));
					
					for(int gp1 = 0, mg = genotypeKeyC[i].size(); gp1 < mg; gp1++){
						for(int gp2 = gp1; gp2 < mg; gp2++){
							for(int d = 0; d < mg; d++){
								// P(gp1)*P(gp2)*P(d)
								tempLoc_true_probs[gp1][gp2][d] = exp(lGenos_base[i][gp1] + 
									lGenos_base[i][gp2] + lGenos_randomDescendant[i][d]);
								if(gp1 != gp2) tempLoc_true_probs[gp2][gp1][d] = tempLoc_true_probs[gp1][gp2][d];
							}
						}
					}
					
					// then for each OBS genotype, calc probability given probs of true genotypes and err rates
					// initialize with 0's
					prob_trio_genos_OBS.push_back(vector <vector <vector <double> > > (genotypeKeyC[i].size(),
                                                      vector <vector <double> > (genotypeKeyC[i].size(),
                                                              vector <double> (genotypeKeyC[i].size(), 0))));
					for(int gp1 = 0, mg = genotypeKeyC[i].size(); gp1 < mg; gp1++){
						for(int gp2 = gp1; gp2 < mg; gp2++){
							for(int d = 0; d < mg; d++){
								// for each possible true genotype
								for(int tgp1 = 0; tgp1 < mg; tgp1++){
									for(int tgp2 = 0; tgp2 < mg; tgp2++){
										for(int td = 0; td < mg; td++){
											prob_trio_genos_OBS[i][gp1][gp2][d] += tempLoc_true_probs[tgp1][tgp2][td] * 
												genotypeErrorRatesC[i][tgp1][gp1] * genotypeErrorRatesC[i][tgp2][gp2] * 
												genotypeErrorRatesC[i][td][d];
										}
									}
								} // end for tgp1
								if(gp1 != gp2) prob_trio_genos_OBS[i][gp2][gp1][d] = prob_trio_genos_OBS[i][gp1][gp2][d];
							}
						}
					}
				} // end for each locus
			} // end scoping

			/*
			 * rearranging prob_trio_genos_OBS to be easier to use for sampling trio genotypes
			 * So for each locus, combo has the index of each trio Genotype combination
			 * probability has the probability of each combo
			 * trioGenosLookup has a vector giving the three genotypes in each combo (gpa, gma, d)
			 */
			vector <vector <int> > combo; // index of the genotype combos
			vector <vector <double> > prob_all; // the probabilities of each genotype combo
			vector <vector <double> > prob_MI; // the probabilities of each MI genotype combo
			vector <vector <double> > prob_no_MI; // the probabilities of each non-MI genotype combo
			vector <vector <vector <int> > > trioGenosLookup; // the genotype combos
			for(int i = 0; i < nLoci; i++){
				vector <int> tempCombo;
				vector <double> tempProb_all;
				vector <double> tempProb_MI;
				vector <double> tempProb_no_MI;
				vector <vector <int> > tempGenos;
				int trioIndex = 0;
				double MI_sum = 0;
				double no_MI_sum = 0;
				for(int gpa = 0, max = prob_trio_genos_OBS[i].size(); gpa < max; gpa++){
					for(int gma = 0; gma < max; gma++){
						for(int d = 0; d < max; d++){
							tempProb_all.push_back(prob_trio_genos_OBS[i][gpa][gma][d]);
							// split probs into MI and no MI
							if(noAllelesInCommonGP(genotypeKeyC[i][gpa], genotypeKeyC[i][gma], 
                              genotypeKeyC[i][d])){
								// there is an MI in the genotype
								tempProb_MI.push_back(prob_trio_genos_OBS[i][gpa][gma][d]);
								tempProb_no_MI.push_back(0);
								MI_sum += prob_trio_genos_OBS[i][gpa][gma][d];
							} else {
								tempProb_MI.push_back(0);
								tempProb_no_MI.push_back(prob_trio_genos_OBS[i][gpa][gma][d]);
								no_MI_sum += prob_trio_genos_OBS[i][gpa][gma][d];
							}
							tempCombo.push_back(trioIndex);
							vector <int> tempVec (3);
							tempVec[0] = gpa;
							tempVec[1] = gma;
							tempVec[2] = d;
							tempGenos.push_back(tempVec);
							trioIndex++;
						}
					}
				}
				// normalize
				for(int j = 0, max = tempCombo.size(); j < max; j++){
					tempProb_MI[j] /= MI_sum;
					tempProb_no_MI[j] /= no_MI_sum;
					tempProb_all[j] /= (MI_sum + no_MI_sum);
				}
				combo.push_back(tempCombo);
				prob_all.push_back(tempProb_all);
				prob_MI.push_back(tempProb_MI);
				prob_no_MI.push_back(tempProb_no_MI);
				trioGenosLookup.push_back(tempGenos);
			}
			
			// forward algorithm to calculate probabilities of all 
			//   states (#MI, #miss poten gp1, #miss poten gp2, #miss d)
			
			vector <double> prob_missing_geno (nLoci, 0); // probability that a genotype at locus i is missing
			for(int i = 0; i < nLoci; i++) prob_missing_geno[i] = missingParamsC[i][0] / 
				(missingParamsC[i][0] + missingParamsC[i][1]);
			// calculate probability of observing an MI at each locus given that all genotypes are present
			vector <double> pMI_noMiss (nLoci, 0);
			calcProbMIperLocus_obs(genotypeKeyC, prob_trio_genos_OBS, pMI_noMiss);
			
			// save the states vector for each step
			vector <vector <vector <vector <vector <double> > > > > all_states;
			fwd_ssGP(miLimit, maxMissingGenos, nLoci, prob_missing_geno, pMI_noMiss, all_states);
			
			// probability of an individual having more than the allowed number of missing genotypes
			double probFail = 1;
			{
				vector <double> probTotalFail (nLoci + 1, 0);
				// can use this even though it's not MIs
				calcProbSumMI(prob_missing_geno, probTotalFail);
				for(int i = 0; i <= maxMissingGenos; i++) probFail -= probTotalFail[i];
			}
			// probability of one or more individuals having more than the allowed number of missing genotypes
			double probAnyFail = 1 - pow(1 - probFail, 3); // just complement of all three succeeding
			
			// calculate probability of being in each strata, given maximum number of missing genotypes per individual
			vector <double> pNumMI_OBS (miLimit + 1, 0);
			for(int mi = 0; mi <= miLimit; mi++){ // #mi
				for(int jgp1 = 0; jgp1 <= maxMissingGenos ; jgp1++){ // #miss gp1
					for(int jgp2 = 0; jgp2 <= maxMissingGenos; jgp2++){ // #miss gp2
						for(int jd = 0; jd <= maxMissingGenos; jd++){ // #miss d
							pNumMI_OBS[mi] += all_states[nLoci][mi][jgp1][jgp2][jd];
						}
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
				for(int jgp1 = 0; jgp1 <= maxMissingGenos ; jgp1++){ // #miss gp1
					for(int jgp2 = 0; jgp2 <= maxMissingGenos; jgp2++){ // #miss gp2
						for(int jd = 0; jd <= maxMissingGenos; jd++){ // #miss d
							vector <int> tempStateVec (4);
							tempStateVec[0] = mi;
							tempStateVec[1] = jgp1;
							tempStateVec[2] = jgp2;
							tempStateVec[3] = jd;
							stateIndex.push_back(iState);
							stateProb.push_back(all_states[nLoci][mi][jgp1][jgp2][jd]);
							statesToSample.push_back(tempStateVec);
							tempSum += all_states[nLoci][mi][jgp1][jgp2][jd];
							iState++;
						}
					}
				}
				// normalize probabilities
				for(int s = 0, max = stateProb.size(); s < max; s++) stateProb[s] /= tempSum;
				
				// vector to use with sampleC below
				vector <vector <int> > possibleMove (9, vector <int> (4, 0)); // all possible moves between states
				// order is mi, miss gp1, miss gp2, miss d
				vector <int> moveIndex (9);
				for(int i = 0; i < 9; i++) moveIndex[i] = i;
				possibleMove[1][0] = 1;
				int tempi = 2;
				for(int m1 = 0; m1 < 2; m1++){ // m1-3: 0 is not missing, 1 is missing
					for(int m2 = 0; m2 < 2; m2++){
						for(int m3 = 0; m3 < 2; m3++){
							if(m1 + m2 + m3 == 0) continue; // already handled no missing above
							possibleMove[tempi][1] = m1;
							possibleMove[tempi][2] = m2;
							possibleMove[tempi][3] = m3;
							tempi++;
						}
					}
				}
				
				for(int n = 0; n < itersPerMIC[mi]; n++){
					if(n % 100 == 0) Rcpp::checkUserInterrupt();
					
					// simulate trio - backwards algorithm
					vector <int> gmaGenos (nLoci);
					vector <int> gpaGenos (nLoci);
					vector <int> dGenos (nLoci);
					
					// first choose state (given number of mi)
					vector <int> curState (statesToSample[sampleC(stateIndex, stateProb, rg)]);
					// choose loci with MIs and missing genotypes
					for(int i = nLoci; i > 0; i--){
						vector <double> moveProb (9, 0); // probabilities of moves
						// no missing, no MI (probability you were in curState, and you stayed there)
						moveProb[0] = all_states[i-1][curState[0]][curState[1]][curState[2]][curState[3]] * 
							pow(1 - prob_missing_geno[i-1], 3) * (1 - pMI_noMiss[i-1]);
						// no missing, MI
						if(curState[0] > 0) moveProb[1] = all_states[i-1][curState[0] - 1][curState[1]][curState[2]][curState[3]] * 
							pow(1 - prob_missing_geno[i-1], 3) * pMI_noMiss[i-1];
						
						// all missing options
						tempi = 2;
						for(int m1 = 0; m1 < 2; m1++){ // m1-3: 0 is not missing, 1 is missing
							if(curState[1] - m1 < 0) continue;
							for(int m2 = 0; m2 < 2; m2++){
								if(curState[2] - m2 < 0){
									tempi += 2;
									continue;
								}
								for(int m3 = 0; m3 < 2; m3++){
									if(curState[3] - m3 < 0) {
										tempi++;
										continue;
									}
									if(m1 + m2 + m3 == 0) continue; // already handled no missing above
									moveProb[tempi] = all_states[i-1][curState[0]][curState[1] - m1][curState[2] - m2][curState[3] - m3] *
										pow(prob_missing_geno[i-1], m1 + m2 + m3) * 
										pow(1 - prob_missing_geno[i-1], 3 - (m1 + m2 + m3));
									tempi++;
								}
							}
						}
						
						// select the move made at this locus
						vector <int> move (possibleMove[sampleC(moveIndex, moveProb, rg)]);
						
						// sample genotypes given move
						vector <int> chosenGenos (3);
						if(move[0] == 1){
							// if MI, sample from MI genos
							chosenGenos = trioGenosLookup[i-1][sampleC(combo[i-1], prob_MI[i-1], rg)];
						} else if (move[1] + move[2] + move[3] == 0){
							// if all observed and no MI, sample from non MI genos
							chosenGenos = trioGenosLookup[i-1][sampleC(combo[i-1], prob_no_MI[i-1], rg)];
						} else {
							// if some are missing, sample from all genos, then add missing
							chosenGenos = trioGenosLookup[i-1][sampleC(combo[i-1], prob_all[i-1], rg)];
							if(move[1] == 1) chosenGenos[0] = -9;
							if(move[2] == 1) chosenGenos[1] = -9;
							if(move[3] == 1) chosenGenos[2] = -9;
						}
						
						gpaGenos[i-1] = chosenGenos[0];
						gmaGenos[i-1] = chosenGenos[1];
						dGenos[i-1] = chosenGenos[2];
						
						// change current state for next locus
						for(int j = 0; j < 4; j++) curState[j] -= move[j];
						
					} // end simulating trio
					
					// calc LLR
					
					// no need to filter by MI b/c only simulating valid MI
					
					// passed the MI filter,so calculate the llr
					double uLLH = 0.0; // log-likelihood unrelated
					double gpLLH = 0.0; // log-likelihood grandparent pair
					for(int j = 0; j < nLoci; j++){ //for each locus
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
					
					double llr = gpLLH - uLLH;
					// determine number under and over threshold
					for(int i = 0, max = llrToTestC.size(); i < max; i++){
						if(llr >= llrToTestC[i]) falsePosPerStrata[mi][i]++;
					}
				} // end for n
			} // end for mi

			Rcpp::NumericVector col_pop;
			Rcpp::NumericVector col_pop2;
			Rcpp::NumericVector col_llr;
			Rcpp::NumericVector col_mi;
			Rcpp::NumericVector col_p;
			Rcpp::NumericVector col_fp;
			Rcpp::NumericVector col_fpSD;
			for(int i = 0, ttmax = llrToTestC.size(); i < ttmax; i++){
				for(int j = 0, tmax = falsePosPerStrata.size(); j < tmax; j++){
					col_pop.push_back(pop);
					col_pop2.push_back(pop2);
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
			fullOut.push_back(Rcpp::DataFrame::create(Rcpp::Named("Pop_descendant") = col_pop,
                                             Rcpp::Named("Pop_baseline") = col_pop2,
                                             Rcpp::Named("llrThreshold") = col_llr,
                                             Rcpp::Named("count_MI") = col_mi,
                                             Rcpp::Named("prob_count_MI") = col_p,
                                             Rcpp::Named("falsePos") = col_fp,
                                             Rcpp::Named("falsePosSD") = col_fpSD
			));
			
			for(int i = 0, max2 = llrToTestC.size(); i < max2; i++){
				// calculate the mean accross strata
				double strata_fpMean = 0;
				for(int j = 0, tmax = falsePosPerStrata.size(); j < tmax; j++) if(itersPerMIC[j] > 
                     0) strata_fpMean += pNumMI_OBS[j] * (falsePosPerStrata[j][i] / itersPerMIC[j]);
				popResults.push_back(pop);
				popAssign.push_back(pop2);
				llrRecord.push_back(llrToTestC[i]);
				falsePosAll.push_back(strata_fpMean);
				SDfalsePosAll.push_back(sqrt(fp_SD_total[i]));
				
			}
		} // end for pop2
	} // end for pop
	
	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Pop_descendant") = popResults,
                                                   Rcpp::Named("Pop_baseline") = popAssign,
                                                   Rcpp::Named("llrThreshold") = llrRecord,
                                                   Rcpp::Named("falsePos") = falsePosAll,
                                                   Rcpp::Named("falsePosSD") = SDfalsePosAll
	);
	fullOut.push_front(results);
	return fullOut;
}



// old IS function, uses baseline trios as sampling distribution
// 
// Rcpp::DataFrame IS_otherPopERRORssGP(Rcpp::List baselineParams,
//                            Rcpp::List unsampledPopParams, Rcpp::List missingParams,
//                            Rcpp::List genotypeKey,
//                            Rcpp::List genotypeErrorRates, Rcpp::NumericVector llrToTest,
//                            int N,
//                            int seed, Rcpp::NumericVector skipBaseline
// ){
// 	//////////////////////
// 	// first, turn all inputs into non-Rcpp c++ classes
// 
// 	vector <vector <vector <double> > > baselineParamsC; //first index is pop, then locus, then allele
// 	for(int i = 0, max = baselineParams.length(); i < max; i++){
// 		baselineParamsC.push_back(rcppListToVectorDouble(baselineParams[i]));
// 	}
// 
// 	vector <vector <vector <double> > > unsampledPopParamsC; //first index is pop, then locus, then allele
// 	for(int i = 0, max = unsampledPopParams.length(); i < max; i++){
// 		unsampledPopParamsC.push_back(rcppListToVectorDouble(unsampledPopParams[i]));
// 	}
// 
// 	vector <vector <double> > missingParamsC; // first index is locus, then missing and not missing parameter
// 	missingParamsC = rcppListToVectorDouble(missingParams);
// 
// 	vector <vector <vector <int> > > genotypeKeyC; //first index is locus, then genotype, values are alleles
// 	genotypeKeyC = convertGenotypeKey(genotypeKey);
// 
// 	vector< vector < vector <double> > > genotypeErrorRatesC; // first index is locus, then actual genotype, then observed
// 	for(int i = 0, max = genotypeErrorRates.length(); i < max; i++){
// 		genotypeErrorRatesC.push_back(rcppMatrixToVectorDouble(genotypeErrorRates[i]));
// 	}
// 
// 	vector <double> llrToTestC;
// 	for(int i = 0, max = llrToTest.length(); i < max; i++) llrToTestC.push_back(llrToTest[i]);
// 	
// 	vector <int> skipBaselineC;
// 	for(int i = 0, max = skipBaseline.length(); i < max; i++) skipBaselineC.push_back(skipBaseline[i]);
// 	
// 	// end conversion of types
// 	////////////////////
// 
// 	int nLoci = genotypeKeyC.size();
// 	
// 	// initiate random number generator
// 	mt19937 rg (seed);
// 
// 	// random number generators for missing genotypes for each locus
// 	// simple Bernoulli
// 	vector <bernoulli_distribution> randMissing;
// 	for(int i = 0; i < nLoci; i++){ 
// 		bernoulli_distribution tempDist (missingParamsC[i][0] / (missingParamsC[i][0] + missingParamsC[i][1]));
// 		randMissing.push_back(tempDist);
// 	}
// 	
// 	vector <int> popResults; // pop that trios are simulated from
// 	vector <int> popAssign; // pop that offspring are assigned to
// 	vector <double> llrRecord;
// 	vector <double> falsePosAll;
// 	vector <double> SDfalsePosAll;
// 
// 	// for each baseline population
// 	for(int pop = 0, max = baselineParamsC.size(); pop < max; pop++){
// 		Rcpp::Rcout<<"Now considering descendants from baseline population "<<pop + 1<<"\n";
// 		
// 
// 		/*
// 		 * calculate genotype likelihoods for a descendant of this population sampled at random
// 		 * Index 1 is locus
// 		 * Index 2 is genotype
// 		 * Value is LOG-likelihood
// 		 */
// 		vector <vector <double> > lGenos_randomDescendant;
// 		for(int i = 0; i < nLoci; i++){
// 			vector <double> tempLH;
// 			for(int j = 0, max2 = genotypeKeyC[i].size(); j < max2; j++){
// 				vector <double> k1 (baselineParamsC[pop][i].size(), 0);
// 				vector <double> k2 (baselineParamsC[pop][i].size(), 0);
// 				k1[genotypeKeyC[i][j][0]]++;
// 				k2[genotypeKeyC[i][j][1]]++;
// 				if(genotypeKeyC[i][j][0] == genotypeKeyC[i][j][1]){ // alleles are the same
// 					// prob sample allele from both
// 					tempLH.push_back(logMultPMF(k1, baselineParamsC[pop][i]) + 
// 						logMultPMF(k2, unsampledPopParamsC[pop][i])
// 					);
// 				} else {
// 					// prob sample allele1 from one, allele2 from other + prob of vice versa
// 					tempLH.push_back(log(
// 						exp(logMultPMF(k1, baselineParamsC[pop][i]) + 
// 						logMultPMF(k2, unsampledPopParamsC[pop][i])) + 
// 						exp(logMultPMF(k2, baselineParamsC[pop][i]) + 
// 						logMultPMF(k1, unsampledPopParamsC[pop][i]))
// 					));
// 				}
// 			}
// 			lGenos_randomDescendant.push_back(tempLH);
// 		}
// 
// 		
// 		for(int pop2 = 0, maxp2 = baselineParamsC.size(); pop2 < maxp2; pop2++){
// 			if(pop2 == pop) continue;
// 			bool skipBaseBool = false;
// 			for(int i = 0, maxSkip = skipBaselineC.size(); i < maxSkip; i++){
// 				if(pop2 == skipBaselineC[i]){
// 					skipBaseBool = true;
// 					break;
// 				}
// 			}
// 			if (skipBaseBool) continue;
// 			
// 			Rcpp::Rcout<<"Calculating error with baseline population "<<pop2 + 1<<"\n";
// 			// calculate genotype log-likelihoods from random sampled genotype
// 			// first index is locus, second is genotype, value is LOG-likelihood
// 			vector <vector <double> > lGenos_base;
// 			vector <vector <double> > lGenos_unsamp;
// 	
// 			for(int i = 0, max3 = genotypeKeyC.size(); i < max3; i++){ // for each locus
// 				// baselineParamsC[pop2][i] are the alpha values for the baseline pop
// 				// unsampledPopParamsC[pop2][i] are the alpha values for the corresponding unsampled pop
// 				// genotypeKeyC[i] is the key for the locus i
// 				vector <double> tempBase;
// 				vector <double> tempUnsamp;
// 				for(int j = 0, max2 = genotypeKeyC[i].size(); j < max2; j++){ // for each genotype
// 					vector <double> k (baselineParamsC[pop2][i].size(), 0);
// 					k[genotypeKeyC[i][j][0]]++;
// 					k[genotypeKeyC[i][j][1]]++;
// 					tempBase.push_back(logMultPMF(k, baselineParamsC[pop2][i]));
// 					tempUnsamp.push_back(logMultPMF(k, unsampledPopParamsC[pop2][i]));
// 				}
// 				lGenos_base.push_back(tempBase);
// 				lGenos_unsamp.push_back(tempUnsamp);
// 			}
// 	
// 			/*	
// 			 *	calculate genotype likelihoods for grandparent(maternal) - grandparent(maternal) - grandchild relationship
// 				given true genotypes
// 				different for each pop b/c it incorporates allele frequencies for both the baseline pop and the unsampled pop
// 			 	Index 1 is locus
// 			 	Index 2 is grandparent1 genotype
// 				Index 3 is grandparent 2 genotype
// 				Index 4 is grandchild genotype
// 				Value is likelihood (NOT log)
// 			 */
// 			vector <vector <vector <vector <double> > > > lGenos_ssGP;
// 			createSSGPvector(genotypeKeyC, lGenos_base, unsampledPopParamsC, pop2, lGenos_ssGP);
// 			
// 			/*	
// 			 *	calculate genotype LOG-likelihoods for grandparent(maternal) - grandparent(maternal) - grandchild relationship
// 			 *	and for all three being unrelated
// 				given OBSERVED genotypes - but only if NO MISSING genotypes in the trio
// 				different for each pop b/c it incorporates allele frequencies for both the baseline pop and the unsampled pop
// 			 	Index 1 is locus
// 			 	Index 2 is grandparent1 genotype
// 				Index 3 is grandparent 2 genotype
// 				Index 4 is grandchild genotype
// 				Value is LOG-likelihood
// 			 */
// 			vector <vector <vector <vector <double> > > > lGenos_ssGP_OBS;
// 			vector <vector <vector <vector <double> > > > lGenos_Unrelated_OBS;
// 			createOBSvector(lGenos_ssGP, genotypeErrorRatesC, lGenos_ssGP_OBS, 
// 	                  lGenos_base, lGenos_unsamp, lGenos_Unrelated_OBS);
// 			
// 			/*
// 			 * rearranging lGenos_ssGP_OBS to be easier to use for sampling trio genotypes
// 			 * So for each locus, combo has the index of each trio Genotype combination
// 			 * probability has the probability of each combo
// 			 * trioGenosLookup has a vector giving the three genotypes in each combo (gpa, gma, d)
// 			 */
// 			vector <vector <int> > combo; // index of the genotype combos
// 			vector <vector <double> > prob; // the probabilities of each genotype combo
// 			vector <vector <vector <int> > > trioGenosLookup; // the genotype combos
// 			for(int i = 0; i < nLoci; i++){
// 				vector <int> tempCombo;
// 				vector <double> tempProb;
// 				vector <vector <int> > tempGenos;
// 				int trioIndex = 0;
// 				for(int gpa = 0, maxGP = lGenos_ssGP_OBS[i].size(); gpa < maxGP; gpa++){
// 					for(int gma = 0; gma < maxGP; gma++){
// 						for(int d = 0; d < maxGP; d++){
// 							tempProb.push_back(exp(lGenos_ssGP_OBS[i][gpa][gma][d]));
// 							tempCombo.push_back(trioIndex);
// 							vector <int> tempVec (3,-9);
// 							tempVec[0] = gpa;
// 							tempVec[1] = gma;
// 							tempVec[2] = d;
// 							tempGenos.push_back(tempVec);
// 							trioIndex++;
// 						}
// 					}
// 				}
// 				combo.push_back(tempCombo);
// 				prob.push_back(tempProb);
// 				trioGenosLookup.push_back(tempGenos);
// 			}
// 			
// 			/*
// 			 * making reference "table" of correction ratio lGenos based on observed genotypes
// 			 * with NO missing data
// 			 * Indices are: locus, gp1 observed geno, gp2 observed geno, d observed geno
// 			 * value is LOG-likelihood
// 			 */
// 			vector <vector <vector <vector <double> > > > CORR_lGenos_OBS;
// 			create_CORR_OBSvector(genotypeErrorRatesC, lGenos_randomDescendant, lGenos_base, CORR_lGenos_OBS);
// 			
// 			// inititate results storage for this pop - one entry for each llr
// 			vector <double> falsePos (llrToTestC.size(), 0);
// 			vector <double> SSfalsePos (llrToTestC.size(), 0); // sum of squares for variance calcs
// 
// 			for(int n = 0; n < N; n++){
// 				if(n % 100 == 0) Rcpp::checkUserInterrupt();
// 				
// 				// simulate trio
// 				vector <int> gmaGenos (nLoci, 0);
// 				vector <int> gpaGenos (nLoci, 0);
// 				vector <int> dGenos (nLoci, 0);
// 				vector <int> sampTrio (3, -1);
// 				for(int i = 0; i < nLoci; i++){ // for each locus
// 					// sample observed genotypes
// 					sampTrio = trioGenosLookup[i][sampleC(combo[i], prob[i], rg)]; // sample a trio of observed genotypes
// 	
// 					// add missing genotypes as appropriate
// 					for(int j = 0; j < 3; j++) if(randMissing[i](rg)) sampTrio[j] = -9;
// 					gpaGenos[i] = sampTrio[0];
// 					gmaGenos[i] = sampTrio[1];
// 					dGenos[i] = sampTrio[2];
// 				}
// 				
// 				// calc LLR
// 				double uLLH = 0.0;
// 				double gpLLH = 0.0;
// 				double CORR_uLLH = 0.0; // for correction ratio for importance sampling - this is prob of sampling "grandparents" from 
// 												// the current "baseline" pop and "descendant" of the current descendant pop
// 				for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ //for each locus
// 					// observed genotypes
// 					int obs_gp1 = gpaGenos[j];
// 					int obs_gp2 = gmaGenos[j];
// 					int obs_d = dGenos[j];
// 	
// 					// if missing genotypes
// 					if(obs_gp1 == -9 || obs_gp2 == -9 || obs_d == -9){
// 						double u_likelihood = 0.0; // likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
// 						double gp_likelihood = 0.0;
// 						double CORR_u_likelihood = 0.0;
// 						double pOA_gp1, pOA_gp2, pOA_d; // prob of observed given actual
// 						for(int gp1 = 0, max3 = genotypeKeyC[j].size(); gp1 < max3; gp1++){ // for each possible gp1 genotype
// 							for(int gp2 = 0; gp2 < max3; gp2++){ // for each possible gp2 genotype
// 								for(int d = 0; d < max3; d++){ // for each possible d genotype
// 									if(obs_gp1 == -9){
// 										pOA_gp1 = 1.0;
// 									} else {
// 										pOA_gp1 = genotypeErrorRatesC[j][gp1][obs_gp1];
// 									}
// 									if(obs_gp2 == -9){
// 										pOA_gp2 = 1.0;
// 									} else {
// 										pOA_gp2 = genotypeErrorRatesC[j][gp2][obs_gp2];
// 									}
// 									if(obs_d == -9){
// 										pOA_d = 1.0;
// 									} else {
// 										pOA_d = genotypeErrorRatesC[j][d][obs_d];
// 									}
// 									// unrelated
// 									// P(gp1|baseline pop)P(gp2|baseline pop)P(d|unsampled pop)P(gp1|obs_gp1)P(gp2|obs_gp2)P(d|obs_d)
// 									u_likelihood += exp(lGenos_base[j][gp1] + lGenos_base[j][gp2] + lGenos_unsamp[j][d]) * 
// 										 pOA_gp1 * pOA_gp2 * pOA_d;
// 									// grandparent pair
// 									gp_likelihood += lGenos_ssGP[j][gp1][gp2][d] *
// 										pOA_gp1 * pOA_gp2 * pOA_d;
// 									// correction for importance sampling
// 									CORR_u_likelihood += exp(lGenos_base[j][gp1] + lGenos_base[j][gp2] + lGenos_randomDescendant[j][d]) * 
// 										 pOA_gp1 * pOA_gp2 * pOA_d;
// 								}
// 							}
// 						}
// 						uLLH += log(u_likelihood);
// 						gpLLH += log(gp_likelihood);
// 						CORR_uLLH += log(CORR_u_likelihood);
// 					} else { // no missing genotypes
// 						uLLH += lGenos_Unrelated_OBS[j][obs_gp1][obs_gp2][obs_d];
// 						gpLLH += lGenos_ssGP_OBS[j][obs_gp1][obs_gp2][obs_d];
// 						CORR_uLLH += CORR_lGenos_OBS[j][obs_gp1][obs_gp2][obs_d];
// 					}
// 				}
// 				// determine number under and over threshold
// 				for(int i = 0, max2 = llrToTestC.size(); i < max2; i++){
// 					if((gpLLH - uLLH) >= llrToTestC[i]){
// 						falsePos[i] += exp(CORR_uLLH - gpLLH); // correction for sampling distribution
// 						SSfalsePos[i] += pow(exp(CORR_uLLH - gpLLH), 2);
// 					}
// 				}
// 			} // end for N
// 			for(int i = 0, max2 = llrToTestC.size(); i < max2; i++){
// 				double fpMean = falsePos[i] / N; // need for variance calc
// 
// 				popResults.push_back(pop);
// 				popAssign.push_back(pop2);
// 				llrRecord.push_back(llrToTestC[i]);
// 				falsePosAll.push_back(fpMean);
// 				SDfalsePosAll.push_back(sqrt(varianceEstim(N, SSfalsePos[i], fpMean)));
// 			}
// 		
// 		} // end for pop2
// 			
// 	} // end for pop
// 
// 	// organize as a DataFrame and return	
// 	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Pop_descendant") = popResults,
//                                                    Rcpp::Named("Pop_baseline") = popAssign,
//                                                    Rcpp::Named("llrThreshold") = llrRecord,
//                                                    Rcpp::Named("falsePos") = falsePosAll,
//                                                    Rcpp::Named("falsePosSD") = SDfalsePosAll
// 																	);
// 	
// 	return results;
// }
