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
//' Stratified sampling Monte Carlo used for estimating false positive 
//' 
//' 
//' This version uses trios from the current "baseline" pop as the importance sampling distribution
//' 
//' @param baselineParams Dirichlet parameters for allele frequencies
//' @param unsampledPopParams Dirichlet parameters for allele frequencies
//' @param missingParams Beta parameters for missing genotypes (failure to genotype rate)
//' @param genotypeKey list of matrix for each locus, col1 is genotype, col2 is allele 1, col3 is allele 2
//' @param genotypeErrorRates list of matrix for each locus, rows are actual genotype, columns are observed,
//'   values are probability
//' @param llrToTest Vector of llr's to test as threshold values
//' @param N number of samples to take
//' @param skipBaseline added unsampled pops to skip as baseline
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::List strat_otherPopERRORsP(Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::List missingParams,
                           Rcpp::List genotypeKey,
                           Rcpp::List genotypeErrorRates, Rcpp::NumericVector llrToTest,
                           Rcpp::NumericVector itersPerMI,
                           int seed, Rcpp::NumericVector skipBaseline
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

	vector <double> itersPerMIC;
	for(int i = 0, max = itersPerMI.length(); i < max; i++) itersPerMIC.push_back(itersPerMI[i]);
	
	vector <int> skipBaselineC;
	for(int i = 0, max = skipBaseline.length(); i < max; i++) skipBaselineC.push_back(skipBaseline[i]);
	
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
	
	vector <int> popResults; // pop that descendants are simulated from
	vector <int> popAssign; // pop that baseline is simulated from
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
			if(skipBaseBool) continue;
			
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

			// first calculate the probability of MI at each locus | unrelated desc from pop and baseline from pop2
			vector <double> pMI (genotypeKeyC.size(), 0.0); // p(MI) observed for each locus
			calcProbMIperLocus_Unrel(genotypeKeyC, genotypeErrorRatesC, lGenos_base, lGenos_randomDescendant, pMI);
			// for(int i=0; i < pMI.size(); i++) Rcpp::Rcout << pMI[i] << " "; // testing
			// Rcpp::Rcout << "\n\n\n";
			
			// then calculate probability of sum(MI) = x for x in [0, num_loci]
			vector <double> pTotalMI (genotypeKeyC.size() + 1, -9); // probability of having x observed MI | true parent
			calcProbSumMI(pMI, pTotalMI);
			// for(int i=0; i < pTotalMI.size(); i++) Rcpp::Rcout << i << ": " << pTotalMI[i] << " \n"; // testing
			// Rcpp::Rcout << "\n\n";
			
			/*
			 * calculate probability of pairs of OBSERVED genotypes given MI at a locus or given no MI at a locus 
			 * indices are: locus, baseline genotype, desc genotype
			 * value is probability / likelihood NOT log
			 */
			vector <vector <vector <double> > > lGenos_OBS_sample_MI;
			vector <vector <vector <double> > > lGenos_OBS_sample_no_MI;
			
			strat_lGenosBuilder_sP(genotypeKeyC, genotypeErrorRatesC, lGenos_base,
                        lGenos_randomDescendant,
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
			createSPvector(genotypeKeyC, lGenos_base, unsampledPopParamsC, pop2, lGenos_sP);
	
			/*	
			 *	calculate genotype LOG-likelihoods for parent - grandchild relationship
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
			
			// some results storage for full results reporting
			vector <vector <double> > falsePosPerStrata; // indices: num of MI, llr threshold
			for(int i=0, tmax = itersPerMIC.size(); i < tmax; i++){
				vector <double> tVec (llrToTestC.size(), 0);
				falsePosPerStrata.push_back(tVec);
			}

			// inititate results storage for this pop - one entry for each llr
			vector <double> fp_SD_total (llrToTestC.size(), 0);

			// save a vector of ints to select loci with MI
			vector <int> lociIndex(nLoci, 0);
			for(int i=0; i < nLoci; i++) lociIndex[i] = i;
			
			for(int mi = 0, maxMI = itersPerMIC.size(); mi < maxMI; mi++){
				Rcpp::checkUserInterrupt();
				
				for(int n = 0; n < itersPerMIC[mi]; n++){
					if(n % 100 == 0) Rcpp::checkUserInterrupt();
					
					// simulate pair
					vector <int> p1Genos (nLoci, 0);
					vector <int> dGenos (nLoci, 0);
					// first choose loci with MI
					vector <bool> lociWithMI (nLoci, false);
					vector <double> temp_pMI (pMI);
					for(int i = 0; i < mi; i++){
						int chosen = sampleC(lociIndex, temp_pMI, rg);
						lociWithMI[chosen] = true;
						temp_pMI[chosen] = 0;
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

					// calc LLR
					double uLLH = 0.0;
					double pLLH = 0.0;

					for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ // for each locus
						// observed genotypes
						int obs_p1 = p1Genos[j];
						int obs_d = dGenos[j];
						
						if(obs_d == -9) continue; // no information
						
						// if parent is missing, there is information IF there is an "unsampled pop"
						// because allele freqs will be different whether it is an offspring or not
						if(obs_p1 == -9){
							double u_likelihood = 0.0; // likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
							double p_likelihood = 0.0;
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
					// determine number under and over threshold
					double llr = pLLH - uLLH;
					for(int i = 0, max2 = llrToTestC.size(); i < max2; i++){
						if(llr >= llrToTestC[i]){
							falsePosPerStrata[mi][i]++;
						}
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
					col_p.push_back(pTotalMI[j]);
					if(itersPerMIC[j] == 0){
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
							0) strata_fpMean += pTotalMI[j] * (falsePosPerStrata[j][i] / itersPerMIC[j]);

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
