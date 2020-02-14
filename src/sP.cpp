#include <Rcpp.h>
#include <vector>
#include <math.h>
#include "misc_math.h"
#include "utils.h"
#include "mend_incompat.h"

using namespace std;

//' calculating llr of parent - offspring pair vs unrelated
//' returning a matrix of descendent, grandparentPopulation, grandparent1, grandparent2, llr
//' @param baseline baseline individuals, col1 is pop, col2 is id, cols 3... are genotypes
//' @param mixture mixture individulas col1 is id, cols 2... are genotypes
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
Rcpp::DataFrame sP(Rcpp::NumericMatrix baseline, Rcpp::NumericMatrix mixture, 
                           Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::List genotypeKey,
                           Rcpp::List genotypeErrorRates, double saveLLR, double MIexcludeProb,
                           bool filterLLR
                           ){

	//////////////////////
	// first, turn all inputs into non-Rcpp c++ classes
	// here, making vectors with rows each being a vector of ints
	vector <vector <int> > baselineC;
	baselineC = rcppMatrixToVectorInt(baseline);
	
	vector <vector <int> > mixtureC;
	mixtureC = rcppMatrixToVectorInt(mixture);
	
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
	vector <int> ParentPop;
	vector <int> parent1;
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
			
			// now calculate the probability of sum(MI) = x | true grandparents
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
		} else{
			Rcpp::Rcout<<"MIexcludeProb was not greater than zero, so likelihoods for all possible combinations"<<
				" will be calculated.\n";
		}

		// get list of all possible parents
		vector <int> baseInds;
		for(int i=0, max = baselineC.size(); i < max; i++) {
			if(baselineC[i][0] == pop) baseInds.push_back(baselineC[i][1]);
		}

		Rcpp::Rcout<<"Number of potential parents: "<< baseInds.size()<<"\n\n";
		// for each individual in the mixture
		for(int m = 0, maxM = mixtureC.size(); m < maxM; m++){
			// cout<<m<<endl; //testing
			for(int i = 0, max = baseInds.size(); i < max; i++){ // for each pair
				if(i % 100 == 0) Rcpp::checkUserInterrupt();
				// if(i % 500 == 0) cout<<"pair number: "<<i<<endl;  //testing
				
				// filter based on MI
				int countMI = 0;
				bool skip = false;
				for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ // for each locus
					int p1Geno = baselineC[baseInds[i]][j+2];
					int dGeno = mixtureC[m][j+1];
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
				
				// passed the MI filter,so calculate the llr
				double uLLH = 0.0; // log-likelihood unrelated
				double pLLH = 0.0; // log-likelihood parent
				for(int j = 0, max2 = genotypeKeyC.size(); j < max2; j++){ //for each locus
					// observed genotypes
					int obs_p1 = baselineC[baseInds[i]][j+2];
					int obs_d = mixtureC[m][j+1];

					// if missing genotypes
					if(obs_p1 == -9 || obs_d == -9){
						double u_likelihood = 0.0; // likelihoods for this locus (NOT log) b/c sum across all possible true genotypes
						double p_likelihood = 0.0;
						double pOA_p1, pOA_d; // prob of observed given actual
						for(int p1 = 0, max3 = genotypeKeyC[j].size(); p1 < max3; p1++){ // for each possible p1 genotype
								for(int d = 0; d < max3; d++){ // for each possible d genotype
									if(obs_p1 == -9){
										pOA_p1 = 1.0;
									} else {
										pOA_p1 = genotypeErrorRatesC[j][p1][obs_p1];
									}
									if(obs_d == -9){
										pOA_d = 1.0;
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
								}
							
						}
						uLLH += log(u_likelihood);
						pLLH += log(p_likelihood);
					} else { // no missing genotypes
						uLLH += lGenos_Unrelated_OBS[j][obs_p1][obs_d];
						pLLH += lGenos_sP_OBS[j][obs_p1][obs_d];
					}
				}
				
				// save results
				if((pLLH - uLLH) >= saveLLR || !filterLLR){
					mixtureInd.push_back(mixtureC[m][0]);
					ParentPop.push_back(pop);
					parent1.push_back(baseInds[i]);
					llr.push_back(pLLH - uLLH);
					mendIncompat.push_back(countMI);
				}

			} // end for each parent

		} // end for each mixture individual
	
	} // end for each baseline pop
	

	// organize as a DataFrame and return	
	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Kid") = mixtureInd,
                                                   Rcpp::Named("Pop") = ParentPop,
                                                   Rcpp::Named("Parent") = parent1,
                                                   Rcpp::Named("llr") = llr,
                                                   Rcpp::Named("MI") = mendIncompat
																	);
	
	return results;
}
