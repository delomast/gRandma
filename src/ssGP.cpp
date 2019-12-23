#include <Rcpp.h>
#include <vector>
#include "misc_math.h"

using namespace std;

//' calculating llr of single-sided grandparent pair vs unrelated
//' returning a matrix of descendent, grandparentPopulation, grandparent1, grandparent2, llr
//' @param baseline baseline individuals, col1 is pop, col2 is id, cols 3... are genotypes
//' @param mixture mixture individulas col1 is id, cols 2... are genotypes
//' @param crossRecords col1 is pop, col2 is id 1, col3 is id2
//' @param baselineParams Dirichlet parameters for allele frequencies
//' @param unsampledPopParams Dirichlet parameters for allele frequencies
//' @param genotypeKey col1 is genotype, col2 is allele 1, col3 is allele 2
//' @keywords internal
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix ssGP(Rcpp::NumericMatrix baseline, Rcpp::NumericMatrix mixture, 
                           Rcpp::NumericMatrix crossRecords, Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::NumericMatrix genotypeKey
                           ){
	// determine if unsampledPops used or not
	bool useCR = true;
	if(crossRecords.nrow() == 0) useCR = false;
	
	// calculate number of MI allowed to consider a pair of potential grandparents
	double miLimit = (baseline.ncol() - 2) * .05; // this is just a placeholder for initial testing
	
	// get number of baseline populations
	int countPops = baselineParams.length();
	
	// make map of baseline individuals and their genotypes - or just the row of the matrix they are in
	//  actually, probably don't need to if individual identifier is their row number
	//
	//
	//
	
	// for each baseline pop
	for(int pop = 0, maxP = countPops; pop < maxP; pop++){
		// check if cross records present
		bool useCRpop = false;
		if (useCR){
			// check that cross records contain records for the current population
			for(int i = 0, max = crossRecords.nrow(); i < max; i++){
				if(crossRecords(i,0) == pop){
					useCRpop = true;
					break;
				}
			}
		}
		
		// either get cross records, or make a list of all possible pairs
		vector <vector <int> > pairs; //all pairs to consider
		vector <int> tempVec;
		tempVec.push_back(-9);
		tempVec.push_back(-9);
		if(useCR){
			for(int i = 0, max = crossRecords.nrow(); i < max; i++){
				if(crossRecords(i,0) == pop){
					tempVec[0] = crossRecords(i,1);
					tempVec[1] = crossRecords(i,2);
					pairs.push_back(tempVec);
				}
			}
		} else {
			vector <int> baseInds;
			for(int i=0, max = baseline.nrow(); i < max; i++) {
				if(baseline(i,0) == pop) baseInds.push_back(baseline(i,1));
			}
			
			pairs = listAllPairs(baseInds);
		}
		
		// testing
		Rcpp::NumericMatrix out(pairs.size(), 2);
		for(int i = 0, max = pairs.size(); i < max; i++){
			out(i,0) = pairs[i][0];
			out(i,1) = pairs[i][1];
		}
		return out;
		// end testing
				
		// calculate genotype likelihoods
		//  first from sampling
		//  then from parent - parent - offspring relationship
		
		// for each individual in the mixture
			// filter pairs based on MI
			// for each potential grandparent pair
			// calculate llr
			// save results
		
		
		
	}
	

	
	
	// organize as a matrix and return	
// 	Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("Kid") = 0,
//                                                    Rcpp::Named("Gpop") = 0,
//                                                    Rcpp::Named("Gma") = 0,
//                                                    Rcpp::Named("Gpa") = 0,
//                                                    Rcpp::Named("llr") = 0
// 																	);
// 	return results;
}
