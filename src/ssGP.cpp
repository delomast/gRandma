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
int ssGP(Rcpp::NumericMatrix baseline, Rcpp::NumericMatrix mixture, 
                           Rcpp::NumericMatrix crossRecords, Rcpp::List baselineParams,
                           Rcpp::List unsampledPopParams, Rcpp::List genotypeKey
                           ){

	// determine if cross records used or not
	bool useCR = true;
	if(crossRecords.nrow() == 0) useCR = false;
	
	// first, turn all inputs into non-Rcpp c++ classes
	// here, making vectors with rows each being a vector of ints
	vector <vector <int> > baselineC;
	baselineC = rcppMatrixToVectorInt(baseline);
	
	vector <vector <int> > mixtureC;
	baselineC = rcppMatrixToVectorInt(mixture);
	
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
	

	// calculate number of MI allowed to consider a pair of potential grandparents
	double miLimit = (baselineC[0].size() - 2) * .05; // this is just a placeholder for initial testing

	
		
	// calculate genotype likelihoods from parent - parent - offspring relationship
	// the same for all baseline pops b/c based on mendelian inheritance
	/*	For lGenos_ppo:
	 	Index 1 is locus
	 	Index 2 is parent1 genotype
		Index 3 is parent 2 genotype
		Index 4 is offspring genotype
		Value is likelihood
	 */
	vector <vector <vector <vector <double> > > > lGenos_ppo;
	
	// initialize with 0's
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ //for each locus
		vector <vector <vector <double> > > tempLocus;
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			vector <vector <double> > tempP1;
			for(int p2 = 0; p2 < max2; p2++){ //for each p2 genotype
				vector <double> tempP2 (max2,0.0);
				tempP1.push_back(tempP2);
			}
			tempLocus.push_back(tempP1);
		}
		lGenos_ppo.push_back(tempLocus);
	}

	// now calculate likelihoods
	
	//something messed up here?
	
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ //for each locus
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			for(int p2 = p1; p2 < max2; p2++){ //for each p2 genotype
				for(int o = 0; o < max2; o++){ //for each offspring genotype
					lGenos_ppo[i][p1][p2][o] = ppoMendelian(genotypeKeyC[i][p1], genotypeKeyC[i][p2], genotypeKeyC[i][o]);
					if(p1 != p2) lGenos_ppo[i][p2][p1][o] = lGenos_ppo[i][p1][p2][o];
				}
			}
		}
	}

//testing
	for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ //for each locus
		for(int p1 = 0, max2 = genotypeKeyC[i].size(); p1 < max2; p1++){ //for each p1 genotype
			for(int p2 = 0; p2 < max2; p2++){ //for each p2 genotype
				for(int o = 0; o < max2; o++){ //for each offspring genotype
					cout<<genotypeKeyC[i][p1][0]<<" "<<genotypeKeyC[i][p1][1]<<" "<<genotypeKeyC[i][p2][0]<<" "<<genotypeKeyC[i][p2][1]<<" "<<genotypeKeyC[i][o][0]<<" "<<genotypeKeyC[i][o][1]<<endl;
					cout<<lGenos_ppo[i][p1][p2][o]<<endl;
					cout<<lGenos_ppo[i][p2][p1][o]<<endl;
				}
			}
		}
	}
	
return 1;
	
	// for each baseline pop
	for(int pop = 0, maxP = baselineParamsC.size(); pop < maxP; pop++){
		// check if cross records present
		bool useCRpop = false;
		if (useCR){
			// check that cross records contain records for the current population
			for(int i = 0, max = crossRecordsC.size(); i < max; i++){
				if(crossRecordsC[i][0] == pop){
					useCRpop = true;
					break;
				}
			}
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
		} else {
			vector <int> baseInds;
			for(int i=0, max = baselineC.size(); i < max; i++) {
				if(baselineC[i][0] == pop) baseInds.push_back(baselineC[i][1]);
			}
			pairs = listAllPairs(baseInds);
		}
		

		// calculate genotype likelihoods from random sampled genotype
		vector <vector <double> > lGenos_base; // one vector for each locus
		vector <vector <double> > lGenos_unsamp; 

		for(int i = 0, max = genotypeKeyC.size(); i < max; i++){ // for each locus
			// baselineParamsC[pop][i] are the alpha values for the baseline pop
			// unsampledPopParamsC[pop][i] are the alpha values for the corresponding unsampled pop
			// genotypeKeyC[i] is the key for the locus i
			vector <double> tempBase;
			vector <double> tempUnsamp;
			for(int j = 0, max2 = genotypeKeyC[i].size(); j < max2; j++){ // for each genotype
				vector <double> k (baselineParamsC[pop][i].size(), 0);
				for(int l = 0, max3 = baselineParamsC[pop][i].size(); l < max3; l++){
					if(genotypeKeyC[i][j][0] == l) k[l]++;
					if(genotypeKeyC[i][j][1] == l) k[l]++;
				}
				tempBase.push_back(logDirichMultPMF(k, baselineParamsC[pop][i]));
				tempUnsamp.push_back(logDirichMultPMF(k, unsampledPopParamsC[pop][i]));
			}
			lGenos_base.push_back(tempBase);
			lGenos_unsamp.push_back(tempUnsamp);
		}

		
		
		
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
 
 

	return 1;
}
