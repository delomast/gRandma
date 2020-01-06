// miscellaneous helper functions

#include <Rcpp.h>
#include <vector>
#include "misc_math.h"
#include "utils.h"

using namespace std;

//' list all possible pairs of individuals from a vector of individual ID's (as integers)
vector <vector <int> > listAllPairs(vector <int> allInds){
	vector <vector <int> > out;
	vector <int> tempVec;
	tempVec.push_back(-9);
	tempVec.push_back(-9);
	for(int i = 0, max = allInds.size() - 1; i < max; i++){
		for(int j = i+1, max2 = allInds.size(); j < max2; j++){
			tempVec[0] = allInds[i];
			tempVec[1] = allInds[j];
			out.push_back(tempVec);
		}
	}
	return out;
}

//' convert a NumericMatrix of ints to to nested vectors, index 1 is row index 2 is column
vector <vector <int> > rcppMatrixToVectorInt(Rcpp::NumericMatrix inputMatrix){
	vector <vector <int> > output;
	for(int i = 0, max = inputMatrix.nrow(); i < max; i++){
		vector <int> tempVec;
		for(int j = 0, max2 = inputMatrix.ncol(); j < max2; j++){
			tempVec.push_back(inputMatrix(i,j));
		}
		output.push_back(tempVec);
	}
	return output;
}

//' convert a List of NumericVectors to nested vectors index 1 is list position, index 2 is vector position
vector <vector <double> > rcppListToVectorDouble(Rcpp::List inputList){
	vector <vector <double> > output;
	for(int i = 0, max = inputList.length(); i < max; i++){
		Rcpp::NumericVector current;
		current = inputList[i];
		vector <double> tempVec;
		for(int j = 0, max2 = current.length(); j < max2; j++){
			tempVec.push_back(current[j]);
		}
		output.push_back(tempVec);
	}
	return output;
}

//' convert the GenotypeKey list to nested vectors
vector <vector <vector <int> > > convertGenotypeKey(Rcpp::List inputList){
	vector <vector <vector <int> > > output;
	for(int i = 0, max = inputList.length(); i < max; i++){ // for each locus
		Rcpp::NumericMatrix current = inputList[i];
		vector <vector <int> > tempLocus;
		for(int j = 0, max2 = current.nrow(); j < max2; j++){ // for each genotype
			vector <int> tempGeno (2, -9);
			tempGeno[0] = current(j,1);
			tempGeno[1] = current(j,2);
			tempLocus.push_back(tempGeno);
		}
		output.push_back(tempLocus);
	}
	return output;
}

//' convert a NumericMatrix to a nexted vector with first index being the row and 
//'   second index being the column
vector < vector <double> > rcppMatrixToVectorDouble(Rcpp::NumericMatrix inputMatrix){
	vector < vector <double> > output;
	for(int i = 0, max = inputMatrix.nrow(); i < max; i++){
		vector <double> tempVec;
		for(int j = 0, max2 = inputMatrix.ncol(); j < max2; j++){
			tempVec.push_back(inputMatrix(i,j));
		}
		output.push_back(tempVec);
	}
	return output;
}

//' takes genotypeKeyC, calculates Parent parent offspring likelihoods for all genotype combinations
//'   and then retuns it as a nested vector with :	 	Index 1 is locus, Index 2 is parent1 genotype
//'   Index 3 is parent 2 genotype, Index 4 is offspring genotype, Value is likelihood
vector <vector <vector <vector <double> > > > createPPOvector(vector <vector <vector <int> > > genotypeKeyC){
	
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
	
	return lGenos_ppo;
}

//' dertermine whether there are any alleles in common between the four alleles in
//'   two potential grandparents (gp1 and gp2) and a potential descendant (d). This
//'   is only for one locus. Returns true if no allles are in common (ie a MI),
//'   and false if there are one or more alleles in common.
//' @param gp1 genotype as vector of alleles
//' @param gp2 genotype as vector of alleles
//' @param d genotype as vector of alleles

bool noAllelesInCommonGP(vector <int> gp1, vector <int> gp2, vector <int> d){
	if(d[0] == gp1[0] || d[0] == gp1[1] || d[1] == gp1[0] || d[1] == gp1[1] ||
    d[0] == gp2[0] || d[0] == gp2[1] || d[1] == gp2[0] || d[1] == gp2[1]){
		return false;
	} else {
		return true;
	}
}

