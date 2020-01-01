#include <Rcpp.h>
#include <math.h>
#include <vector>

using namespace std;



//' pmf of Dirichlet-multinomial distribution
//' @param n the total number of draws
//' @param k a vector of the number of observations in each category
//' @param a a vector of the parameters (alpha) describing the Dirichlet
//' @noRd
//' @keywords internal

double logDirichMultPMF(vector <double> k, vector <double> a){
	int numAlpha = a.size();
	int numK = k.size(); // just prevent compiler warning b/c I am nitpicky
	if (numK != numAlpha) Rcpp::stop("internal error: k not equal to a in logdirichMultPMF");

	double n = 0;
	for(int i = 0; i < numAlpha; i++) n += k[i];
	double sumAlpha = 0;	
	for(int i = 0; i < numAlpha; i++) sumAlpha += a[i];
	
	//start density calculation
	double dens = lgamma(n + 1) + lgamma(sumAlpha) - lgamma(n + sumAlpha);
	for(int i = 0; i < numAlpha; i++){
		dens += lgamma(k[i] + a[i]) - lgamma(k[i] + 1) - lgamma(a[i]);
	}
	return dens;
}

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

//' caluculate likelihood of parent 1, parent 2, and offspring genotypes given Mendalian inheritance
//' @param p1 parent 1 genotypes (two ints - diploid)
//' @param p2 parent 2 genotypes (two ints - diploid)
//' @param o offspring genotypes (two ints - diploid)
double ppoMendelian(vector <int> p1, vector <int> p2, vector <int> o){
	
	//something messed up here?
	
	bool o0p1 = (o[0] == p1[0] || o[0] == p1[1]); // is o[0] in p1
	bool o0p2 = (o[0] == p2[0] || o[0] == p2[1]); // is o[0] in p2
	bool o1p1 = (o[1] == p1[0] || o[1] == p1[1]); // is o[1] in p1
	bool o1p2 = (o[1] == p2[0] || o[1] == p2[1]); // is o[1] in p2
	if(o0p1 + o1p2 < 2 && o1p1 + o0p2 < 2) return 0;
	
	// just make some Punnet squares!
	double sum = 0;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			if((o[0] = p1[i] && o[1] == p2[j]) || (o[1] = p1[i] && o[0] == p2[j])) sum++;
		}
	}
	return sum / 4;
}

/*

// wrote these and then remembered want to make this work for both
//  SNPs and microhaplotypes, keeping them in for now in case they
//  turn out to be useful later on

// calculate log of the beta function
double logBeta(double a, double b) {
	return lgamma(a) + lgamma(b) - lgamma(a+b);
}

// calculate log of the binomial coefficient (choose)
double logChoose(double n, double k){
	return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
}

// pmf of beta-binomial distribution
// [[Rcpp::export]]
double logBetaBinomPMF(double n, double k, double a, double b){
	return logChoose(n, k) + logBeta(k + a, n - k + b) - logBeta(a, b);
}
*/
