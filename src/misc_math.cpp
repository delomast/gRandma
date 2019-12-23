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
	if (k.size() != numAlpha) Rcpp::stop("internal error: k not equal to a in logdirichMultPMF");

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
