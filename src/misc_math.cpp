// miscellaneous functions performing calculations

#include <Rcpp.h>
#include <math.h>
#include <vector>
#include "misc_math.h"

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


//' caluculate likelihood of parent 1, parent 2, and offspring genotypes given Mendalian inheritance
//' @param p1 parent 1 genotypes (two ints - diploid)
//' @param p2 parent 2 genotypes (two ints - diploid)
//' @param o offspring genotypes (two ints - diploid)
double ppoMendelian(vector <int> p1, vector <int> p2, vector <int> o){

	bool o0p1 = (o[0] == p1[0] || o[0] == p1[1]); // is o[0] in p1
	bool o0p2 = (o[0] == p2[0] || o[0] == p2[1]); // is o[0] in p2
	bool o1p1 = (o[1] == p1[0] || o[1] == p1[1]); // is o[1] in p1
	bool o1p2 = (o[1] == p2[0] || o[1] == p2[1]); // is o[1] in p2
	if((o0p1 + o1p2) < 2 && (o1p1 + o0p2) < 2) return 0;
	
	// just make some Punnet squares!
	double sum = 0;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			if((o[0] == p1[i] && o[1] == p2[j]) || (o[1] == p1[i] && o[0] == p2[j])) sum++;
		}
	}
	return sum / 4.0;
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
