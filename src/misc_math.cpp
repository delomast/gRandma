// miscellaneous functions performing calculations
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <random>
#include "misc_math.h"

using namespace std;

/*
// pmf of Dirichlet-multinomial distribution
// @param k a vector of the number of observations in each category
// @param a a vector of the parameters (alpha) describing the Dirichlet
 
double logDirichMultPMF(const vector <double>& k, const vector <double>& a){
	if (k.size() != a.size()) Rcpp::stop("internal error: k not equal to a in logdirichMultPMF");

	int numAlpha = a.size();
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

*/

// pmf of multinomial distribution
// @param k a vector of the number of observations in each category
// @param a a vector of the parameters (alpha) describing relative probabilities of each group
double logMultPMF(const vector <double>& k, const vector <double>& a){
	if (k.size() != a.size()) Rcpp::stop("internal error: k not equal to a in logMultPMF");

	int numA = a.size();
	double n = 0;
	for(int i = 0; i < numA; i++) n += k[i];
	double l_sumA = 0;	
	for(int i = 0; i < numA; i++) l_sumA += a[i];
	l_sumA = log(l_sumA);

	//start density calculation
	double dens = lgamma(n + 1);
	for(int i = 0; i < numA; i++){
		dens += (k[i] * (log(a[i]) - l_sumA)) - lgamma(k[i] + 1);
	}
	return dens;
}
 
 
// caluculate likelihood of parent 1, parent 2, and offspring genotypes given Mendalian inheritance
// @param p1 parent 1 genotypes (two ints - diploid)
// @param p2 parent 2 genotypes (two ints - diploid)
// @param o offspring genotypes (two ints - diploid)
 
double ppoMendelian(const vector <int>& p1, const vector <int>& p2, const vector <int>& o){

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

// random sample from categorical distribution
// items are items to select from, represented by ints
// probs are relative probabilities (normalized within the function)

int sampleC(const vector <int>& items, const vector <double>& probs, mt19937& rNum) {
	
	// input checking and remove items with 0 probability
	if(items.size() != probs.size()) Rcpp::stop("invalid input to sampleC, dimensions differ");
	
	vector <double> newProbs;
	vector <int> newItems;
	
	// now change to section of the interval from 0-cumulSum and remove items with prob of 0
	double cumulSum = 0;
	for(int i=0, max = probs.size(); i<max; i++){
		if(probs[i] < 0) Rcpp::stop("internal error: invalid probs to sampleC");
		if(probs[i] > 0){
			cumulSum += probs[i];
			newProbs.push_back(cumulSum);
			newItems.push_back(items[i]);
		}
	}
	if (newProbs.size() == 0) Rcpp::stop("internal error: all 0 probs to sampleC");
	// end input check
	
	// get random uniform
	uniform_real_distribution <double> rUnif(0.0,cumulSum);
	double rN = rUnif(rNum);
	// assign result
	for(int i=0, max=newProbs.size(); i<max; i++){
		if(rN < newProbs[i]){
			return newItems[i];
		}
	}
	// Rcpp::Rcout << rN << "\n";
	// Rcpp::Rcout << numeric_limits<double>::min() << "\n";
	// Rcpp::Rcout << numeric_limits<double>::denorm_min() << "\n";
	// for(int i=0, max=newProbs.size(); i<max; i++) Rcpp::Rcout << newProbs[i] << "\n";
	Rcpp::stop("internal error: in sampleC, no result chosen");
	return 0;
}

// n is number of samples
// ss is sum of squares of samples
// mean is mean of samples
double varianceEstim(double n, double ss, double mean){
	return ((ss / (n - 1)) - ((pow(mean, 2) * n) / (n - 1))) / n;
}
