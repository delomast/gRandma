#' log-likelihood of relationship of one locus with true genotypes
#' 
#'  
probGP <- function(potKid, potGma, potGpa, altFreq, refFreq){
	llh <- pmfBetaBinomial(potGma, 2, altFreq, refFreq, log = TRUE) + pmfBetaBinomial(potGpa, 2, altFreq, refFreq, log = TRUE)
	probSum <- 0
	for(f in 0:2){
		for(m in 0:2){
			probSum <- probSum +
				pmfBetaBinomial(f, 2, altFreq, refFreq) * condProbTrio(potGma, potGpa, m) * condProbTrio(f, m, potKid)
		}
	}
	llh <- llh + log(probSum)

	return(llh)
}

#' log-likelihood of relationship of one locus with genotyping error
#' 
#' 
probGPerror <- function(potKid, potGma, potGpa, altFreq, refFreq, epsPerAllele){
	lh <- 0
	#cylcling over all possible combinations of true genotypes
	for(k in 0:2){
		for(m in 0:2){
			for(p in 0:2){
				
				lh <- lh + exp(probGP(k, m, p, altFreq, refFreq)) * p_obsG(potKid, k, epsPerAllele) * 
					p_obsG(potGma, m, epsPerAllele) * p_obsG(potGpa, p, epsPerAllele)
				
				
			}
		}
	}
	return(log(lh))
}

#' prob of genotype observed given true and genotyping error rate
#' 
p_obsG <- function(obsG, trueG, epsPerAllele){
	if(obsG == -9) return(1) # to sum over missing data... ?
	if(obsG == 0 || obsG == 2){
		if(trueG == obsG){
			return((1 - epsPerAllele)^2)
		} else if (trueG == 1) {
			return((1-epsPerAllele)*epsPerAllele)
		} else {
			return(epsPerAllele^2)
		}
	} else { #obs is 1
		if(trueG == 1){
			return((1 - epsPerAllele)^2 + epsPerAllele^2)
		} else {
			return((1-epsPerAllele)*epsPerAllele*2)
		}
	}
}

#' log-likelihood of relationship of multiple loci with genotyping error
#' 
multWithError <- function(kid, gMa, gPa, altFreqList, refFreqList, epsList){
	llh <- 0
	for(i in 1:length(kid)){
		potGpa <- gPa[i]
		potGma <- gMa[i]
		potKid <- kid[i]
		epsPerAllele <- epsList[i]
		altFreq <- altFreqList[i]
		refFreq <- refFreqList[i]
		
		llh <- llh + probGPerror(potKid, potGma, potGpa, altFreq, refFreq, epsPerAllele)
		
	}
	return(llh)
}


probUerror <- function(kid, gMa, gPa, altFreqList, refFreqList, epsList){
	llh <- 0
	for(i in 1:length(kid)){
		potGpa <- gPa[i]
		potGma <- gMa[i]
		potKid <- kid[i]
		epsPerAllele <- epsList[i]
		altFreq <- altFreqList[i]
		refFreq <- refFreqList[i]
		
		lh <- 0
		#cylcling over all possible combinations of true genotypes
		for(k in 0:2){
			for(m in 0:2){
				for(p in 0:2){
					
					lh <- lh + pmfBetaBinomial(m, 2, altFreq, refFreq) * pmfBetaBinomial(p, 2, altFreq, refFreq) *
						pmfBetaBinomial(k, 2, altFreq, refFreq) * p_obsG(potKid, k, epsPerAllele) * 
						p_obsG(potGma, m, epsPerAllele) * p_obsG(potGpa, p, epsPerAllele)

				}
			}
		}
		
		llh <- llh + log(lh)
	}

	return(llh)
}
