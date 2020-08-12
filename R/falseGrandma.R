#' estimate per-comparison error rates
#' 
#' some more here
#' 
#' and more here
#' 
#' @param gmaData the gmaData object containing your baseline populations and potential descendents. This 
#'   input is created by \code{createGmaInput}
#' @param relationship the relationship you want to test for: "ssGP" - single sided grandparentage 
#'   (a pair of either two maternal grandparents OR two paternal grandparents); "sP" - single parent inference 
#' @param llrToTest a vector of llr thresholds to estimate error rates for
#' @param N the number of Monte Carlo samples to take to estimate the error rates (ignored for stratified methods)
#' @param seed a positive integer to use as a seed for random number generation.
#' @param itersperMI the number of iterations per Mendelian incompatibility, in order of 0, 1, ... 
#'   (ignored for non-stratified methods)
#' @param errorType the type of error estimate to make
#' @param MIexcludeProb the maximum probability of exclusion for a true relationship due to 
#'   Mendelian incompatibilities. If \code{0}, then no filtering is performed
#'   based on Mendelian incompatibilities.
#' @param maxMissingGenos the maximum number of missing genotypes a sample can have before you would 
#'   choose to omit it from analysis
#' @param method strat for stratified, IS for importance sampling. Only used for ssGP. Do not use method = "old", 
#'   this is currently for internal testing and will be soon removed. 
#' 
#' @export
falseGrandma <- function(gmaData, relationship = c("ssGP", "sP"), 
								 llrToTest, N = 10000, seed = sample.int(.Machine$integer.max, 1), itersPerMI = NULL, 
								 errorType = c("falseNegative", "pairwise", "Unrel", "Aunt", "HalfAunt", "ParCous",
								 				  "True_GAunt", "True_Unrel", "True_HGAunt", "True_GpCous", 
								 				  "GAunt_Unrel", "HGAunt_Unrel", "GpCous_Unrel", "GAunt", "GAunt_HGAunt", 
								 				  "Gaunt_GpCous", "HGAunt", "HGAunt_GpCous", "GpCous"),
								 MIexcludeProb = .0001, maxMissingGenos = NULL,
								 method = c("strat", "IS", "old")){
	method <- match.arg(method)
	rel <- match.arg(relationship)
	tRel <- match.arg(errorType)
	if(seed < 0) {
		warning("seed must not be negative, setting seed to 1")
		seed <- 1
	}
	if(is.null(maxMissingGenos)) maxMissingGenos <- ceiling(.1 * length(gmaData$genotypeKey))
	if(maxMissingGenos %% 1 != 0) stop("maxMissingGenos must be an integer")
	if(method == "IS" && (tRel != "Unrel" || rel != "ssGP")) stop("method of IS is only an option for ssGP and Unrel")
	
	useUnsamp <- FALSE
	skipBaseline <- c() # so unsampled pops aren't tested as a baseline
	if(!is.null(gmaData$unsampledPopsParams)){
		useUnsamp <- TRUE
		if(tRel == "pairwise"){
			# need to perform pairwise with unsampled pops
			addedUnsampledPopsNames <- paste0("UnsampledPop_", names(gmaData$unsampledPopsParams))
			# just in case there is already a pop with that name
			index <- 1
			while(any(addedUnsampledPopsNames %in% names(gmaData$baselineParams))){
				addedUnsampledPopsNames[addedUnsampledPopsNames %in% names(gmaData$baselineParams)] <-
					paste0(addedUnsampledPopsNames[addedUnsampledPopsNames %in% names(gmaData$baselineParams)],
							 "_added_", index)
				index <- index + 1
			}
			rm(index)
			addedUPparams <- gmaData$unsampledPopsParams
			names(addedUPparams) <- addedUnsampledPopsNames
			gmaData$baselineParams <- c(gmaData$baselineParams, addedUPparams)
			skipBaseline <- addedUnsampledPopsNames
		}
	}
	
	# turn pop names into ints for speed and avoid headache of dealing with strings
	# 0 index for c++
	allPops <- names(gmaData$baselineParams)
	popKey <- data.frame(popName = allPops, popInt = 0:(length(allPops)-1), stringsAsFactors = FALSE)
	names(gmaData$baselineParams) <- as.character(popKey[match(names(gmaData$baselineParams), popKey[,1]),2])
	gmaData$baselineParams <- gmaData$baselineParams[order(as.numeric(names(gmaData$baselineParams)))]
	if(length(skipBaseline) > 0) skipBaseline <- popKey[match(skipBaseline, popKey[,1]),2]
	skipBaseline <- as.numeric(skipBaseline)
	
	cat("\nPopulations will be evaluated in the following order:\n")
	for(i in 1:nrow(popKey)){
		cat("\t", i, popKey[i,1], "\n")
	}
	cat("\n")
	
	if(useUnsamp) names(gmaData$unsampledPopsParams) <- as.character(popKey[match(names(gmaData$unsampledPopsParams), popKey[,1]),2])
	
	# now we have to take into account cases where some (or all) baseline pops don't have "unsampledPopParams"
	new <- list()
	for(i in 1:length(allPops)){
		if(useUnsamp && as.character(i-1) %in% names(gmaData$unsampledPopsParams)){
			new[[i]] <- gmaData$unsampledPopsParams[[as.character(i-1)]]
		} else {
			new[[i]] <- gmaData$baselineParams[[i]]
		}
	}
	gmaData$unsampledPopsParams <- new
	
	# change things to numeric matrices when pass them to c++
	gmaData$genotypeKey <- lapply(gmaData$genotypeKey, as.matrix)
	
	ssGP_err_rels <- c("Unrel", "True_GAunt", "True_Unrel", "True_HGAunt", "True_GpCous", 
							 "GAunt_Unrel", "HGAunt_Unrel", "GpCous_Unrel", "GAunt", "GAunt_HGAunt", 
							 "Gaunt_GpCous", "HGAunt", "HGAunt_GpCous", "GpCous")
	
	if(rel == "ssGP"){
		if(tRel %in% c(ssGP_err_rels, "pairwise") && method == "strat"){
			if(is.null(itersPerMI)) stop("itersPerMI must be input for this option.")
			if(any((itersPerMI %% 1) != 0)) stop("all itersPerMI must be integers")
			if(any(itersPerMI == 1)) warning("some itersPerMI are 1, SD will be undefined.")
			if(any(itersPerMI < 0)) stop("some itersPerMI are negative.")
		}
		if(tRel == "pairwise"){
			errResults <- strat_otherPopERRORssGP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
											gmaData$missingParams, gmaData$genotypeKey,
											gmaData$genotypeErrorRates, llrToTest, itersPerMI,
											round(seed), skipBaseline, MIexcludeProb,
											maxMissingGenos)
		} else {
			if(method == "IS"){
				errResults <- list(ERRORssGP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
												gmaData$missingParams, gmaData$genotypeKey,
												gmaData$genotypeErrorRates, llrToTest, round(N), round(seed),
												MIexcludeProb, maxMissingGenos)
					)
			} else if (tRel == "falseNegative"){
				# just running the IS function, doesn't really add significant comp time
				errResults <- ERRORssGP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
													  gmaData$missingParams, gmaData$genotypeKey,
													  gmaData$genotypeErrorRates, llrToTest, round(N), round(seed),
													  MIexcludeProb, maxMissingGenos)
				errResults <- list(errResults[,c(1,2,5,6)])
				
			} else if (tRel %in% ssGP_err_rels) {
				# method is strat, false positive of some sort
				errResults <- strat_ERRORssGP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
														gmaData$missingParams, gmaData$genotypeKey,
														gmaData$genotypeErrorRates, llrToTest,
														itersPerMI,
														round(seed), c(0:length(ssGP_err_rels))[which(ssGP_err_rels == tRel)], 
														MIexcludeProb, maxMissingGenos)
			} else {
				stop("relationship and error type combination not recognized")
			}
		}
		
	} else if(rel == "sP"){
		if(tRel %in% c("Unrel", "Aunt", "HalfAunt", "ParCous", "pairwise")){
			if(is.null(itersPerMI)) stop("itersPerMI must be input for this option.")
			if(any((itersPerMI %% 1) != 0)) stop("all itersPerMI must be integers")
			if(any(itersPerMI == 1)) warning("some itersPerMI are 1, SD will be undefined.")
			if(any(itersPerMI < 0)) stop("some itersPerMI are negative.")
		}
		if(tRel == "pairwise"){

			errResults <- strat_otherPopERRORsP(gmaData$baselineParams,
	                          gmaData$unsampledPopsParams, gmaData$missingParams,
	                          gmaData$genotypeKey,
	                          gmaData$genotypeErrorRates, llrToTest,
	                          itersPerMI,
	                          round(seed), skipBaseline, MIexcludeProb, maxMissingGenos)

		} else if(tRel == "falseNegative"){
			errResults <- list(falseNeg_ERRORsP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
						gmaData$missingParams, gmaData$genotypeKey,
		         gmaData$genotypeErrorRates, llrToTest, round(N), round(seed),
						MIexcludeProb, maxMissingGenos)
				)

		} else if (tRel %in% c("Unrel", "Aunt", "HalfAunt", "ParCous")) {
			
			if(method == "old"){
				# for testing
				errResults <- old_strat_ERRORsP(gmaData$baselineParams,
													 gmaData$unsampledPopsParams, gmaData$missingParams,
													 gmaData$genotypeKey,
													 gmaData$genotypeErrorRates, llrToTest,
													 itersPerMI,
													 round(seed), c(0,1,2,3)[which(c("Unrel", "Aunt", "HalfAunt", "ParCous") == tRel)],
													 MIexcludeProb)
			} else {
				errResults <- strat_ERRORsP(gmaData$baselineParams,
									gmaData$unsampledPopsParams, gmaData$missingParams,
									gmaData$genotypeKey,
									gmaData$genotypeErrorRates, llrToTest,
									itersPerMI,
									round(seed), c(0,1,2,3)[which(c("Unrel", "Aunt", "HalfAunt", "ParCous") == tRel)],
									MIexcludeProb, maxMissingGenos)
			}
		} else {
			stop("relationship and error type combination not recognized")
		}

	} else {
		stop("This relationship not set up at this time")
	}
	
	# turn pop names back to strings
	for(i in 1:length(errResults)){
		if(tRel == "pairwise") errResults[[i]][,2] <- popKey[match(errResults[[i]][,2], popKey[,2]),1]
		errResults[[i]][,1] <- popKey[match(errResults[[i]][,1], popKey[,2]),1]
	}
	
	return(errResults)
}
