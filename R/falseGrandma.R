#' estimate per-comparison error rates
#' 
#' some more here
#' 
#' and more here
#' 
#' @param gmaData the gmaData object containing your baseline populations and potential descendents. This 
#'   input is created by \code{createGmaInput}
#' @param relationship the relationship you want to test for: "ssGP" - single sided grandparentage 
#'   (a pair of either two maternal grandparents OR two paternal grandparents)
#' @param llrToTest a vector of llr thresholds to estimate false positive and false negative rates for
#' @param N the number of Monte Carlo samples to take to estimate the error rates
#' @param seed a positive integer to use as a seed for random number generation. If \code{NULL}, a seed is 
#'   chosen based on the current system time 
#' @param pairwise TRUE to estimate false positive rate for pairwise combinations of the baseline populations.
#'   In otherwords, the rate at which mixture individuals from one population have llr's equal to or above 
#'   the threshold when being tested against the baseline of another population.
#' 
#' @export
falseGrandma <- function(gmaData, relationship = c("ssGP", "sP"), 
								 llrToTest, N = 10000, seed = NULL, itersPerMI = NULL, 
								 errorType = c("falseNegative", "pairwise", "Unrel", "Aunt", "HalfAunt", "ParCous"),
								 MIexcludeProb = .0001, maxMissingGenos = round(0.1 * (length(itersPerMI) - 1), 0),
								 method = c("strat", "IS")){
	method <- match.arg(method)
	rel <- match.arg(relationship)
	tRel <- match.arg(errorType)
	if(is.null(seed)) seed <- ceiling(as.numeric(format(Sys.time(), "%S")) * 
												 	as.numeric(format(Sys.time(), "%j")) * 
												 	as.numeric(format(Sys.time(), "%M")))
	if(seed < 0) {
		warning("seed must not be negative, setting seed to 0")
		seed <- 0
	}
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
	
	if(rel == "ssGP"){
		if(tRel == "pairwise"){
			errResults <- otherPopERRORssGP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
							gmaData$missingParams, gmaData$genotypeKey,
							 gmaData$genotypeErrorRates, llrToTest, round(N), round(seed), skipBaseline)
		} else {
			if(method == "IS"){
				errResults <- list(ERRORssGP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
												gmaData$missingParams, gmaData$genotypeKey,
												gmaData$genotypeErrorRates, llrToTest, round(N), round(seed))
					)
			} else {
				# method is strat
				if(is.null(itersPerMI)) stop("itersPerMI must be input for this option.")
				if(any((itersPerMI %% 1) != 0)) stop("all itersPerMI must be integers")
				if(any(itersPerMI == 1)) warning("some itersPerMI are 1, SD will be undefined.")
				if(any(itersPerMI < 0)) stop("some itersPerMI are negative.")
				
				errResults <- strat_ERRORssGP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
														gmaData$missingParams, gmaData$genotypeKey,
														gmaData$genotypeErrorRates, llrToTest,
														itersPerMI,
														round(seed), c(0,1)[1], MIexcludeProb, maxMissingGenos)
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
	                          round(seed), skipBaseline, MIexcludeProb)

		} else if(tRel == "falseNegative"){
			errResults <- list(falseNeg_ERRORsP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
						gmaData$missingParams, gmaData$genotypeKey,
		         gmaData$genotypeErrorRates, llrToTest, round(N), round(seed),
						MIexcludeProb)
				)

		} else if (tRel %in% c("Unrel", "Aunt", "HalfAunt", "ParCous")) {
			errResults <- strat_ERRORsP(gmaData$baselineParams,
									gmaData$unsampledPopsParams, gmaData$missingParams,
									gmaData$genotypeKey,
									gmaData$genotypeErrorRates, llrToTest,
									itersPerMI,
									round(seed), c(0,1,2,3)[which(c("Unrel", "Aunt", "HalfAunt", "ParCous") == tRel)],
									MIexcludeProb)
		} else {
			stop("relationship and error type combination not recognized")
		}

	} else {
		stop("This relationship not set up at this time")
	}
	
	# turn pop names back to strings
	if(!is.data.frame(errResults)){
		for(i in 1:length(errResults)){
			if(tRel == "pairwise") errResults[[i]][,2] <- popKey[match(errResults[[i]][,2], popKey[,2]),1]
			errResults[[i]][,1] <- popKey[match(errResults[[i]][,1], popKey[,2]),1]
		}
	} else {
		if(tRel == "pairwise") errResults[,2] <- popKey[match(errResults[,2], popKey[,2]),1]
		errResults[,1] <- popKey[match(errResults[,1], popKey[,2]),1]
	}
	
	return(errResults)
}
