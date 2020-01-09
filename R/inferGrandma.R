#' Calculate log-likelihood ratios to infer relationships
#' 
#' some more here
#' 
#' and some more here
#' 
#' @param gmaData the gmaData object containing your baseline populations and potential descendents. This 
#'   input is created by \code{createGmaInput}
#' @param relationship the relationship you want to test for: "ssGP" - single sided grandparentage 
#'   (a pair of either two maternal grandparents OR two paternal grandparents)
#' @param crossRecords dataframe or matrix: column 1 is population, column 2 is individual
#'  1 identifier, column 3 is individual 2 identifier
#' @param minLLR the minimum LLR to include in the results. Only used if \code{filterLLR} is \code{TRUE}.
#' @param filterLLR \code{TRUE} to filter results based on minLLR, \code{FALSE} to not.
#' @param MIexcludeProb the maximum probability of exclusion for a true grandparent pair due to 
#'   Mendelian incompatibilities. If \code{0}, then no filtering of potential grandparent pairs is performed
#'   based on Mendelian incompatibilities.
#' @export
inferGrandma <- function(gmaData, relationship = c("ssGP", "test1", "test2"), crossRecords = NULL, minLLR = 0,
								 filterLLR = TRUE, MIexcludeProb = .0001){
	rel <- match.arg(relationship)
	
	useUnsamp <- FALSE
	if(!is.null(gmaData$unsampledPopsParams)) useUnsamp <- TRUE
	
	# need to calc llr for all individuals in the mixture compared to all individuals in the baseline
	# this is really just a wrapper for subfunctions for each relationship
	
	# turn pop names and indiv names into ints for speed and avoid headache of dealing with strings
	# 0 index for c++
	allPops <- unique(gmaData$baseline[,1])
	popKey <- data.frame(popName = allPops, popInt = 0:(length(allPops)-1), stringsAsFactors = FALSE)
	gmaData$baseline[,1] <- as.numeric(popKey[match(gmaData$baseline[,1], popKey[,1]),2])
	names(gmaData$baselineParams) <- as.character(popKey[match(names(gmaData$baselineParams), popKey[,1]),2])
	gmaData$baselineParams <- gmaData$baselineParams[order(as.numeric(names(gmaData$baselineParams)))]
	
	cat("\nBaseline populations will be evaluated in the following order:\n")
	for(i in 1:nrow(popKey)){
		cat("\t", i, popKey[i,1], "\n")
	}
	cat("\n")
	
	if(useUnsamp) {
		names(gmaData$unsampledPopsParams) <- as.character(popKey[match(names(gmaData$unsampledPopsParams), popKey[,1]),2])
	}
	# now we have to take into account cases where some (or all) baseline pops 
	# don't have "unsampledPopParams"
	new <- list()
	for(i in 1:length(allPops)){
		if(useUnsamp && as.character(i-1) %in% names(gmaData$unsampledPopsParams)){
			new[[i]] <- gmaData$unsampledPopsParams[[as.character(i-1)]]
		} else {
			new[[i]] <- gmaData$baselineParams[[i]]
		}
	}
	gmaData$unsampledPopsParams <- new
	if(!is.null(crossRecords)) crossRecords[,1] <- as.numeric(popKey[match(crossRecords[,1], popKey[,1]),2])

	indivKey <- data.frame(indivName = c(gmaData$baseline[,2], gmaData$mixture[,1]), 
								  indivInt = 0:(nrow(gmaData$baseline) + nrow(gmaData$mixture) - 1), stringsAsFactors = FALSE)
	gmaData$baseline[,2] <- as.numeric(indivKey[match(gmaData$baseline[,2], indivKey[,1]),2])
	gmaData$mixture[,1] <- as.numeric(indivKey[match(gmaData$mixture[,1], indivKey[,1]),2])
	if(!is.null(crossRecords)) {
		crossRecords[,2] <- as.numeric(indivKey[match(crossRecords[,2], indivKey[,1]),2])
		crossRecords[,3] <- as.numeric(indivKey[match(crossRecords[,3], indivKey[,1]),2])
	}

	# make sure loci are in order in baseline, mixture, genotypeKey, genotypeErrorRates,
	#   baselineParams, and unsampledPopParams
	############# todo
	
	# turn NA genotypes into -9 for c++ to easily recognize
	gmaData$baseline[,3:ncol(gmaData$baseline)] <- apply(gmaData$baseline[,3:ncol(gmaData$baseline), drop = FALSE], 2, convertMissing)
	gmaData$mixture[,2:ncol(gmaData$mixture)] <- apply(gmaData$mixture[,2:ncol(gmaData$mixture), drop = FALSE], 2, convertMissing)
	
	if(rel == "ssGP"){
		if(is.null(crossRecords)){
			crossRecords <- matrix(0, 0, 0)
		} else{
			crossRecords <- as.matrix(crossRecords)
		}
		# change things to numeric matrices when pass them to c++
		gmaData$genotypeKey <- lapply(gmaData$genotypeKey, as.matrix)
		
		results <- ssGP(as.matrix(gmaData$baseline), as.matrix(gmaData$mixture), crossRecords, 
							 gmaData$baselineParams, gmaData$unsampledPopsParams, gmaData$genotypeKey,
							 gmaData$genotypeErrorRates, minLLR, MIexcludeProb, filterLLR)

	}
	
	# decode everything from ints back into strings
	results[,2] <- popKey[match(results[,2], popKey[,2]),1]
	results[,1] <- indivKey[match(results[,1], indivKey[,2]),1]
	results[,3] <- indivKey[match(results[,3], indivKey[,2]),1]
	results[,4] <- indivKey[match(results[,4], indivKey[,2]),1]
		
		
	return(results) #testing
	
}
