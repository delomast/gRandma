#' Calculate log-likelihood ratios to infer relationships
#' 
#' some more here
#' 
#' and some more here
#' 
#' @param gmaData
#' @param relationship
#' @param crossRecords dataframe or matrix: column 1 is population, column 2 is individual
#'  1 identifier, column 3 is individual 2 identifier
#' 
#' 
#' 
#' @export
inferGrandma <- function(gmaData, relationship = c("ssGP", "test1", "test2"), crossRecords = NULL){
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
	
	if(useUnsamp) {
		names(gmaData$unsampledPopsParams) <- as.character(popKey[match(names(gmaData$unsampledPopsParams), popKey[,1]),2])
	}
	# now we have to take into account cases where some (or all) baseline pops 
	# don't have "unsampledPopParams"
	new <- list()
	for(i in 1:length(allPops)){
		if(useUnsamp && as.character(i-1) %in% names(gmaData$unsampledPopsparams)){
			new[[i]] <- gmaData$unsampledPopsparams[[as.character(i-1)]]
		} else {
			new[[i]] <- gmaData$baselineParams[[i]]
		}
	}
	gmaData$unsampledPopParams <- new

	indivKey <- data.frame(indivName = c(gmaData$baseline[,2], gmaData$mixture[,1]), 
								  indivInt = 0:(nrow(gmaData$baseline) + nrow(gmaData$mixture) - 1), stringsAsFactors = FALSE)
	gmaData$baseline[,2] <- as.numeric(indivKey[match(gmaData$baseline[,2], indivKey[,1]),2])
	gmaData$mixture[,1] <- as.numeric(indivKey[match(gmaData$mixture[,1], indivKey[,1]),2])
	
	# make sure loci are in order in baseline, mixture, genotypeKey, baselineParams, and unsampledPopParams
	############# todo
	
	if(rel == "ssGP"){
		if(is.null(crossRecords)){
			crossRecords <- matrix(0, 0, 0)
		} else{
			crossRecords <- as.matrix(crossRecords)
		}
		# change things to numeric matrices when pass them to c++
		gmaData$genotypeKey <- lapply(gmaData$genotypeKey, as.matrix)
		results <- ssGP(as.matrix(gmaData$baseline), as.matrix(gmaData$mixture), crossRecords, 
							 gmaData$baselineParams, gmaData$unsampledPopParams, gmaData$genotypeKey)
		return(results) #testing
	}
	
	# decode everything from ints back into strings
	
}
