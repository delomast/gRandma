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
#' 
#' @export
falseGrandma <- function(gmaData, relationship = c("ssGP", "test1", "test2"), 
								 llrToTest, N = 10000, seed = NULL){
	rel <- match.arg(relationship)
	if(is.null(seed)) seed <- ceiling(as.numeric(format(Sys.time(), "%S")) * 
												 	as.numeric(format(Sys.time(), "%j")) * 
												 	as.numeric(format(Sys.time(), "%M")))
	
	useUnsamp <- FALSE
	if(!is.null(gmaData$unsampledPopsParams)) useUnsamp <- TRUE
	
	# turn pop names into ints for speed and avoid headache of dealing with strings
	# 0 index for c++
	allPops <- unique(gmaData$baseline[,1])
	popKey <- data.frame(popName = allPops, popInt = 0:(length(allPops)-1), stringsAsFactors = FALSE)
	names(gmaData$baselineParams) <- as.character(popKey[match(names(gmaData$baselineParams), popKey[,1]),2])
	gmaData$baselineParams <- gmaData$baselineParams[order(as.numeric(names(gmaData$baselineParams)))]
	
	cat("\nBaseline populations will be evaluated in the following order:\n")
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
	
	if(rel == "ssGP"){

		# change things to numeric matrices when pass them to c++
		gmaData$genotypeKey <- lapply(gmaData$genotypeKey, as.matrix)
		
		errResults <- ERRORssGP(gmaData$baselineParams, gmaData$unsampledPopsParams, 
										gmaData$missingParams, gmaData$genotypeKey,
                           gmaData$genotypeErrorRates, llrToTest, round(N), round(seed))

	}
	
	# turn pop names back to strings
	errResults[,1] <- popKey[match(errResults[,1], popKey[,2]),1]
	
	return(errResults)
}
