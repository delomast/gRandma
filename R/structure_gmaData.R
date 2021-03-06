#' implementing a structure to hold data
#' 
#' baseline: the dataframe of the baseline to use for inferring relationships, with genotypes
#'   recoded to integers
#' mixture: the dataframe of the mixture to use for inferring relationships, with genotypes
#'   recoded to integers
#' baselineParams: A list (on entry per baseline population) of lists of vectors 
#'   contining alphas for Dirichlet distributions for allele frequencies
#' unsamParams:  A list (on entry per baseline population) of lists of vectors 
#'   contining alphas for Dirichlet distributions for allele frequencies in the unsampled breeding 
#'   population where the corresponding baseline would be found
#' lociErr: a list of matrices (one for each locus). Each matrix gives the probability
#'   if the true genotype is X (rows), the observed genotype is Y. So, rows should always sum to 1.
#'   entries are named for the locus, rows and columns are named for the genotype. Row and column 
#'   names are just added to make human reading/editing easier. the order is what's used by the 
#'   internal functions. 
#' @keywords internal
#' @noRd

construct_grandma <- function(x){

	# checking required values
	gmaDataNames <- c("baseline", "mixture", "unsampledPops", "genotypeErrorRates", "genotypeKeys", 
							"alleleKeys", "baselineParams", "unsampledPopsParams", "missingParams")
	if(sum(names(x) != gmaDataNames) > 0) stop("Wrong names or order to make a gmaData object")
	if(!is.data.frame(x$baseline)) stop("baseline must be a dataframe")
	if(!is.null(x$mixture) && !is.data.frame(x$mixture)) stop("mixture must be a dataframe or NULL")
	if(!is.null(x$unsampledPops) && !is.data.frame(x$unsampledPops)) stop("unsampledPops must be a dataframe or NULL")
	if(!is.list(x$genotypeErrorRates)) stop("genotypeErrorRates must be a list")
	if(!is.list(x$genotypeKeys)) stop("genotypeKeys must be a list")
	if(!is.list(x$alleleKeys)) stop("alleleKeys must be a list")
	if(!is.list(x$baselineParams)) stop("baselineParams must be a list")
	if(!is.null(x$unsampledPopsParams) && !is.list(x$unsampledPopsParams)) stop("unsampledPopsParams must be a list or NULL")
	if(!is.list(x$missingParams)) stop("missingParams must be a list")
	
	class(x) <- "gmaData"
	return(x)
}

#' print method for gmaData
#' @param x a gmaData object
#' @param ... ignored
#' @export
print.gmaData <- function(x, ...){
	cat("\ngmaData object with\n\t")
	cat(length(x$baselineParams), " baseline population(s)\n\t")
	if(is.null(x$unsampledPopsParams)){
		cat("No unsampledPops\n\t")
	} else{
		cat(length(x$unsampledPopsParams), " unsampled population(s)\n\t")
	}
	cat(nrow(x$mixture), " mixture individuals\n\t")
	cat(length(x$genotypeErrorRates), " loci\n")
	cat("The number of loci by number of unique alleles is:\n")
	print(table(sapply(x$genotypeKeys, function(y){
		return(length(unique(y[,2])))
	})))
	cat("\n\n")
	
}

#' An example dataset for input
#' 
#' A dataset containing microhaps and SNPs for a large 
#' number of individuals across several populations
#' @format a dataframe
"data_mh_snp"

