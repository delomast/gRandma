# CHNK data from gRandma demo
load("CHNK_for_testing.rda")

gmaInputFall <- createGmaInput(baseline = filterFall$baseline, mixture = filterFall$mixture,
										 unsampledPops = NULL,perAlleleError = .005,
										 dropoutProb = 0)

gmaInputSawt <- createGmaInput(baseline = filterSawt$baseline, mixture = filterSawt$mixture,
										 unsampledPops = NULL, perAlleleError = .005,
										 dropoutProb = .005)
commonCol <- intersect(colnames(filterSawt$baseline), colnames(filterFall$baseline))
bothBase <- rbind(filterSawt$baseline[,commonCol], filterFall$baseline[,commonCol])
gmaInputBoth <- createGmaInput(baseline = bothBase,
										 unsampledPops = NULL, perAlleleError = .005,
										 dropoutProb = .005)

fN <- falseGrandma(gmaInputSawt, relationship = c("sP"), 
								 llrToTest = c(1,5,10,15,20), N = 1000, seed = 7, itersPerMI = NULL, 
								 errorType = c("falseNegative"))
fN
fPunrel <- falseGrandma(gmaInputSawt, relationship = c("sP"), 
								 llrToTest = c(1,5,10,15,20), seed = 7, 
						 itersPerMI = c(rep(2000, 15), rep(0, 238 - 15)), 
								 errorType = c("Aunt"), method = "old")
system.time(
fPunrel2 <- falseGrandma(gmaInputSawt, relationship = c("sP"), 
								llrToTest = c(1,5,10,15,20), seed = 7, 
								itersPerMI = rep(10000, 15), 
								errorType = c("Aunt"))
)
fPunrel[[1]]
fPunrel2[[1]]
fPunrel[[2]][1:10,]
fPunrel2[[2]][1:10,]

fPpair <- falseGrandma(gmaInputBoth, relationship = c("sP"), 
								llrToTest = c(10), seed = 7, 
								itersPerMI = rep(2000, 15), 
								errorType = c("pairwise"))
fPpair[[1]]


falseGrandma(gmaInputSawt, relationship = c("sP"), 
				 llrToTest = c(1,5,10,15,20), seed = 7, 
				 itersPerMI = c(rep(2000, 15), rep(0, 238 - 15)), 
				 errorType = c("Aunt"))[[1]]

temp <- falseGrandma(gmaInputSawt, relationship = c("sP"), 
				 llrToTest = c(1,5,10,15,20), seed = 7, 
				 itersPerMI = c(rep(20, 15), rep(0, 238 - 15)), 
				 errorType = c("Unrel"))

temp <- falseGrandma(gmaInputBoth, relationship = c("sP"), 
				 llrToTest = c(1,5,10,15,20), seed = 7, 
				 itersPerMI = rep(20, 238), 
				 errorType = c("pairwise"))

assign <- inferGrandma(gmaInputSawt, relationship = c("sP"), crossRecords = NULL, minLLR = 0,
								 filterLLR = TRUE, MIexcludeProb = .0001)


gmaInputSawt$baseline[1:200,2]
crossRec <- data.frame("OtsSAWT14S", gmaInputSawt$baseline[1:100,2], gmaInputSawt$baseline[101:200,2])
assignGP <- inferGrandma(gmaInputSawt, relationship = c("ssGP"), crossRecords = crossRec, minLLR = 0,
								 filterLLR = TRUE, MIexcludeProb = .0001)
head(assignGP)

fpStrat <- falseGrandma(gmaInputSawt, relationship = c("ssGP"), 
								llrToTest = c(1,5,10,15,20), seed = 7, 
								itersPerMI = c(rep(20, 15), rep(0, 238 - 15)), 
								errorType = c("Unrel"))
colnames(filterSawt$baseline)[1:10]
gmaInputSawtSmall <- createGmaInput(baseline = filterSawt$baseline[,1:102],
										 unsampledPops = NULL, perAlleleError = .005,
										 dropoutProb = .005)


system.time(
strat <- falseGrandma(gmaInputSawt, relationship = c("ssGP"), 
				 llrToTest = seq(1,10,1), seed = 7, 
				 itersPerMI = rep(10000, 10), 
				 errorType = c("Unrel"))
)

is <- falseGrandma(gmaInputSawt, relationship = c("ssGP"), 
				 llrToTest = seq(1,10,1), seed = 7, 
				 N = 6000, method = "IS", 
				 errorType = c("Unrel"))

is

strat[[1]]
plot(is[[1]]$falsePos, strat[[1]]$falsePosUnrel)
abline(0,1)



plot(is[[1]]$falsePos, strat[[1]]$falsePosUnrel)
abline(0,1)


for(r in c("Unrel", "True_GAunt", "True_Unrel", "True_HGAunt", "True_GpCous", 
			  "GAunt_Unrel", "HGAunt_Unrel", "GpCous_Unrel", "GAunt", "GAunt_HGAunt", 
			  "Gaunt_GpCous", "HGAunt", "HGAunt_GpCous", "GpCous")){
	print(
		falseGrandma(gmaInputSawt, relationship = c("ssGP"), 
						 llrToTest = c(1), seed = 7, 
						 itersPerMI = rep(100, 10), 
						 errorType = r)
	)
}

strat_postfunk <- falseGrandma(gmaInputSawt, relationship = c("ssGP"), 
							 llrToTest = seq(1,10,1), seed = 7, 
							 itersPerMI = rep(1000, 10), 
							 errorType = c("Unrel"))

identical(strat_prefunk, strat_postfunk)

pair_ssGP <- falseGrandma(gmaInputFall, relationship = c("ssGP"), 
				 llrToTest = seq(1,10,1), seed = 7, 
				 itersPerMI = rep(1000, 10), 
				 errorType = c("pairwise"))

pair_ssGP[[1]]
pair_ssGP[[2]]
