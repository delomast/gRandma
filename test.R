# CHNK data from gRandma demo
load("CHNK_for_testing.rda")

gmaInputFall <- createGmaInput(baseline = filterFall$baseline, mixture = filterFall$mixture,
										 unsampledPops = NULL,perAlleleError = .005,
										 dropoutProb = 0)

gmaInputSawt <- createGmaInput(baseline = filterSawt$baseline, mixture = filterSawt$mixture,
										 unsampledPops = NULL, perAlleleError = .005,
										 dropoutProb = .005)

fN <- falseGrandma(gmaInputSawt, relationship = c("sP"), 
								 llrToTest = c(1,5,10,15,20), N = 1000, seed = 7, itersPerMI = NULL, 
								 errorType = c("falseNegative"))
fN
fPunrel <- falseGrandma(gmaInputSawt, relationship = c("sP"), 
								 llrToTest = c(1,5,10,15,20), seed = 7, 
						 itersPerMI = c(rep(200, 15), rep(0, 238 - 15)), 
								 errorType = c("Unrel"))
fPunrel[[1]]
fPunrel[[2]][1:20,]

assign <- inferGrandma(gmaInputSawt, relationship = c("sP"), crossRecords = NULL, minLLR = 0,
								 filterLLR = TRUE, MIexcludeProb = .0001)

gmaInputSawt$baseline[1:200,2]
crossRec <- data.frame("OtsSAWT14S", gmaInputSawt$baseline[1:100,2], gmaInputSawt$baseline[101:200,2])
assignGP2 <- inferGrandma(gmaInputSawt, relationship = c("ssGP"), crossRecords = crossRec, minLLR = 0,
								 filterLLR = TRUE, MIexcludeProb = .0001)
identical(assignGP, assignGP2)
all.equal(assignGP, assignGP2)
head(assignGP)
head(assignGP2)

