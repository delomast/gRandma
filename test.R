


pmfBetaBinomial(1,2,1000,1000)

potGma <- c(0,0,0,0)
potGpa <- c(2,2,2,2)
potGma <- 0
potGpa <- 2
potKid <- 2

refFreq <- 100
altFreq <- 1


probU <- pmfBetaBinomial(potGma, 2, altFreq, refFreq) * pmfBetaBinomial(potGpa, 2, altFreq, refFreq) *
	pmfBetaBinomial(potKid, 2, altFreq, refFreq)
probGP <- pmfBetaBinomial(potGma, 2, altFreq, refFreq) * pmfBetaBinomial(potGpa, 2, altFreq, refFreq)
probSum <- 0
for(f in 0:2){
	for(m in 0:2){
		probSum <- probSum +
			pmfBetaBinomial(f, 2, altFreq, refFreq) * condProbTrio(potGma, potGpa, m) * condProbTrio(f, m, potKid)
	}
}
probGP <- probGP * probSum

probGP/probU

probGP(c(1,1,1,1),c(0,0,0,0), c(2,2,2,2), 1,10) - probU(c(1,1,1,1),c(0,0,0,0), c(2,2,2,2), 1,10)
probGP(c(1,1,1,1),c(0,0,0,0), c(2,2,2,2), 1,100) - probU(c(1,1,1,1),c(0,0,0,0), c(2,2,2,2), 1,100)


exp(probGP(c(1,1,1,1),c(0,0,0,0), c(2,2,2,2), 1,10))

probGPerror(1, 2, 2, 100, 1, .005)


multWithError(kid = c(1,1,1,1), gMa = c(0,0,0,0), gPa = c(2,2,2,2), altFreqList = rep(1,4), refFreqList = rep(100,4), epsList = rep (.005,4)) -
probUerror(kid = c(1,1,1,1), gMa = c(0,0,0,0), gPa = c(2,2,2,2), altFreqList = rep(1,4), refFreqList = rep(100,4), epsList = rep (.005,4))


#######################################################################################
### testing with known grandparents
testData <- read.csv("DWOR_CHNK.txt", sep = "\t", stringsAsFactors = FALSE)
head(testData)

genos <- testData[,15:ncol(testData)]
str(genos)
genos <- apply(genos,2, function(x){
	y <- rep(-9, length(x))
	y[x == "AA"] <- 0
	y[x == "AB"] <- 1
	y[x == "BA"] <- 1
	y[x == "BB"] <- 2
	# y[x == "00"] <- -9
	return(as.numeric(y))
	
})
rownames(genos) <- testData$Individual.Name

boolKeep <- apply(genos[testData$Pedigree == "OtsDWOR10S",],2, function(x) sum(x != -9)/length(x) > .7)
genos <- genos[,boolKeep]

boolKeep <- apply(genos[testData$Pedigree == "OtsDWOR18S",],2, function(x) sum(x != -9)/length(x) > .7)
genos <- genos[,boolKeep]

genos <- genos[,colnames(genos) != "Ots_SEXY3.1"]
head(genos)

genos[1:10,1:10]
apply(genos[testData$Pedigree == "OtsDWOR10S",],1,function(x) sum(x != -9)/length(x))
boolKeep <- apply(genos,1,function(x) sum(x != -9)/length(x) >= .9)
genos <- genos[boolKeep,]

## filter out nonvariable SNPs
boolKeep <- apply(genos, 2, function(x) sum(x %in% c(-9, 0)) != length(x) && sum(x %in% c(-9, 2)) != length(x) )
genos <- genos[,boolKeep]

#working with 2010 and 2018
#make list of all crosses

baseline <- testData[testData$Pedigree == "OtsDWOR10S",]
baseline <- baseline[baseline$Individual.Name %in% rownames(genos),]
baseline <- baseline[,1:14]
#check CrossedAdditional
sum(!is.na(baseline$CrossedAdditional))
#0
baseline <- baseline[,colnames(baseline) != "CrossedAdditional"]
crossCols <- colnames(baseline)[grepl("Cross", colnames(baseline))]
crosses <- matrix(0,0,2)
for(i in 1:nrow(baseline)){
	indiv <- baseline[i,"Individual.Name"]
	partners <- baseline[i,crossCols]
	partners <- partners[!is.na(partners) & partners != ""]
	if(length(partners) < 1) next
	tempCross <- cbind(indiv, partners)
	rem <- c()
	for(j in 1:nrow(tempCross)){
		if(sum(crosses[,1] == tempCross[j,2] & crosses[,2] == tempCross[j,1]) > 0) rem <- c(rem,j)
	}
	if(length(rem > 0)) tempCross <- tempCross[-rem,]
	crosses <- rbind(crosses, tempCross)
}
nrow(baseline)
nrow(crosses)

# now allele counts in baseline pop and add .5
## this is posterior mean with beta(.5,.5) as the prior
refFreqList <- apply(genos[baseline$Individual.Name,],2,function(x) sum(x == 1) + 2*sum(x == 0) + .5)
altFreqList <- apply(genos[baseline$Individual.Name,],2,function(x) sum(x == 1) + 2*sum(x == 2) + .5)

# define mixture
mixture <- testData$Individual.Name[testData$Pedigree == "OtsDWOR18S"]
mixture <- mixture[mixture %in% rownames(genos)]

#now grandparentage
# results <- matrix(nrow = length(mixture), ncol=nrow(crosses))
# for(k in 1:length(mixture)){
# 	kidGenos <- genos[mixture[k],]
# 	for(c in 1:nrow(crosses)){
# 		Gma <- genos[crosses[c,1],]
# 		Gpa <- genos[crosses[c,2],]
# 		
# 		results[k,c] <- multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
# 			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
# 	}
# }
#### not run, takes too long

### takes too long, selecting just some known grandaparents and some random non-grandparents

mixture <- testData[testData$Pedigree == "OtsDWOR18S",1:14]
mixture <- mixture[mixture$Individual.Name %in% rownames(genos),]
mixture <- mixture[!is.na(mixture$GenMa),]
mixture$trueGpa <- ""
mixture$trueGma <- ""
for(i in 1:nrow(mixture)){
	if(mixture$GenMa[i] %in% testData$Individual.Name){
		# tempGpa <- testData$GenPa[match(mixture$GenMa, testData$Individual.Name),]
		# tempGma <- testData$GenMa[match(mixture$GenMa, testData$Individual.Name),]
		mixture$trueGpa[i] <- testData$GenPa[match(mixture$GenMa[i], testData$Individual.Name)]
		mixture$trueGma[i] <- testData$GenMa[match(mixture$GenMa[i], testData$Individual.Name)]
	}
}

knownGPs <- mixture[mixture$trueGma %in% baseline$Individual.Name & mixture$trueGpa %in% baseline$Individual.Name,]

llr_knownGPs <- c()
for(i in sample(1:nrow(knownGPs), 50, replace = FALSE)){
	
	kidGenos <- genos[knownGPs[i,"Individual.Name"],]
	Gma <- genos[knownGPs[i,"trueGma"],]
	Gpa <- genos[knownGPs[i,"trueGpa"],]
		
	llr_knownGPs <- c(llr_knownGPs,
							multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
							)
	print(i)
}
print("switch relationship")
llr_notGPs <- c()
for(i in sample(1:nrow(knownGPs), 50, replace = FALSE)){
	
	kidGenos <- genos[knownGPs[i,"Individual.Name"],]
	randGpar <- sample(baseline$Individual.Name[!(baseline$Individual.Name %in% c(knownGPs[i,"trueGma"], knownGPs[i,"trueGpa"]))],2, replace = FALSE)
	Gma <- genos[randGpar[1],]
	Gpa <- genos[randGpar[2],]
	
	res <- multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
	if(res > 2) break
	llr_notGPs <- c(llr_notGPs,
							res
							)
	print(i)
}

boxplot(llr_knownGPs, llr_notGPs)
plot(c(llr_knownGPs, llr_notGPs), col = c(rep("blue", 50), rep("red", 50)), pch = 16, main = paste(ncol(genos), "SNPs"))
abline(h=0)
saved <- c(llr_knownGPs, llr_notGPs)

# run overnight

llr_knownGPs_all <- c()
for(i in 1:nrow(knownGPs)){
	
	kidGenos <- genos[knownGPs[i,"Individual.Name"],]
	Gma <- genos[knownGPs[i,"trueGma"],]
	Gpa <- genos[knownGPs[i,"trueGpa"],]
		
	llr_knownGPs_all <- c(llr_knownGPs_all,
							multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
							)
	print(i)
}
print("switch relationship")
llr_notGPs_lots <- c()
for(i in sample(1:nrow(knownGPs), 1000, replace = TRUE)){
	
	kidGenos <- genos[knownGPs[i,"Individual.Name"],]
	randGpar <- sample(baseline$Individual.Name[!(baseline$Individual.Name %in% c(knownGPs[i,"trueGma"], knownGPs[i,"trueGpa"]))],2, replace = FALSE)
	Gma <- genos[randGpar[1],]
	Gpa <- genos[randGpar[2],]
	
	res <- multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
	llr_notGPs_lots <- c(llr_notGPs_lots,
							res
							)
	print(i)
}

save(llr_knownGPs_all, llr_notGPs_lots, file = "overnight_saved.rda")

load("overnight_saved.rda")

boxplot(llr_knownGPs_all, llr_notGPs_lots)
plot(c(llr_knownGPs_all, llr_notGPs_lots), col = c(rep("blue", length(llr_knownGPs_all)), rep("red", length(llr_notGPs_lots))), 
	  pch = 16, main = paste(ncol(genos), "SNPs"))
abline(h=0)
library(tidyverse)

dataPlot <- data.frame(llr = c(llr_knownGPs_all, llr_notGPs_lots), type = c(rep("trueGP", length(llr_knownGPs_all)), rep("notGP", length(llr_notGPs_lots))))
fewSNPsPlot <- ggplot(dataPlot, aes(x=type, y=llr)) + geom_violin(draw_quantiles = c(.25, .5, .75), fill="#A4A4A4") + 
	theme(legend.position = "none") + geom_hline(yintercept = 0) + geom_point() + ggtitle(paste("DWOR CHNK", ncol(genos), "SNPs"))
fewSNPsPlot

###############################################
#### Now look with Sockeye - more SNPs


### testing with known grandparents
testData <- read.csv("all_sockeye_BCB.txt", sep = "\t", stringsAsFactors = FALSE)
head(testData)

genos <- testData[,28:ncol(testData)]
str(genos)
genos <- apply(genos,2, function(x){
	y <- rep(-9, length(x))
	y[x == "AA"] <- 0
	y[x == "AB"] <- 1
	y[x == "BA"] <- 1
	y[x == "BB"] <- 2
	# y[x == "00"] <- -9
	return(as.numeric(y))
	
})
rownames(genos) <- testData$Individual.Name

genos[1:10,1:10]
boolKeep <- apply(genos,1,function(x) sum(x != -9)/length(x) >= .9)
genos <- genos[boolKeep,]

## filter out nonvariable SNPs
boolKeep <- apply(genos, 2, function(x) sum(x %in% c(-9, 0)) != length(x) && sum(x %in% c(-9, 2)) != length(x) )
genos <- genos[,boolKeep]
str(genos)

## remove failed individuals from testData
testData <- testData[testData$Individual.Name %in% rownames(genos),]

#make list of all crosses
table(substr(rownames(genos), 1,10))
# OneEAGL08B OneEAGL10B OneEAGL11B OneEAGL12B OneEAGL13B OneEAGL14B 
#         82        580       1124        890       1365       1362 
# OneEAGL15B OneEAGL16B OneEAGL17B OneNOAA08B OneNOAA10B OneNOAA11B 
#       1254       1343       1383          5        361       1033 
# OneNOAA12B OneNOAA13B OneNOAA14B OneNOAA15B OneNOAA16B OneNOAA17B 
#       1357       1326       1337       1427       1383       1398
## lets work with mixture of SY2017
mixtureInds <- rownames(genos)[grepl("17BCB", rownames(genos))]

mixture <- testData[testData$Individual.Name %in% mixtureInds, 1:27]
mixture <- mixture[!is.na(mixture$GenMa),]
mixture$trueGpa <- ""
mixture$trueGma <- ""
for(i in 1:nrow(mixture)){
	if(mixture$GenMa[i] %in% testData$Individual.Name){
		# tempGpa <- testData$GenPa[match(mixture$GenMa, testData$Individual.Name),]
		# tempGma <- testData$GenMa[match(mixture$GenMa, testData$Individual.Name),]
		mixture$trueGpa[i] <- testData$GenPa[match(mixture$GenMa[i], testData$Individual.Name)]
		mixture$trueGma[i] <- testData$GenMa[match(mixture$GenMa[i], testData$Individual.Name)]
	}
}

baselinePeds <- unique(substr(c(mixture$trueGma, mixture$trueGpa), 1,12))
baselinePeds <- baselinePeds[baselinePeds != ""]

baseline <- testData[testData$Pedigree %in% baselinePeds,]
baseline <- baseline[baseline$Individual.Name %in% rownames(genos),]
baseline <- baseline[,1:27]
#check CrossedAdditional
sum(!is.na(baseline$CrossedAdditional))
#0
baseline <- baseline[,colnames(baseline) != "CrossedAdditional"]
crossCols <- colnames(baseline)[grepl("Cross", colnames(baseline))]
crosses <- matrix(0,0,2)
for(i in 1:nrow(baseline)){
	indiv <- baseline[i,"Individual.Name"]
	partners <- baseline[i,crossCols]
	partners <- partners[!is.na(partners) & partners != ""]
	if(length(partners) < 1) next
	tempCross <- cbind(indiv, partners)
	rem <- c()
	for(j in 1:nrow(tempCross)){
		if(sum(crosses[,1] == tempCross[j,2] & crosses[,2] == tempCross[j,1]) > 0) rem <- c(rem,j)
	}
	if(length(rem > 0)) tempCross <- tempCross[-rem,]
	crosses <- rbind(crosses, tempCross)
}
nrow(baseline)
nrow(crosses)

# now allele counts in baseline pop and add .5
## this is posterior mean with beta(.5,.5) as the prior
refFreqList <- apply(genos[baseline$Individual.Name,],2,function(x) sum(x == 1) + 2*sum(x == 0) + .5)
altFreqList <- apply(genos[baseline$Individual.Name,],2,function(x) sum(x == 1) + 2*sum(x == 2) + .5)


### selecting just some known grandaparents and some random non-grandparents

knownGPs <- mixture[mixture$trueGma %in% baseline$Individual.Name & mixture$trueGpa %in% baseline$Individual.Name,]

llr_knownGPs <- c()
for(i in sample(1:nrow(knownGPs), 50, replace = FALSE)){
	
	kidGenos <- genos[knownGPs[i,"Individual.Name"],]
	Gma <- genos[knownGPs[i,"trueGma"],]
	Gpa <- genos[knownGPs[i,"trueGpa"],]
		
	llr_knownGPs <- c(llr_knownGPs,
							multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
							)
	print(i)
}
print("switch relationship")
llr_notGPs <- c()
for(i in sample(1:nrow(knownGPs), 50, replace = FALSE)){
	
	kidGenos <- genos[knownGPs[i,"Individual.Name"],]
	randGpar <- sample(baseline$Individual.Name[!(baseline$Individual.Name %in% c(knownGPs[i,"trueGma"], knownGPs[i,"trueGpa"]))],2, replace = FALSE)
	Gma <- genos[randGpar[1],]
	Gpa <- genos[randGpar[2],]
	
	res <- multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
	if(res > 2) break
	llr_notGPs <- c(llr_notGPs,
							res
							)
	print(i)
}

boxplot(llr_knownGPs, llr_notGPs)
plot(c(llr_knownGPs, llr_notGPs), col = c(rep("blue", 50), rep("red", 50)), pch = 16, main = paste(ncol(genos), "SNPs"))
abline(h=0)
saved <- c(llr_knownGPs, llr_notGPs)

# runs for a while

llr_knownGPs_all_sock <- c()
for(i in 1:min(nrow(knownGPs), 400)){
	
	kidGenos <- genos[knownGPs[i,"Individual.Name"],]
	Gma <- genos[knownGPs[i,"trueGma"],]
	Gpa <- genos[knownGPs[i,"trueGpa"],]
		
	llr_knownGPs_all_sock <- c(llr_knownGPs_all_sock,
							multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
							)
	print(i)
}
print("switch relationship")
llr_notGPs_lots_sock <- c()
for(i in sample(1:nrow(knownGPs), 400, replace = TRUE)){
	
	kidGenos <- genos[knownGPs[i,"Individual.Name"],]
	randGpar <- sample(baseline$Individual.Name[!(baseline$Individual.Name %in% c(knownGPs[i,"trueGma"], knownGPs[i,"trueGpa"]))],2, replace = FALSE)
	Gma <- genos[randGpar[1],]
	Gpa <- genos[randGpar[2],]
	
	res <- multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
	llr_notGPs_lots_sock <- c(llr_notGPs_lots_sock,
							res
							)
	print(i)
}

save(llr_knownGPs_all_sock, llr_notGPs_lots_sock, file = "overnight_saved_sock.rda")

load("overnight_saved_sock.rda")

dataPlot2 <- data.frame(llr = c(llr_knownGPs_all_sock, llr_notGPs_lots_sock), 
								type = c(rep("trueGP", length(llr_knownGPs_all_sock)), rep("notGP", length(llr_notGPs_lots_sock))))
manySNPsPlot <- ggplot(dataPlot2, aes(x=type, y=llr)) + geom_violin(draw_quantiles = c(.25, .5, .75), fill="#A4A4A4") + 
	theme(legend.position = "none") + geom_hline(yintercept = 0) + geom_point() + ggtitle(paste("BCB SOCK", ncol(genos), "SNPs"))
manySNPsPlot

library(ggridges)

manyMarkersRidge <- ggplot(dataPlot2, aes(x = llr, y = type)) + 
		# geom_density_ridges(stat = "binline", bins = 20, draw_baseline = FALSE) + #geom_point() + 
	geom_density_ridges(aes(point_color = type), rel_min_height = .01, jittered_points = TRUE, scale = 3) +
	geom_vline(xintercept = max(dataPlot2$llr[dataPlot2$type == "notGP"]), colour = "red") + 
	geom_vline(xintercept = min(dataPlot2$llr[dataPlot2$type == "trueGP"]), colour = "blue") + 
	scale_discrete_manual(aesthetics = "point_colour", values =c("red", "blue")) + ggtitle(paste("BCB SOCK", ncol(genos), "SNPs"))

fewMarkersRidge <- ggplot(dataPlot, aes(x = llr, y = type)) + 
		# geom_density_ridges(stat = "binline", bins = 20, draw_baseline = FALSE) + #geom_point() + 
	geom_density_ridges(aes(point_color = type), rel_min_height = .01, jittered_points = TRUE, scale = 3) +
	geom_vline(xintercept = max(dataPlot$llr[dataPlot$type == "notGP"]), colour = "red") + 
	geom_vline(xintercept = min(dataPlot$llr[dataPlot$type == "trueGP"]), colour = "blue") + 
	scale_discrete_manual(aesthetics = "point_colour", values =c("red", "blue")) + ggtitle(paste("DWOR CHNK", "92", "SNPs"))

library(cowplot)

plot_grid(plotlist = list(fewMarkersRidge, manyMarkersRidge), nrow = 2, labels =NULL)

temp <- Gma
temp[1:100] <- -9

system.time(
multWithError(kid = kidGenos, gMa = temp, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
probUerror(kid = kidGenos, gMa = temp, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
)
   # user  system elapsed 
   # 2.36    0.66    3.04 
# ouch! need to rewrite in c++

logDirichMultPMF(c(1,0,1,0,0), c(10000,10,10,10,10))

testData <- read.csv("all_sockeye_BCB.txt", sep = "\t", stringsAsFactors = FALSE)
testData <- testData[testData$Pedigree == "OneEAGL14BCB",]

save(testData, file = "testDataFilter.rda") #saving for upload to github for access elswhere


load("testDataFilter.rda")

genos <- testData[,28:ncol(testData)]
base <- genos[,1]
base[base == "00"] <- NA
base2 <- as.data.frame(cbind("pop", paste0("indiv", 1:length(base)), substr(base,1,1), substr(base,2,2)), stringsAsFactors = FALSE)
colnames(base2) <- c("population", "id", "locus1", "locus1.A2")
dropout <- data.frame(locus = "locus1", allele = c("A", "B"), dropoutRate = .01, stringsAsFactors = FALSE)
# dropout <- data.frame(locus = "locus1", allele = c("A", "B"), dropoutRate = 0, stringsAsFactors = FALSE)
snpError <- data.frame(locus = "locus1", order = 1, numBases = 2, error = .005, stringsAsFactors = FALSE)

#SNP
createGmaInput(baseline = base2, mixture = base2[,2:ncol(base2)], unsampledPops = NULL, perSNPerror = snpError, dropoutProb = dropout)

#microhap
base2$locus1[!is.na(base2$locus1)] <- paste0(base2$locus1[!is.na(base2$locus1)], sample(c("C", "G"), sum(!is.na(base2$locus1)), replace = TRUE))
base2$locus1.A2[!is.na(base2$locus1)] <- paste0(base2$locus1.A2[!is.na(base2$locus1)], sample(c("C", "G"), sum(!is.na(base2$locus1)), replace = TRUE))
uAllele <- unique(na.omit(c(base2$locus1, base2$locus1.A2)))
dropout <- data.frame(locus = "locus1", allele = uAllele, dropoutRate = .01, stringsAsFactors = FALSE)
snpError <- data.frame(locus = "locus1", order = 1:unique(nchar(uAllele)), numBases = 2, error = .005, stringsAsFactors = FALSE)


base2$locus1[!is.na(base2$locus1)] <- paste0(base2$locus1[!is.na(base2$locus1)], sample(c("C", "G", "T"), sum(!is.na(base2$locus1)), replace = TRUE))
base2$locus1.A2[!is.na(base2$locus1)] <- paste0(base2$locus1.A2[!is.na(base2$locus1)], sample(c("C", "G", "T"), sum(!is.na(base2$locus1)), replace = TRUE))
uAllele <- unique(na.omit(c(base2$locus1, base2$locus1.A2)))
dropout <- data.frame(locus = "locus1", allele = uAllele, dropoutRate = .01, stringsAsFactors = FALSE)
snpError <- data.frame(locus = "locus1", order = 1:unique(nchar(uAllele)), numBases = c(2,2,3), error = .005, stringsAsFactors = FALSE)


base2 <- base2[1:10,]
base2 <- base2[c(1,5),]

test <- createGmaInput(baseline = base2, mixture = base2[,2:ncol(base2)], unsampledPops = NULL, perSNPerror = snpError, dropoutProb = dropout)

str(test)
table(test$baseline$locus1)
head(test$baseline)
head(test$mixture)
test$unsampledPops
test$genotypeErrorRates
rownames(test$genotypeErrorRates[[1]])
rowSums(test$genotypeErrorRates[[1]])
test$genotypeKeys
test$baselineParams
test$unsampledPopsParams
test$alleleKeys

test

res <- inferGrandma(test, "ssGP")
res[res$Kid == 2,]



Rcpp::cppFunction("double testCpp(int b){
double a = 5.5;			 
a++;
return a;
}
			 ")

testCpp(12323)


Rcpp::evalCpp("false + !false")

str(knownGPs)

mixture <- knownGPs$Individual.Name
baseline <- c(knownGPs$trueGpa, knownGPs$trueGma)

mixture <- genos[mixture,]
baseline <- genos[baseline,]

boolKeep <- apply(rbind(mixture, baseline), 2, function(x) any(!is.na(x)))
mixture <- mixture[,boolKeep]
baseline <- baseline[,boolKeep]

unique(substr(rownames(baseline),1,12))
baselineIn <- data.frame(substr(rownames(baseline),1,12), rownames(baseline), stringsAsFactors = FALSE)

for(i in 1:ncol(baseline)){
	baseline[baseline[,i] == -9,i] <- NA
	baseline[baseline[,i] == 0,i] <- "AA"
	baseline[baseline[,i] == 1,i] <- "AB"
	baseline[baseline[,i] == 2,i] <- "BB"
	baselineIn <- cbind(baselineIn, substr(baseline[,i],1,1), substr(baseline[,i],2,2), stringsAsFactors = FALSE)
	colnames(baselineIn)[(ncol(baselineIn) - 1) : ncol(baselineIn)] <- c(paste0("Locus_", i), paste0("Locus_", i, ".A2"))
}

mixtureIn <- data.frame(rownames(mixture), stringsAsFactors = FALSE)
for(i in 1:ncol(mixture)){
	mixture[mixture[,i] == -9,i] <- NA
	mixture[mixture[,i] == 0,i] <- "AA"
	mixture[mixture[,i] == 1,i] <- "AB"
	mixture[mixture[,i] == 2,i] <- "BB"
	mixtureIn <- cbind(mixtureIn, substr(mixture[,i],1,1), substr(mixture[,i],2,2), stringsAsFactors = FALSE)
	colnames(mixtureIn)[(ncol(mixtureIn) - 1) : ncol(mixtureIn)] <- c(paste0("Locus_", i), paste0("Locus_", i, ".A2"))
}

colnames(baselineIn)[3:ncol(baselineIn)] == colnames(mixtureIn)[2:ncol(mixtureIn)]

dropout <- data.frame()
snpError <- data.frame()
for (i in seq(3, ncol(baselineIn) - 1, 2)){
	uAllele <- unique(na.omit(c(baselineIn[,i], baselineIn[,i+1], mixtureIn[,i-1], mixtureIn[,i])))
	dropout <- rbind(dropout, data.frame(
		locus = colnames(baselineIn)[i],
		allele = uAllele,
		dropoutRate = .005,
		stringsAsFactors = FALSE
	), stringsAsFactors = FALSE)
	
	snpError <- rbind(snpError, data.frame(
		locus = colnames(baselineIn)[i],
		order = 1:unique(nchar(uAllele)),
		numBases = 2,
		error = .005,
		stringsAsFactors = FALSE
	), stringsAsFactors = FALSE)
	
	
}

test <- createGmaInput(baseline = baselineIn, mixture = mixtureIn[1:10,], unsampledPops = NULL, perSNPerror = snpError, dropoutProb = dropout)
