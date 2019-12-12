


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

# now allele counts in baseline pop
refFreqList <- apply(genos[baseline$Individual.Name,],2,function(x) sum(x == 1) + 2*sum(x == 0))
altFreqList <- apply(genos[baseline$Individual.Name,],2,function(x) sum(x == 1) + 2*sum(x == 2))

# define mixture
mixture <- testData$Individual.Name[testData$Pedigree == "OtsDWOR18S"]
mixture <- mixture[mixture %in% rownames(genos)]

#now grandparentage
results <- matrix(nrow = length(mixture), ncol=nrow(crosses))
for(k in 1:length(mixture)){
	kidGenos <- genos[mixture[k],]
	for(c in 1:nrow(crosses)){
		Gma <- genos[crosses[c,1],]
		Gpa <- genos[crosses[c,2],]
		
		results[k,c] <- multWithError(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep(.005,length(kidGenos))) -
			probUerror(kid = kidGenos, gMa = Gma, gPa = Gpa, altFreqList = altFreqList, refFreqList = refFreqList, epsList = rep (.005,length(kidGenos)))
	}
}

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
