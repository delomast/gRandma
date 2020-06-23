# creating example data file

# loading STHD mix of mh and snp loci
data_mh_snp <- read.table("S:\\Eagle Fish Genetics Lab\\Tom\\microhap infrastructure\\EFGLmh\\example_snp_mh.txt", 
								  sep = "\t", header = TRUE, na.strings = c("0", "00", "000", ""), colClasses = "character")
# remove metadata
data_mh_snp <- data_mh_snp[,-(3:8)]
# remove some failed inds
gS <- apply(data_mh_snp[,3:ncol(data_mh_snp)], 1, function(x) sum(is.na(x)) / length(x))
data_mh_snp <- data_mh_snp[gS <= .1,]
table(data_mh_snp$Pedigree)
# remove small pop
data_mh_snp <- data_mh_snp[data_mh_snp$Pedigree != "OmyEFSW19S",]
# change names to be generic 
data_mh_snp$Individual.Name <- paste0("Ind", 1:nrow(data_mh_snp))
data_mh_snp$Pedigree <- paste0("Pop_", as.numeric(as.factor(data_mh_snp$Pedigree)))
table(data_mh_snp$Pedigree)
colnames(data_mh_snp)[1:2] <- c("Pop", "Ind")
nLoci <- (ncol(data_mh_snp) - 2) / 2
colnames(data_mh_snp)[seq(3,ncol(data_mh_snp) - 1, 2)] <- paste0("Locus_", 1:nLoci)
colnames(data_mh_snp)[seq(4,ncol(data_mh_snp), 2)] <- paste0("Locus_", 1:nLoci, ".A2")
# remove loci with all fails or no variation
to_remove <- c()
for(i in seq(3,ncol(data_mh_snp) - 1, 2)){
	a <- c(data_mh_snp[,i], data_mh_snp[,i+1])
	a <- unique(a)
	a <- a[!is.na(a)]
	if(length(a) < 2) to_remove <- c(to_remove, i)
}
to_remove <- c(to_remove, to_remove + 1)
data_mh_snp <- data_mh_snp[,-to_remove]
# add data to package
usethis::use_data(data_mh_snp, overwrite = TRUE, internal = FALSE)

# add vignette
usethis::use_vignette("Load_in_data_and_gmaData_structure")
