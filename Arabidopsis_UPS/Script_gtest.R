## ---------------- ##
## G-test filtering ##
## ---------------- ##

## Load the data from the Arabidopsis + UPS experiment (with MBR)
# Metadata
MetadataMQ <- read.delim("Arabidopsis_UPS/DATA/metadataMQ.txt")
# Peptide data loading (MQ export with MBR)
peptidesMQ <- read.delim("Arabidopsis_UPS/DATA/MBR/peptides.txt")
# 23125 peptides
## Replace "0" by "NA" in intensity values

rownames(peptidesMQ) <- mapply(paste,peptidesMQ$Leading.razor.protein,1:length(peptidesMQ$Leading.razor.protein),MoreArgs=list(sep="_"))

peptidesMQ[,grep(pattern = "Intensity.",colnames(peptidesMQ))][peptidesMQ[,grep(pattern = "Intensity.",
                                                                                colnames(peptidesMQ))]==0]=NA


## Remove reverse peptide sequences and contaminants.
peptidesMQ_clean<-subset(peptidesMQ, subset = peptidesMQ$Reverse!="+" & 
                           peptidesMQ$Potential.contaminant!="+")
# 22829 peptides (296 peptides were removed)
## Define the quantitative data matrix
qData = peptidesMQ_clean[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean))]
## For each peptide, get a table counting missing values by condition.
tableNA.qData <- apply(is.na(qData),1,table,MetadataMQ$Condition)
tableNA.qData[[1]]
tableNA.qData[[3]]
## For each peptide, are intensities in some condition missing ?
id.mix <- unlist(lapply(tableNA.qData,function(res) nrow(res)>1))
head(id.mix)
length(id.mix[id.mix])
# There are 9797 peptides for which some intensities are missing.
## Extract p-values of the g-test for the 9797 peptides which some intensities
## values are missing.
res.g.test <- cbind(rownames=as.data.frame(rownames(qData)[id.mix]),
                    p.val=apply(is.na(qData[id.mix,]),1,
                                function(tab) return(ProteoMM::g.test(x=tab,y=MetadataMQ$Condition)$p.value)))

p.res.g.test=res.g.test[,2]
names(p.res.g.test) <- rownames(res.g.test)

res.g.test[p.res.g.test<0.05,]
str(res.g.test)
head(res.g.test)
sum(p.res.g.test<0.05)
#2212 without correction for multiple tests

p.adjust(p.res.g.test, method = "BH")
sum(p.adjust(p.res.g.test, method = "BH")<0.05)
#488 after correction for multiple tests



## Perform g-test after "our" filtering criterion
## Remove reverse peptide sequences and contaminants.
peptidesMQ_filter<-subset(peptidesMQ, subset = peptidesMQ$Reverse!="+" & 
                            peptidesMQ$Potential.contaminant!="+")
# 22829 peptides (296 peptides were removed)
## At least 1 quantified value in each condition
peptidesMQ_filter<-peptidesMQ_filter[which(apply(is.na(peptidesMQ_filter[,grep(pattern = "Intensity.Point1",colnames(peptidesMQ_filter))]),1,sum)<3
                                   & apply(is.na(peptidesMQ_filter[,grep(pattern = "Intensity.Point2",colnames(peptidesMQ_filter))]),1,sum)<3
                                   & apply(is.na(peptidesMQ_filter[,grep(pattern = "Intensity.Point3",colnames(peptidesMQ_filter))]),1,sum)<3
                                   & apply(is.na(peptidesMQ_filter[,grep(pattern = "Intensity.Point4",colnames(peptidesMQ_filter))]),1,sum)<3
                                   & apply(is.na(peptidesMQ_filter[,grep(pattern = "Intensity.Point5",colnames(peptidesMQ_filter))]),1,sum)<3
                                   & apply(is.na(peptidesMQ_filter[,grep(pattern = "Intensity.Point6",colnames(peptidesMQ_filter))]),1,sum)<3
                                   & apply(is.na(peptidesMQ_filter[,grep(pattern = "Intensity.Point7",colnames(peptidesMQ_filter))]),1,sum)<3),]
# 14321 peptides (8508 peptides removed)
## Define the quantitative data matrix
qData.f = peptidesMQ_filter[,grep(pattern = "Intensity.",colnames(peptidesMQ_filter))]
## For each peptide, get a table counting missing values by condition.


tableNA.qData.f <- apply(is.na(qData.f),1,table,MetadataMQ$Condition)
tableNA.qData.f[[1]]
tableNA.qData.f[[3]]
## For each peptide, are intensities in some condition missing ?
id.mix.f <- unlist(lapply(tableNA.qData.f,function(res) nrow(res)>1))
#id.mix is empty
head(id.mix.f)
length(id.mix.f[id.mix.f])
# There are 6250 peptides for which some intensities are missing.
## Extract p-values of the g-test for the 6250 peptides which some intensities
## values are missing.
res.g.test.f <- cbind(rownames=as.data.frame(rownames(qData.f)[id.mix.f]),
                    p.val=apply(is.na(qData.f[id.mix.f,]),1,
                                function(tab) return(ProteoMM::g.test(x=tab,y=MetadataMQ$Condition)$p.value)))



p.res.g.test.f=res.g.test.f[,2]
names(p.res.g.test.f) <- rownames(res.g.test.f)

res.g.test.f[p.res.g.test.f<0.05,]
str(res.g.test.f)
head(res.g.test.f)
sum(p.res.g.test.f<0.05)
#84 without correction for multiple tests

p.adjust(p.res.g.test.f, method = "BH")
sum(p.adjust(p.res.g.test.f, method = "BH")<0.05)
#0 after correction for multiple tests



# Analysis

length(grep("ARATH",rownames(res.g.test[p.res.g.test<0.05,])))
# 1767 peptides ARATH found significant by raw g.test

length(grep("HUMAN",rownames(res.g.test[p.res.g.test<0.05,])))
# 445 peptides HUMAN found significant by g.test corrected for multiple tests

#-> raw g.test found significant 1767/445=3.97 ARATH than HUMAN
# g.test increases false positives. Need for an adapted g.test for specific experimental designs such as the ones than we analyzed.

p.adjust(p.res.g.test, method = "BH")
sum(p.adjust(p.res.g.test, method = "BH")<0.05)
#488 after correction for multiple tests

length(grep("ARATH",rownames(res.g.test[p.adjust(p.res.g.test, method = "BH")<0.05,])))
# 280 peptides ARATH found significant by raw g.test

length(grep("HUMAN",rownames(res.g.test[p.adjust(p.res.g.test, method = "BH")<0.05,])))
# 208 peptides HUMAN found significant by raw g.test corrected for multiple tests

#-> corrected g.test found significant 280/208=1.35 ARATH than HUMAN
# g.test increases false positives. Need for an adapted g.test for specific experimental designs such as the ones than we analyzed.


nrow((res.g.test[res.g.test[p.res.g.test<0.05,1],])[rownames(res.g.test[res.g.test[p.res.g.test<0.05,1],]) %in% names(tableNA.qData.f),])
#84 peptides found significant by raw g.test after filtering

length(grep("ARATH",rownames((res.g.test[res.g.test[p.res.g.test<0.05,1],])[rownames(res.g.test[res.g.test[p.res.g.test<0.05,1],]) %in% names(tableNA.qData.f),])))
# After filtering 76 peptides ARATH found significant by raw g.test

length(grep("HUMAN",rownames((res.g.test[res.g.test[p.res.g.test<0.05,1],])[rownames(res.g.test[res.g.test[p.res.g.test<0.05,1],]) %in% names(tableNA.qData.f),])))
# After filtering 8 peptides HUMAN found significant by g.test corrected for multiple tests

#-> raw g.test found significant 76/8=9.5 ARATH than HUMAN
# g.test increases false positives. Need for an adapted g.test for specific experimental designs such as the ones than we analyzed.

#No significant peptides for corrected g.test after filtering.




