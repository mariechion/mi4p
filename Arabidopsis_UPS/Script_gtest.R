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
res.g.test[res.g.test[,2]<0.05,]
str(res.g.test)
head(res.g.test)


## Perform g-test after "our" filtering criterion
## At least 1 quantified value in each condition
peptidesMQ_filter<-peptidesMQ[which(apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point1",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point2",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point3",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point4",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point5",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point6",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point7",colnames(peptidesMQ))]),1,sum)<3),]
# 14492 peptides (8633 peptides removed)
## Remove reverse peptide sequences and contaminants.
peptidesMQ_filter<-subset(peptidesMQ_filter, subset = peptidesMQ_filter$Reverse!="+" & 
                           peptidesMQ_filter$Potential.contaminant!="+")
# 14321 peptides (171 peptides removed)
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
res.g.test[res.g.test[,2]<0.05,]