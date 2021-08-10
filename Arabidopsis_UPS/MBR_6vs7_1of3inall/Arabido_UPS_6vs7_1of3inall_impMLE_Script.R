library("doMC")
registerDoMC(detectCores(all.tests = FALSE, logical = TRUE)-1)

# Load functions
source("Functions/multi_impute.R")
source("Functions/rubin2_all.R")
source("Functions/proj_matrix.R")
source("Functions/mi4limma.R")
source("Functions/limmaCompleteTest_mod.R")

# Random number generation
set.seed(17)

# Metadata
MetadataMQ <- data.frame(Sample.name = c(grep(pattern = "Intensity.Point6", x = colnames(peptidesMQ_clean), value=T), 
                                         grep(pattern = "Intensity.Point7", x = colnames(peptidesMQ_clean), value=T)),
                         Condition = factor(rep(c("Point6", "Point7"), c(3,3))),
                         Bio.Rep = 1:6)

# Peptide data loading (MQ export with MBR)
peptidesMQ <- read.delim("Arabidopsis_UPS/DATA/MBR/peptides.txt")

# -- PRE-PROCESSING -- #
peptidesMQ[,grep(pattern = "Intensity.",colnames(peptidesMQ))][peptidesMQ[,grep(pattern = "Intensity.",
                                                                                colnames(peptidesMQ))]==0]=NA

# At least 1 quantified value in each condition
peptidesMQ_clean<-peptidesMQ[which(apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point6",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point7",colnames(peptidesMQ))]),1,sum)<3),]

peptidesMQ_clean<-subset(peptidesMQ_clean, subset = peptidesMQ_clean$Reverse!="+" & 
                           peptidesMQ_clean$Potential.contaminant!="+")

peptidesMQ_clean[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean))]<-
  log2(peptidesMQ_clean[,grep(pattern = "Intensity.", colnames(peptidesMQ_clean))])

data.pept <- subset(peptidesMQ_clean, select = c(grep(pattern = "Intensity.Point6", x = colnames(peptidesMQ_clean), value=T), 
                                                 grep(pattern = "Intensity.Point7", x = colnames(peptidesMQ_clean), value=T)))

sum(is.na(data.pept))/
  prod(dim(data.pept))*100
# 9.504985

# -- NORMALISATION -- #
data.pept.norm <- preprocessCore::normalize.quantiles(as.matrix(data.pept))

# -- IMPUTATION -- #
data.pept.imp <- multi.impute(data = data.pept.norm,
                              conditions = MetadataMQ$Condition,
                              method = "MLE")

# -- ESTIMATION -- #
VarRubin.mat <- rubin2.all(data = data.pept.imp, metacond = MetadataMQ$Condition,
                           is.parallel = T)

# -- PROJECTION -- #
VarRubin.S2 <- proj_matrix(VarRubin.matrix = VarRubin.mat, metadata = MetadataMQ)

# -- MODERATED T-TEST -- #
res.mi4limma <- mi4limma(qData = apply(data.pept.imp,1:2,mean), 
                         sTab = MetadataMQ, 
                         VarRubin = sqrt(VarRubin.S2))

res.dapar <- limmaCompleteTest.mod(qData = apply(data.pept.imp,1:2,mean), 
                                   sTab = MetadataMQ)

# -- BENJAMINI-HOCHBERG 1% FDR -- #
dapar_5fmol_vs_10fmol <- adjust.p(res.dapar$res.l$P_Value$Point6_vs_Point7_pval, alpha = 0.01)
mi4limma_5fmol_vs_10fmol <- adjust.p(res.mi4limma$P_Value$Point6_vs_Point7_pval, alpha = 0.01)
