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
MetadataMQ <- read.delim("Yeast_UPS/DATA/metadataMQ.txt")
# Peptide data loading (MQ export with MBR)
peptidesMQ <- read.delim("Yeast_UPS/DATA/peptides.txt")

# -- PRE-PROCESSING -- #
peptidesMQ[,grep(pattern = "Intensity.",colnames(peptidesMQ))][peptidesMQ[,grep(pattern = "Intensity.",
                                                                                colnames(peptidesMQ))]==0]=NA

# At least 1 quantified value in each condition
peptidesMQ_clean<-peptidesMQ[which(apply(is.na(peptidesMQ[,grep(pattern = "Intensity.TG.0.5fmol",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.TG.1fmol",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.TG.2.5fmol",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.TG.5fmol",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.TG.10fmol",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.TG.25fmol",colnames(peptidesMQ))]),1,sum)<3),]

peptidesMQ_clean<-subset(peptidesMQ_clean, subset = peptidesMQ_clean$Reverse!="+" & 
                           peptidesMQ_clean$Potential.contaminant!="+")

sum(is.na(peptidesMQ_clean[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean))]))/
  prod(dim(peptidesMQ_clean[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean))]))*100 
#  5.694311%

peptidesMQ_clean[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean))]<-
  log2(peptidesMQ_clean[,grep(pattern = "Intensity.", colnames(peptidesMQ_clean))])


# -- NORMALISATION -- #
data.pept.norm <- preprocessCore::normalize.quantiles(as.matrix(peptidesMQ_clean[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean))]))
peptidesMQ_clean_norm <- peptidesMQ_clean
peptidesMQ_clean_norm[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean_norm))] <- data.pept.norm


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
dapar_0.5fmol_vs_25fmol <- adjust.p(res.dapar$res.l$P_Value$` 0.5fmol_vs_25fmol_pval`, alpha = 0.01)
mi4limma_0.5fmol_vs_25fmol <- adjust.p(res.mi4limma$P_Value$` 0.5fmol_vs_25fmol_pval`, alpha = 0.01)

dapar_1fmol_vs_25fmol <- adjust.p(res.dapar$res.l$P_Value$`1fmol_vs_25fmol_pval`, alpha = 0.01)
mi4limma_1fmol_vs_25fmol <- adjust.p(res.mi4limma$P_Value$`1fmol_vs_25fmol_pval`, alpha = 0.01)

dapar_2.5fmol_vs_25fmol <- adjust.p(res.dapar$res.l$P_Value$`2.5fmol_vs_25fmol_pval`, alpha = 0.01)
mi4limma_2.5fmol_vs_25fmol <- adjust.p(res.mi4limma$P_Value$`2.5fmol_vs_25fmol_pval`, alpha = 0.01)

dapar_5fmol_vs_25fmol <- adjust.p(res.dapar$res.l$P_Value$`25fmol_vs_5fmol_pval`, alpha = 0.01)
mi4limma_5fmol_vs_25fmol <- adjust.p(res.mi4limma$P_Value$`25fmol_vs_5fmol_pval`, alpha = 0.01)

dapar_10fmol_vs_25fmol <- adjust.p(res.dapar$res.l$P_Value$`10fmol_vs_25fmol_pval`, alpha = 0.01)
mi4limma_10fmol_vs_25fmol <- adjust.p(res.mi4limma$P_Value$`10fmol_vs_25fmol_pval`, alpha = 0.01)
