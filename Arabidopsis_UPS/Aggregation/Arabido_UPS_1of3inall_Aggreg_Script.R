library(tidyverse)
library(doMC)
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
MetadataMQ <- read.delim("Arabidopsis_UPS/DATA/metadataMQ.txt")
# Peptide data loading (MQ export with MBR)
peptidesMQ <- read.delim("Arabidopsis_UPS/DATA/noMBR/peptides.txt")

# -- PRE-PROCESSING -- #
peptidesMQ[,grep(pattern = "Intensity.",colnames(peptidesMQ))][peptidesMQ[,grep(pattern = "Intensity.",
                                                                                colnames(peptidesMQ))]==0]=NA
# At least 1 quantified value in each condition
peptidesMQ_clean<-peptidesMQ[which(apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point1",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point2",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point3",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point4",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point5",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point6",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point7",colnames(peptidesMQ))]),1,sum)<3),]

peptidesMQ_clean<-subset(peptidesMQ_clean, subset = peptidesMQ_clean$Reverse!="+" & 
                           peptidesMQ_clean$Potential.contaminant!="+" &
                           peptidesMQ_clean$Unique..Groups.=="yes")

# log2 transformation
peptidesMQ_clean[,grep(pattern="Intensity.", x=colnames(peptidesMQ_clean))] <- log2(peptidesMQ_clean[,grep(pattern="Intensity.", x=colnames(peptidesMQ_clean))])


# -- NORMALISATION -- #
data.pept.norm <- preprocessCore::normalize.quantiles(as.matrix(peptidesMQ_clean[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean))]))
peptidesMQ_clean_norm <- peptidesMQ_clean
peptidesMQ_clean_norm[,grep(pattern = "Intensity.",colnames(peptidesMQ_clean_norm))] <- data.pept.norm

# -- PEPTIDE-LEVEL IMPUTATION -- #
imp.list <- sapply(as.character(1:6),  function(aaa) {
  imp.int = imp4p::impute.mle(tab = as.matrix(data.pept.norm), conditions = as.factor(MetadataMQ$Condition))
  tab.final = peptidesMQ_clean_norm
  tab.final[,grep(pattern = "Intensity.",colnames(tab.final))] = imp.int
  return(tab.final)
},
simplify = F, USE.NAMES = T)

# For each imputed dataset, sum by proteinGroups all unique peptides' intensities
db  = imp.list %>% 
  bind_rows(.id = 'Imput') %>%
  group_by(Imput, Leading.razor.protein) %>% 
  summarise(across(starts_with("Intensity."), ~ log2(sum(2^(.x))) ) )

data.imp.ag <- array(unlist(lapply(levels(factor(db$Imput)), function(aaa){db[db$Imput==aaa,-c(1,2)]})),
                     dim = c(nrow(lapply(levels(factor(db$Imput)), function(aaa){db[db$Imput==aaa,-c(1,2)]})[[1]]), 
                             ncol(lapply(levels(factor(db$Imput)), function(aaa){db[db$Imput==aaa,-c(1,2)]})[[1]]), 
                             length(lapply(levels(factor(db$Imput)), function(aaa){db[db$Imput==aaa,-c(1,2)]}))))

# -- ESTIMATION -- #
VarRubin.mat <- rubin2.all(data = data.imp.ag, metacond = as.factor(MetadataMQ$Condition),
                           is.parallel = F)

# -- PROJECTION -- #
VarRubin.S2 <- proj_matrix(VarRubin.matrix = VarRubin.mat, metadata = MetadataMQ)

# -- MODERATED T-TEST -- #
res.mi4limma <- mi4limma(qData = apply(data.imp.ag,1:2,mean), 
                         sTab = MetadataMQ, 
                         VarRubin = sqrt(VarRubin.S2))

res.dapar <- limmaCompleteTest.mod(qData = apply(data.imp.ag,1:2,mean), 
                                   sTab = MetadataMQ)

# -- BENJAMINI-HOCHBERG 1% FDR -- #
dapar_0.05fmol_vs_10fmol <- adjust.p(res.dapar$res.l$P_Value$Point1_vs_Point7_pval, alpha = 0.01)
mi4limma_0.05fmol_vs_10fmol <- adjust.p(res.mi4limma$P_Value$Point1_vs_Point7_pval, alpha = 0.01)

dapar_0.25fmol_vs_10fmol <- adjust.p(res.dapar$res.l$P_Value$Point2_vs_Point7_pval, alpha = 0.01)
mi4limma_0.25fmol_vs_10fmol <- adjust.p(res.mi4limma$P_Value$Point2_vs_Point7_pval, alpha = 0.01)

dapar_0.5fmol_vs_10fmol <- adjust.p(res.dapar$res.l$P_Value$Point3_vs_Point7_pval, alpha = 0.01)
mi4limma_0.5fmol_vs_10fmol <- adjust.p(res.mi4limma$P_Value$Point3_vs_Point7_pval, alpha = 0.01)

dapar_1.25fmol_vs_10fmol <- adjust.p(res.dapar$res.l$P_Value$Point4_vs_Point7_pval, alpha = 0.01)
mi4limma_1.25fmol_vs_10fmol <- adjust.p(res.mi4limma$P_Value$Point4_vs_Point7_pval, alpha = 0.01)

dapar_2.5fmol_vs_10fmol <- adjust.p(res.dapar$res.l$P_Value$Point5_vs_Point7_pval, alpha = 0.01)
mi4limma_2.5fmol_vs_10fmol <- adjust.p(res.mi4limma$P_Value$Point5_vs_Point7_pval, alpha = 0.01)

dapar_5fmol_vs_10fmol <- adjust.p(res.dapar$res.l$P_Value$Point6_vs_Point7_pval, alpha = 0.01)
mi4limma_5fmol_vs_10fmol <- adjust.p(res.mi4limma$P_Value$Point6_vs_Point7_pval, alpha = 0.01)
