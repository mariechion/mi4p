#' A single simulated dataset
#' 
#' This dataset was simulated using the default values of the 
#' values of the options of the protdatasim function and the set.seed value set 
#' to 4619.
#' 
#' @name datasim
#' @docType data
#' @format A data frame with 200 observations on the following 11 variables.
#' \describe{ 
#' \item{id.obs}{a numeric vector} 
#' \item{X1}{a numeric vector} 
#' \item{X2}{a numeric vector} 
#' \item{X3}{a numeric vector} 
#' \item{X4}{a numeric vector} 
#' \item{X5}{a numeric vector} 
#' \item{X6}{a numeric vector} 
#' \item{X7}{a numeric vector} 
#' \item{X8}{a numeric vector} 
#' \item{X9}{a numeric vector} 
#' \item{X10}{a numeric vector} 
#' }
#' @author M. Chion, Ch. Carapito and F. Bertrand.
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}. \doi{doi:10.1371/journal.pcbi.1010420}.
#' @source We simulated the data.
#' @keywords datasets
#' @examples
#' 
#' data(datasim)
#' str(datasim)
#' 
NULL

#' A list of simulated datasets.
#' 
#' This list of 100 datasets was simulated using the default 
#' values of the options of the protdatasim function and the set.seed value set 
#' to 4619.
#' 
#' @name norm.200.m100.sd1.vs.m200.sd1.list
#' @docType data
#' @format The format is: List of 100 data.frames.
#' \describe{ 
#' \item{data.frame}{200 obs. of 11 variables
#' \describe{ 
#' \item{id.obs}{int [1:200] 1 2 3 4 5 6 7 8 9 10 ...}  
#' \item{X1}{num [1:200] 99.6 99.9 100.2 99.8 100.4 ...}  
#' \item{X2}{num [1:200] 97.4 101.3 100.3 100.2 101.7 ...}  
#' \item{X3}{num [1:200] 100.3 100.9 99.1 101.2 100.6 ...}
#' \item{X4}{num [1:200] 99.4 99.2 98.5 99.1 99.5 ...}  
#' \item{X5}{num [1:200] 98.5 99.7 100 100.2 100.7 ...}  
#' \item{X6}{num [1:200] 200 199 199 200 199 ...}  
#' \item{X7}{num [1:200] 200 200 202 199 199 ...}  
#' \item{X8}{num [1:200] 202 199 200 199 201 ...}  
#' \item{X9}{num [1:200] 200 200 199 201 200 ...}  
#' \item{X10}{num [1:200] 200 198 200 201 199 ...}
#' }}
#' \item{attr(*, "metadata")}{'data.frame': 10 obs. of 3 variables:
#' \describe{ 
#' \item{Sample.name}{chr [1:10] "X1" "X2" "X3" "X4" ...}
#' \item{Condition}{Factor w/ 2 levels "A","B": 1 1 1 1 1 2 2 2 2 2} 
#' \item{Bio.Rep}{int [1:10] 1 2 3 4 5 6 7 8 9 10}
#' }}
#' \item{...}{...}
#' }
#' @author M. Chion, Ch. Carapito and F. Bertrand.
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}.  \doi{doi:10.1371/journal.pcbi.1010420}.
#' @source We simulated the data.
#' @keywords datasets
#' @examples
#' 
#' data(norm.200.m100.sd1.vs.m200.sd1.list)
#' str(norm.200.m100.sd1.vs.m200.sd1.list)
#' 
NULL

#' Extract of the abundances of Exp1_R25_pept dataset
#' 
#' The data frame \code{qData} contains the first 500 rows of six columns that are the quantitation of peptides for the six replicates. They were obtained using the code \code{exprs(Exp1_R25_pept)[1:500,]}.
#' 
#' The \code{DAPARdata}'s \code{Exp1_R25_pept} dataset is the final outcome of a quantitative mass spectrometry-based proteomic analysis of two samples containing different concentrations of 48 human proteins (UPS1 standard from Sigma-Aldrich) within a constant yeast background (see Giai Gianetto et al. (2016) for details). It contains the abundance values of the different human and yeast peptides identified and quantified in these two conditions. The two conditions represent the measured abundances of peptides when respectively 25 fmol and 10 fmol of UPS1 human proteins were mixed with the yeast extract before mass spectrometry analyses. This results in a concentration ratio of 2.5. Three technical replicates were acquired for each condition.
#' 
#' @name qData
#' @docType data
#' @format The format is: num [1:500, 1:6] 24.8 24.7 24.6 NA 24.5 ...  -
#' attr(*, "dimnames")=List of 2 ..$ : chr [1:500] "0" "1" "2" "3" ...  ..$ :
#' chr [1:6] "Intensity_C_R1" "Intensity_C_R2" "Intensity_C_R3"
#' "Intensity_D_R1" ...
#' @references 
#' Cox J., Hein M.Y., Luber C.A., Paron I., Nagaraj N., Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep, 13(9):2513-26.
#' 
#' Giai Gianetto, Q., Combes, F., Ramus, C., Bruley, C., Coute, Y., Burger, T. (2016). Calibration plot for proteomics: A graphical tool to visually check the assumptions underlying FDR control in quantitative experiments. Proteomics, 16(1), 29-32.
#' @source The \code{DAPARdata} package.
#' @keywords datasets
#' @examples
#' 
#' data(qData)
#' str(qData)
#' pairs(qData)
#'  
NULL


#' Experimental design for the Exp1_R25_pept dataset
#' 
#' The data frame \code{sTab}  contains the experimental design and gives few informations about the samples. They were obtained using the code \code{pData(Exp1_R25_pept)}.
#' 
#' The \code{DAPARdata}'s \code{Exp1_R25_pept} dataset is the final outcome of a quantitative mass spectrometry-based proteomic analysis of two samples containing different concentrations of 48 human proteins (UPS1 standard from Sigma-Aldrich) within a constant yeast background (see Giai Gianetto et al. (2016) for details). It contains the abundance values of the different human and yeast peptides identified and quantified in these two conditions. The two conditions represent the measured abundances of peptides when respectively 25 fmol and 10 fmol of UPS1 human proteins were mixed with the yeast extract before mass spectrometry analyses. This results in a concentration ratio of 2.5. Three technical replicates were acquired for each condition.
#' 
#' @name sTab
#' @docType data
#' @format A data frame with 6 observations on the following 3 variables.
#' \describe{ 
#' \item{Sample.name}{a character vector}
#' \item{Condition}{a character vector} 
#' \item{Bio.Rep}{a numeric vector} 
#' }
#' @references
#' Cox J., Hein M.Y., Luber C.A., Paron I., Nagaraj N., Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep, 13(9):2513-26.
#' 
#' Giai Gianetto, Q., Combes, F., Ramus, C., Bruley, C., Coute, Y., Burger, T. (2016). Calibration plot for proteomics: A graphical tool to visually check the assumptions underlying FDR control in quantitative experiments. Proteomics, 16(1), 29-32.
#' @source The \code{DAPARdata} package.
#' @keywords datasets
#' @examples
#' 
#' data(sTab)
#' str(sTab)
#' 
NULL

#' mm_peptides - peptide-level intensities for mouse
#'
#' A dataset containing the protein and petide information and peptide-level
#' intensities for 6 samples: 3 CG and 3 mCG groups. There are 69 proteins.
#' The columns are as follows:
#'
#' \itemize{
#'   \item Sequence - peptide sequence - randomly chosen from a larger list of
#'         sequences
#'   \item MatchedID - numeric ID that links proteins in the two datasets,
#'         unnecessary if datasets are for the same species
#'   \item ProtID - protein ID, artificial protein ID, eg. Prot1, Prot2, ...
#'   \item GeneID - gene ID, artificial gene ID, eg. Gene1, Gene2, ...
#'   \item ProtName - artificial Protein Name
#'   \item ProtIDLong - long protein ID, full protein name, here artificially
#'         simulated
#'   \item GeneIDLong - long gene ID, full gene name, here artificially
#'         simulated
#'   \item CG1 - raw intensity column for sample 1 in CG group
#'   \item CG2 - raw intensity column for sample 2 in CG group
#'   \item CG3 - raw intensity column for sample 3 in CG group
#'   \item mCG1 - raw intensity column for sample 1 in mCG group
#'   \item mCG2 - raw intensity column for sample 2 in mCG group
#'   \item mCG3 - raw intensity column for sample 3 in mCG group
#' }
#'
#' @docType data
#' @keywords datasets
#' @name mm_peptides
#' @usage data(mm_peptides)
#' @format A data frame with 1102 rows and 13 colummns, compiring 7 columns of
#'         metadata and 6 columns of peptide intensities. 69 proteins.
NULL
