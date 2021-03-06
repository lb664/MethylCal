#' NFkB-related genes in celiac disease
#'
#' Dataset containing human genomic control DNA measured at eight 
#' distinct actual methylation percentages (0\%, 12.5\%, 25\%, 37.5\%,
#' 50\%, 62.5\%, 87.5\% and 100\%) obtained by mixing different ratios
#' of unmethylated and fully methylated human control DNA, 13 controls
#' (C), 17 celiac patients at the time of diagnosis (D) and the same
#' patients after two years of treatment with a gluten-free diet (T).
#' The patient data was pyrosequenced in three runs.
#'
#' @docType data
#'
#' @usage data(Celiac_data)
#'
#' @format A data frame consisting of eight NFkB-related genes (MAP3K7,
#' TNFAIP3, RELA, NFKBIA, TRADD, MAP3K14, MALT1 and MAP3K7IP1/TAB1) 
#' containing 35 CpG positions and 58 variables:
#' \describe{
#'   \item{Target}{character: DMR/CpG island name}
#'   \item{Chrom}{integer: Chromosome}
#'   \item{CpG_pos}{integer: CpG position}
#'   \item{0.0\%, 12.5\%, 25.0\%, 37.5\%, 50.0\%, 62.5\%, 87.5\%, 100.0\%}{numeric: Actual Methylation Percentages (AMPs)}
#'   \item{10C, 11C, 12C, 13C, 01C, 02C, 03C, 04C, 05C, 06C, 07C, 08C, 09C}{numeric: Healthy controls}
#'   \item{10D, 11D, 12D, 13D, 14D, 15D, 16D, 17D, 01D, 02D, 03D, 04D, 05D, 06D, 07D, 08D, 09D}{numeric: Celiac patients at the time of diagnosis}
#'   \item{10T, 11T, 12T, 13T, 14T, 15T, 16T, 17T, 01T, 02T, 03T, 04T, 05T, 06T, 07T, 08T, 09T}{numeric: Same patients after two years of treatment with a gluten-free diet}
#' }
#'
#' @keywords Celiac dataset
#'
#' @references Fernandez-Jimenez N, Castellanos-Rubio A, Plaza-Izurieta
#' L, Irastorza I, Elcoroaristizabal X, Jauregi-Miguel A, Jauregi-Miguel
#' A, Lopez-Euba T, Tutau C, de Pancorbo MM, Vitoria JC, Ramon Bilbao J.
#' Coregulation and modulation of NFkB-related genes in celiac disease:
#' uncovered aspects of gut mucosal inflammation. Hum Mol Genet. 2014;
#' 23(5):1298-1310. (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3919015/}{PMC})
#'
#' @references Fernandez-Jimenez N, Plaza-Izurieta L, Lopez-Euba T, 
#' Jauregi-Miguel A, Ramon Bilbao J. Cubic regression-based degree of 
#' correction predicts the performance of whole bisulfitome amplified 
#' DNA methylation analysis. Epigenetics. 2012; 7(12):1349-54. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/23154537}{PubMed})
#'
#' @examples
#' data(Celiac_data)
#' head(Celiac_data)
"Celiac_data"