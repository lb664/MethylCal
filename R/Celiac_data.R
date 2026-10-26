#' NFkB-related genes in celiac disease
#'
#' Dataset containing human genomic control DNA measured at eight distinct
#' actual methylation percentages (0\%, 12.5\%, 25\%, 37.5\%, 50\%, 62.5\%,
#' 87.5\% and 100\%) obtained by mixing different ratios of unmethylated and
#' fully methylated human control DNA, 13 controls (C), 17 celiac patients at
#' the time of diagnosis (D) and the same patients after two years of
#' treatment with a gluten-free diet (T). The patient data were pyrosequenced
#' in three runs.
#'
#' @docType data
#'
#' @usage data(Celiac_data)
#'
#' @format A data frame consisting of eight NFkB-related genes (MAP3K7,
#'   TNFAIP3, RELA, NFKBIA, TRADD, MAP3K14, MALT1 and MAP3K7IP1/TAB1)
#'   containing 35 CpG positions and 58 variables:
#' \describe{
#'   \item{Target}{character: DMR/CpG island name}
#'   \item{Chrom}{integer: Chromosome}
#'   \item{CpG_pos}{integer: CpG position}
#'   \item{0.0\%, 12.5\%, 25.0\%, 37.5\%}{numeric: Actual Methylation Percentages (AMPs)}
#'   \item{50.0\%, 62.5\%, 87.5\%, 100.0\%}{numeric: Actual Methylation Percentages (AMPs)}
#'   \item{10C, 11C, 12C, 13C, 01C}{numeric: Healthy controls}
#'   \item{02C, 03C, 04C, 05C, 06C}{numeric: Healthy controls (continued)}
#'   \item{07C, 08C, 09C}{numeric: Healthy controls (continued)}
#'   \item{10D, 11D, 12D, 13D, 14D}{numeric: Celiac patients at the time of diagnosis}
#'   \item{15D, 16D, 17D, 01D, 02D}{numeric: Celiac patients at diagnosis (continued)}
#'   \item{03D, 04D, 05D, 06D, 07D}{numeric: Celiac patients at diagnosis (continued)}
#'   \item{08D, 09D}{numeric: Celiac patients at diagnosis (continued)}
#'   \item{10T, 11T, 12T, 13T, 14T}{numeric: Same patients after two years of treatment with a gluten-free diet}
#'   \item{15T, 16T, 17T, 01T, 02T}{numeric: Same patients after treatment (continued)}
#'   \item{03T, 04T, 05T, 06T, 07T}{numeric: Same patients after treatment (continued)}
#'   \item{08T, 09T}{numeric: Same patients after treatment (continued)}
#' }
#'
#' @keywords Celiac dataset
#'
#' @references \insertRef{FernandezJimenez2014}{MethylCal}
#' @references \insertRef{FernandezJimenez2012}{MethylCal}
#'
#' @examples
#' \dontrun{
#'    # Load the Celiac dataset and inspect the first rows.
#'  
#'    data(Celiac_data)
#'    head(Celiac_data)
#' }
"Celiac_data"
