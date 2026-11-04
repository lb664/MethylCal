#' Patients potentially affected by the Beckwith-Wiedemann syndrome (BWS)
#'
#' Dataset containing human genomic control DNA (Promega), with five distinct
#' actual methylation percentages (0\%, 25\%, 50\%, 75\% and 100\%) obtained
#' by mixing different ratios of unmethylated and fully methylated human
#' control DNA, DNA from 15 healthy controls, and DNA from 18 patients.
#'
#' @docType data
#'
#' @usage data(BWS_data)
#'
#' @format A data frame consisting of nine independent assays in two CpG
#'   islands located on the SDHC gene promoter (SDHC_CpG:17 and SDHC_CpG:27) and
#'   seven imprinted Differentially Methylated Regions (PLAGL1, GRB10, MEST,
#'   NRAP, H19/IGF2, KCNQ1OT1 and MEG3) with 115 CpG positions and 41
#'   variables:
#' \describe{
#'   \item{Target}{character: DMR/CpG island name}
#'   \item{Chrom}{integer: Chromosome}
#'   \item{CpG_pos}{integer: CpG position}
#'   \item{0.0\%, 25.0\%, 50.0\%, 75.0\%, 100.0\%}{numeric: Actual Methylation Percentages (AMPs)}
#'   \item{NC15, NC13, NC06, NC05, NC03}{numeric: Healthy controls}
#'   \item{NC01, HC37, HC34, HC33, HC32}{numeric: Healthy controls (continued)}
#'   \item{HC31, HC29, HC25, HC23, HC22}{numeric: Healthy controls (continued)}
#'   \item{B5A42, B5B37, B5C09, B5C06, B5C05}{numeric: Patients potentially affected by the Beckwith-Wiedemann Syndrome (BWS)}
#'   \item{B5C41, B5C04, B5C31, B5C03, B5C28}{numeric: Patients potentially affected by BWS (continued)}
#'   \item{B5C02, B5C19, B5C17, B5C14, B5C13}{numeric: Patients potentially affected by BWS (continued)}
#'   \item{B5C12, B5C01, B5B38}{numeric: Patients potentially affected by BWS (continued)}
#' }
#'
#' @keywords Beckwith-Wiedemann syndrome dataset
#'
#' @references \insertRef{Ochoa2019}{MethylCal}
#'
#' @examples
#' \dontrun{
#'    # Load the BWS dataset and inspect the first rows.
#'  
#'    data(BWS_data)
#'    head(BWS_data)
#' }
"BWS_data"
