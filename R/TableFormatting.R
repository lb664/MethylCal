#' Exporting results utility function
#'
#' Reformatting calibration output.
#'
#' @param data Input data frame obtained either using MethylCal's
#' Bayesian calibration method (Ochoa et al., 2019) or Moskalev's cubic
#' polynomial calibration (Moskalev et al., 2011).
#' @param n_Control Number of controls samples.
#' @param n_Case Number of case samples.
#' @param dir User-specified directory where the calibrated methylation 
#' level are saved. If the directory is not specified, table is
#' saved in the current working directory.
#' @param writing If \code{writing = TRUE} (default), the file of the 
#' corrected methylation levels of the case/control samples are saved
#' in \code{dir}.
#'
#' @keywords Reformatting, calibration input data frame
#'
#' @import INLA
#' @import lattice
#' @import latticeExtra
#'
#' @export
#'
#' @return Short formatted output data frame. It includes the name of 
#' the DRM/CpG island/gene ("Target"), chromosome ("Chrom"), CpG position
#' ("CpG_group"), uncalibrated methylation levels as well as the
#' corresponding levels based either on MethylCal's Bayesian calibration
#' method or Moskalev's cubic polynomial calibration.
#'
#' @references Ochoa E, Zuber V, Fernandez-Jimenez N, Bilbao JR, Clark 
#' GR, Maher ER and Bottolo L. MethylCal: Bayesian calibration of methylation
#' levels. Submitted. 2019.
#' @references Moskalev EA, Zavgorodnij MG, Majorova SP, Vorobjev IA, 
#' Jandaghi P, Bure IV, Hoheisel JD. Correction of PCR-bias in quantitative
#' DNA methylation studies by means of cubic polynomial regression. 
#' Nucleic Acids Res. 2011; 39(11):e77. (\href{https://www.ncbi.nlm.nih.gov/pubmed/21486748}{PubMed})
#'
#' @examples
#' data(BWS_data)
#' AMP = c(0, 25, 50, 75, 100)
#' data = Formatting(BWS_data, AMP = AMP, n_Control = 15, n_Case = 18)
#' corr_data = MoskalevCorrection(data, Target = "KCNQ1OT1", n_Control = 15, n_Case = 18)
#' output_data = TableFormatting(corr_data, n_Control = 15, n_Case = 18)
#'
#' data(Celiac_data)
#' AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#' data = Formatting(Celiac_data, AMP = AMP, n_Control = 13, n_Case = (2 * 17))
#' corr_data = MethylCalCorrection(data, Target = "NFKBIA", n_Control = 13, n_Case = (2 * 17))
#' output_data = TableFormatting(corr_data, n_Control = 13, n_Case = (2 * 17), writing = FALSE)


TableFormatting = function(data, n_Control = 0, n_Case = 0, dir = NULL, writing = TRUE)
{
    Sys_info = Sys.info()
    
    if (is.null(dir))
    {
        dir = getwd()
    }
        
    if (n_Control == 0 & n_Case == 0) 
    {
        stop("No Controls or Cases specified")
    } else {
        Control = colnames(data)[(10 : (9 + n_Control))]
        Control_ix = sort(Control, index.return = T)$ix
        if (n_Case > 0)
        {
            Case = colnames(data)[(10 + n_Control) : (9 + n_Control + n_Case)]
            Case_ix = sort(Case, index.return = T)$ix
        } else {
            Case = NULL
        }
    }
    
    data = data[!is.na(data$CpG_pos), ]
    n_data = nrow(data)
    p_data = ncol(data)
    n_CpG = length(unique(data$CpG_pos))
    n_CpG_plot = min(6, n_CpG)
    n_AMP = length(unique(data$AMP_group_label))
    AMP = sort(unique(data$x))
    AMP_level = as.character(data$AMP_group_label)[1 : n_AMP]
    
    Control_Case = c(Control, Case)
    n_Control = length(Control)
    n_Case = length(Case)
    n_Control_Case = length(Control_Case)
    data_Control_Case_idx = 10 : (10 + n_Control_Case - 1)
    data_corrected_idx = (p_data - n_Control_Case + 1) : p_data
    data_Control_Case = round(data[seq(1, n_data, n_AMP), data_Control_Case_idx], 3)
    data_corrected = round(data[seq(1, n_data, n_AMP), data_corrected_idx], 3)
    
    data_tmp = matrix(data = NA, nrow = n_CpG, ncol = 3)
    colnames(data_tmp) = c("Target", "Chrom", "CpG_position")
    data_tmp[, 1] = levels(droplevels(data$Target[1]))
    data_tmp[, 2] = data$Chrom[1]
    data_tmp[, 3] = unique(data$CpG_pos)
    data_tmp = cbind(cbind(data_tmp, data_Control_Case), data_corrected)
    
    if (writing == T)
    {
        if (any(colnames(data) == "xCpG_group"))
        {
            name_tab = paste(paste(dir, "/", "MethylCalCorrection_", levels(droplevels(data$Target[1])), sep = ""), "txt", sep = ".")
            
        } else {
            name_tab = paste(paste(dir, "/", "MoskalevCorrection_", levels(droplevels(data$Target[1])), sep = ""), "txt", sep = ".")
        }

    utils::write.table(data_tmp, file = name_tab, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t") 
    }
    
    return(data_tmp)
}
