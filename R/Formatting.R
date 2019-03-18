#' Reformatting utility function
#'
#' Reformatting data frame for MethylCal analysis.
#'
#' @param data Input data frame formatted. It should contain the
#' following columns: "Target" (name of the target DMR/CpG island/gene 
#' name), "Chrom" (Chromosome), "CpG pos" (CpG position) along with
#' the columns for the measurements of the apparent level of 
#' methylation for each Actual Methylation Percentage (AMP) and for 
#' each CpG. Missing values are allowed and should be marked as NA.
#' This utility function can process several Targets if included in
#' the input data frame.
#' @param AMP Numeric vector of Actual Methylation Percentages (AMPs).
#' @param n_Control Number of controls samples.
#' @param n_Case Number of case samples.
#'
#' @keywords Reformatting, input data frame
#'
#' @import INLA
#' @import lattice
#' @import latticeExtra
#'
#' @export
#'
#' @return Long formatted output data frame. It includes the name of 
#' the DRM/CpG island/gene ("Target"), chromosome ("Chrom"), CpG position
#' ("CpG_group"), at which CpG the measurement has been taken ("CpG_group")
#' within each DRM/CpG island/gene, a label attached to each CpG 
#' ("CpG_group_label"), at which AMP the measurement has been taken 
#' ("AMP_group"), a label attached to each CpG ("CpG_group_label")
#' as well as the apparent ("y") and actual ("x") methylation levels.
#'
#' If control/cases samples are included, methylation levels are stored
#' in the long formatted data frame at each row corresponding to 
#' CpG_group = 1 and AMP_group = 1.
#'
#' @references Ochoa E, Zuber V, Fernandez-Jimenez N, Bilbao JR, Clark 
#' GR, Maher ER, Bottolo L. MethylCal: Bayesian calibration of methylation
#' levels. Submitted. 2019.
#'
#' @examples
#' data(BWS_data)
#' AMP = c(0, 25, 50, 75, 100)
#' data = Formatting(BWS_data, AMP = AMP)
#' data = Formatting(BWS_data, AMP = AMP, n_Control = 15)
#' data = Formatting(BWS_data, AMP = AMP, n_Control = 15, n_Case = 18)
#'
#' data(Celiac_data)
#' AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#' data = Formatting(Celiac_data, AMP = AMP)
#' data = Formatting(Celiac_data, AMP = AMP, n_Control = 13)
#' data = Formatting(Celiac_data, AMP = AMP, n_Control = 13, n_Case = (2 * 17))


Formatting = function(data, AMP = NULL, n_Control = 0, n_Case = 0) 
{
    
    if (is.null(AMP) & ((n_Control + n_Case) == 0))
    {
        stop("Please provide either the vector of AMPs or the number of cases and controls")
    }
    
    Sys_info = Sys.info()
    
    n_data = nrow(data)
    p_data = ncol(data)
    col_names = colnames(data)
    n_CpG = n_data
    
    if (is.null(AMP))
    {
        AMP = as.vector(as.numeric(unlist(strsplit(col_names[4 : (p_data - n_Control - n_Case)], "%"))))
    } else {
        idx = which(!is.na(AMP))
        AMP = AMP[idx]
    }
    
    n_AMP = length(AMP)
    
    data_tmp = matrix(data = NA, nrow = n_CpG * n_AMP, ncol = (9 + n_Control + n_Case))
    data_tmp[, 1] = rep(data[, 1], each = n_AMP)
    data_tmp[, 2] = rep(data[, 2], each = n_AMP)
    data_tmp[, 3] = rep(data[, 3], each = n_AMP)
    Target = unique(data[, 1])
    n_Target = length(Target)
    
    for (t in 1 : n_Target) 
    {
        if (t == 1)
        {
            tmp = rep(seq(1 : sum(data[, 1] == Target[t])), each = n_AMP)
        } else {
            tmp = c(tmp, rep(seq(1 : sum(data[, 1] == Target[t])), each = n_AMP))
        }
    }
    
    data_tmp[, 4] = tmp
    data_tmp[, 5] = paste("CpG", tmp, sep = "_")
    data_tmp[, 6] = rep(1 : n_AMP, times = n_CpG)
    data_tmp[, 7] = rep(paste(AMP, "%", sep = ""), times = n_CpG)
    data_tmp[, 8] = as.vector(t(as.matrix(data[, 4 : (4 + n_AMP - 1)])))
    data_tmp[, 9] = rep(AMP, n_CpG)
    
    if (n_Control > 0 | n_Case > 0) 
    { 
        for (j in 1 : (n_Control + n_Case))
        {
            data_tmp[seq(1, (n_CpG * n_AMP), by = n_AMP), (9 + j)] = data[, (3 + n_AMP + j)]
        }
    }
    
    data_tmp = data.frame(data_tmp, stringsAsFactors = TRUE)
    factorToNumeric = function(f) as.numeric(levels(f))[as.integer(f)]
    
    if (n_Control == 0 & n_Case == 0) 
    {
        idx = c(2, 3, 4, 6, 8, 9)
    } else {
        idx = c(2, 3, 4, 6, 8, 9, (10 : (9 + n_Control + n_Case)))
    }
    
    data_tmp[, idx] = lapply(data_tmp[, idx], factorToNumeric)
    
    if (n_Control == 0 & n_Case == 0) 
    {
        colnames(data_tmp) = c("Target", "Chrom", "CpG_pos", "CpG_group", "CpG_group_label", "AMP_group", "AMP_group_label", "y", "x")
    } else {
        colnames(data_tmp) = c("Target", "Chrom", "CpG_pos", "CpG_group", "CpG_group_label", "AMP_group", "AMP_group_label", "y", "x", col_names[(4 + n_AMP) : ((3 + n_AMP) + n_Control + n_Case)])
    }
    
    return(data_tmp)
}
