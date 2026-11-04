#' Export results utility function
#'
#' Reformatting calibration output.
#'
#' @param data Input data frame obtained either using MethylCal's Bayesian
#'   calibration method \insertCite{Ochoa2019}{MethylCal} or Moskalev's cubic
#'   polynomial calibration \insertCite{Moskalev2011}{MethylCal}.
#' @param n_Control Number of control samples.
#' @param n_Case Number of case samples.
#' @param dir User-specified directory where the calibrated methylation levels
#'   are saved. If the directory is not specified, the table is saved in the
#'   current working directory.
#' @param write_out If \code{write_out = TRUE} (default), the file of the
#'   corrected methylation levels of the case/control samples are saved in
#'   \code{dir}.
#' @param writing Deprecated alias of \code{write_out}, kept for backward
#'   compatibility with MethylCal version 0.1.0.
#' @keywords Reformatting, calibration input data frame
#'
#' @export
#'
#' @return Short formatted output data frame. It includes the name of the
#'   DMR/CpG island/gene ("Target"), chromosome ("Chrom"), CpG position
#'   ("CpG_group"), uncalibrated methylation levels as well as the
#'   corresponding levels based either on MethylCal's Bayesian calibration
#'   method or Moskalev's cubic polynomial calibration.
#'
#' @references \insertRef{Ochoa2019}{MethylCal}
#' @references \insertRef{Moskalev2011}{MethylCal}
#'
#' @examples
#' \dontrun{
#'    # Export corrected methylation levels to a tab-delimited text file. The first
#'    # example uses Moskalev correction on BWS; the second uses MethylCal
#'    # correction on the Celiac dataset with file writing disabled.
#'  
#'    data(BWS_data)
#'    AMP = c(0, 25, 50, 75, 100)
#'    data = Formatting(BWS_data, AMP = AMP, n_Control = 15, n_Case = 18)
#'    corr_data = MoskalevCorrection(data, Target = "KCNQ1OT1", n_Control = 15,
#'                                   n_Case = 18)
#'    output_data = TableFormatting(corr_data, n_Control = 15, n_Case = 18)
#'  
#'    data(Celiac_data)
#'    AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#'    data = Formatting(Celiac_data, AMP = AMP, n_Control = 13, n_Case = (2 * 17))
#'    corr_data = MethylCalCorrection(data, Target = "NFKBIA", n_Control = 13,
#'                                    n_Case = (2 * 17))
#'    output_data = TableFormatting(corr_data, n_Control = 13, n_Case = (2 * 17),
#'                                  write_out = FALSE)
#'    
#' }


TableFormatting = function(data, n_Control = 0, n_Case = 0, dir = NULL, write_out = TRUE,
                           writing = NULL)
{
    if (!is.null(writing))
    {
        warning("Argument 'writing' is deprecated; use 'write_out' instead.", call. = FALSE)
        write_out = writing
    }

    if (is.null(data))
    {
        warning("TableFormatting: input data is NULL ",
                "(correction step may have failed). Returning NULL.")
        return(invisible(NULL))
    }
    Sys_info = Sys.info()
    if (is.null(dir)) dir = getwd()
    dir = sub("/+$", "", dir)   # strip trailing slashes to avoid double // in filenames
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    target_safe = gsub("[/\\:*?\"<>|]", "_", as.character(data$Target[1]))

    if (n_Control == 0 & n_Case == 0)
    {
        stop("No Controls or Cases specified")
    } else {
        Control  = colnames(data)[(10 : (9 + n_Control))]
        Control_ix = sort(Control, index.return = TRUE)$ix
        if (n_Case > 0)
        {
            Case = colnames(data)[(10 + n_Control) : (9 + n_Control + n_Case)]
            Case_ix = sort(Case, index.return = TRUE)$ix
        } else {
            Case = NULL
        }
    }

    data = data[!is.na(data$CpG_pos), ]
    n_data = nrow(data)
    p_data = ncol(data)
    n_CpG  = length(unique(data$CpG_pos))
    n_AMP  = length(unique(data$AMP_group_label))
    AMP    = sort(unique(data$x))

    Control_Case     = c(Control, Case)
    n_Control        = length(Control)
    n_Case           = length(Case)
    n_Control_Case   = length(Control_Case)
    data_Control_Case_idx = 10 : (10 + n_Control_Case - 1)
    data_corrected_idx    = (p_data - n_Control_Case + 1) : p_data
    data_Control_Case     = round(data[seq(1, n_data, n_AMP), data_Control_Case_idx], 3)
    data_corrected        = round(data[seq(1, n_data, n_AMP), data_corrected_idx], 3)

    data_tmp = matrix(data = NA, nrow = n_CpG, ncol = 3)
    colnames(data_tmp) = c("Target", "Chrom", "CpG_position")
    data_tmp[, 1] = target_safe
    data_tmp[, 2] = data$Chrom[1]
    data_tmp[, 3] = unique(data$CpG_pos)
    data_tmp = cbind(cbind(data_tmp, data_Control_Case), data_corrected)

    if (write_out == TRUE)
    {
        if (any(colnames(data) == "xCpG_group"))
        {
            name_tab = paste0(dir, "/MethylCalCorrection_",
                              target_safe, ".txt")
        } else {
            name_tab = paste0(dir, "/MoskalevCorrection_",
                              target_safe, ".txt")
        }
        utils::write.table(data_tmp, file = name_tab, quote = FALSE,
                           row.names = FALSE, col.names = TRUE, sep = "\t")
    }

    return(data_tmp)
}
