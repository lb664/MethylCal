#' Moskalev correction of case/controls samples
#'
#' Correction of case/controls samples by using Moskalev's cubic polynomial
#' regression.
#'
#' @param data Formatted input data frame obtained from the function
#' \code{\link{Formatting}}.
#' @param Target Name of the target DMR/CpG island/gene to be visualised.
#' @param n_Control Number of controls samples.
#' @param n_Case Number of case samples.
#' @param level_Control Level of significance of the differential
#' methylation test between case and control corrected samples. Default
#' value is \code{level_Control = 0.9986501} which correspond to z = 3
#' quantile for normally distributed data.
#' @param opt_BoxPlot Boxplot option: If \code{opt_BoxPlot = 0} (default), 
#' BoxPlot is centered on the median, whereas if \code{opt_BoxPlot = 1} 
#' it is centered on the mean.
#' @param dir In Unix-specific OS, user-specified directory where the
#' plots in \code{\link[grDevices]{pdf}} format are saved. If the directory
#' is not specified, figures are saved in the current working directory.
#' @param printing If \code{printing = TRUE} (default), the corrected 
#' methylation levels of the case/control samples using Moskalev's 
#' cubic polynomial regression (Moskalev et al., 2011) are printed on
#' the screen.
#' @param plotting If \code{plotting = TRUE} (default), the corrected 
#' methylation levels for the control samples as well as its (1-\code{level_Control)}\%
#' confidence interval are depicted. BoxPlot of the corrected methylation
#' levels are also shown for each control sample. If cases are included,
#' a second figure presents a BoxPlot for each case sample as well as
#' the controls' (1-\code{level_Control)}\% confidence interval and 
#' the hyper- and hypo-methylated cases (top red triangles).
#' @param cex_par Number indicating the amount by which plotting text
#' and symbols should be scaled relative to the default (\code{cex_par = 1}).
#'
#' @keywords Case and control samples, Moskalev calibration
#'
#' @import INLA
#' @import lattice
#' @import latticeExtra
#'
#' @export
#'
#' @return This function returns the corrected methylation level for
#' the control and (if selected) case samples using Moskalev's cubic
#' polynomial regression. Based a parametric t-test at (1-\code{level_Control)}\%,
#' hyper- and hypo-methylated cases are also flagged.
#'
#' In Unix-specific OS, figures are saved in the current directory,
#' unless otherwise specified by the user, in \code{\link{pdf}} format.
#' In Windows OS, figures are printed on the screen.
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
#' data = Formatting(BWS_data, AMP = AMP, n_Control = 15)
#' corr_data = MoskalevCorrection(data, Target = "KCNQ1OT1", n_Control = 15)
#'
#' data(Celiac_data)
#' AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#' data = Formatting(Celiac_data, AMP = AMP, n_Control = 13, n_Case = (2 * 17))
#' corr_data = MoskalevCorrection(data, Target = "NFKBIA", n_Control = 13, n_Case = (2 * 17), 
#' opt_BoxPlot = 1)


MoskalevCorrection = function(data, Target = NULL, n_Control = 0, n_Case = 0, level_Control = 0.9986501, opt_BoxPlot = 0, dir = NULL, printing = TRUE, plotting = TRUE, cex_par = 1.25)
{
    
    if (length(unique(data$Target)) > 1)
    {
        if (is.null(Target))
        {
            stop("Please provide the name of the Target for the analysis")
        } else if (sum(data$Target == Target) == 0) {
            stop("Target name not recognised")
        }
    } else {
        Target = unique(data$Target)
    }
    
    f_Moskalev = function(x, y, par)
    {
        m = par[4] * x^3 + par[3] * x^2 + par[2] * x^1 + par[1] * x^0 - y
        m = m ^2
        return(m)
    }
    
    plotsegraph = function(loc, mean, min, max, wiskwidth, lwd = 1, col = "back") {
        w = wiskwidth /2
        graphics::segments(x0 = loc, x1 = loc, y0 = min, y1 = max, lwd = lwd, col = col)
        graphics::segments(x0 = loc - w, x1 = loc + w, y0 = min, y1 = min, lwd = lwd, col = col)
        graphics::segments(x0 = loc - w, x1 = loc + w, y0 = max, y1 = max, lwd = lwd, col = col)
    }
    
    Sys_info = Sys.info()
    
    if (is.null(dir))
    {
        dir = getwd()
    }
    
    cat(rep("\n", 1))
    print(paste("Target", Target, sep = " "))
    cat(rep("\n", 1))
    
    level = 0.95
    n_predict = 200
    optimize_interval = c(-10, 110)
    
    idx = which(data$Target == Target)
    data = data[idx, ]
    col_names = colnames(data)
    data$x1 = data$x ^1
    data$x2 = data$x ^2
    data$x3 = data$x ^3
    
    data = data[!is.na(data$CpG_pos), ]
    n_data = nrow(data)
    p_data = ncol(data)
    n_CpG = length(unique(data$CpG_pos))
    n_CpG_plot = min(6, n_CpG)
    n_AMP = length(unique(data$AMP_group_label))
    AMP = sort(unique(data$x))
    AMP_level = as.character(data$AMP_group_label)[1 : n_AMP]
    n_predict = n_predict + 1
    x_pred = seq(from = 0, to = 100, length.out = n_predict)
    q = stats::qnorm(level + (1 - level) /2)
    q_Control = stats::qnorm(level_Control + (1 - level_Control) /2)
    
    if (n_Control == 0 & n_Case == 0) 
    {
        stop("No case or control samples specified")
    } else {
        if (all(is.na(col_names[(10 : (9 + n_Control))])))
        {
            stop("Input data contain neither case nor control samples")
        } else {
            Control = colnames(data)[(10 : (9 + n_Control))]
            Control_ix = sort(Control, index.return = TRUE)$ix
            if (n_Case > 0)
            {
                Case = colnames(data)[(10 + n_Control) : (9 + n_Control + n_Case)]
                Case_ix = sort(Case, index.return = TRUE)$ix
            } else {
                Case = NULL
            }
        }
    }
    
    Control_Case = c(Control, Case)
    n_Control = length(Control)
    n_Case = length(Case)
    n_Control_Case = length(Control_Case)
    
    CpG_group = rep(NA, n_CpG * n_predict)
    c_lm = rep(NA, 4)
    std_c_lm = rep(NA, 4)
    x_pred_plot = rep(NA, n_CpG)
    fitted_c_lm = rep(NA, n_CpG * n_AMP)
    predict_c_lm = rep(NA, n_CpG * n_predict)
    predict_c_lm_CI_low = rep(NA, n_CpG * n_AMP)
    predict_c_lm_CI_high = rep(NA, n_CpG * n_AMP)
    
    for (s in 1 : n_CpG)
    {
        CpG_idx = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        
        if (!all(is.na(data_CpG$y)))
        {
            fit_lm = stats::lm(y ~ x1 + x2 + x3, data = data_CpG)
            c_lm = stats::coef(summary(fit_lm))[, 1]
            std_c_lm = stats::coef(summary(fit_lm))[, 2]
            
            fitted_c_lm[(1 : n_AMP) + n_AMP * (s - 1)] = stats::fitted(fit_lm)
            x_pred_plot[(1 : n_predict) + n_predict * (s - 1)] = x_pred
            predict_c_lm[(1 : n_predict) + n_predict * (s - 1)] = c_lm[1] + c_lm[2] * x_pred + c_lm[3] * x_pred ^2 + c_lm[4] * x_pred ^3
            # predict_c_lm_CI_low[(1 : n_AMP) + n_AMP * (s - 1)] = (c_lm[1] - q * std_c_lm[1]) + (c_lm[2] - q * std_c_lm[2]) * data_CpG$x + (c_lm[3] - q * std_c_lm[3]) * data_CpG$x ^2 + (c_lm[4] - q * std_c_lm[4]) * data_CpG$x ^3
            # predict_c_lm_CI_high[(1 : n_AMP) + n_AMP * (s - 1)] = (c_lm[1] + q * std_c_lm[1]) + (c_lm[2] + q * std_c_lm[2]) * data_CpG$x + (c_lm[3] + q * std_c_lm[3]) * data_CpG$x ^2 + (c_lm[4] + q * std_c_lm[4]) * data_CpG$x ^3
            predict_c_lm_CI_low[(1 : n_AMP) + n_AMP * (s - 1)] = stats::predict(fit_lm, data_CpG, interval = "confidence", level = level)[, 2]
            predict_c_lm_CI_high[(1 : n_AMP) + n_AMP * (s - 1)] = stats::predict(fit_lm, data_CpG, interval = "confidence", level = level)[, 3]
        }
    
    }
    
    Moskalev_colnames = rep(NA, n_Control_Case)
    
    for (Control_Case_idx in 1 : n_Control_Case)
    {
        Moskalev_colnames[Control_Case_idx] = paste(Control_Case[Control_Case_idx], "Moskalev_corrected", sep = "_")
    }
    
    data_Control_Case = matrix(data = NA, nrow = n_data, ncol = n_Control_Case)
    colnames(data_Control_Case) = Moskalev_colnames
    data = cbind(data, data_Control_Case)
    
    if (printing == TRUE)
    {
        print("Corrected case/control samples:")
    }
    
    for (i in 1 : n_Control_Case)
    {
        Control_Case_name = Control_Case[i]
        Control_Case_idx = which(colnames(data) == Control_Case_name)
        data_colnames_idx = which(colnames(data) == Moskalev_colnames[i])
        
        if (printing == TRUE)
        {
            print(paste("Corrected Target:", Target, "-", "Case Control:", Control_Case_name, sep = " "))
        }
        
        for (s in 1 : n_CpG)
        {
            CpG_idx = which(data$CpG_group == s)
            data_CpG = data[CpG_idx, ]
            y_correct = NA
            
            if (!all(is.na(data_CpG$y)))
            {
                if (is.na(data[CpG_idx, Control_Case_idx][1])) {
                    y_correct = NA
                    data[CpG_idx[1], data_colnames_idx] = y_correct
                } else {
                    fit_lm = stats::lm(y ~ x1 + x2 + x3, data = data_CpG)
                    par = stats::coef(summary(fit_lm))[, 1]
                    y_correct = stats::optimize(f_Moskalev, y = data_CpG[1, Control_Case_idx], par = par, interval = optimize_interval)$minimum
                    y_correct[y_correct < 0] = 0
                    y_correct[y_correct > 100] = 100
                    data[CpG_idx[1], data_colnames_idx] = y_correct
                }
            }
            
            if (printing == TRUE)
            {
                print(paste("Original:", round(data[CpG_idx[1], Control_Case_idx], 2), "-", "Moskalev corrected:", round(y_correct, 2), sep = " "))
            }
        
        }
    
    }
    
    if (printing == TRUE)
    {
        cat(rep("\n", 1))
    }
    
    if (plotting == TRUE)
    {
        title_name = paste(levels(droplevels(data$Target[1])), "Moskalev etal. 2011b: Controls", sep = " - ")
        
        original_data_idx = match(Control_Case_name, colnames(data))
        Moskalev_data_idx = match(Moskalev_colnames, colnames(data))
        original_data = data[seq(1, n_data, n_AMP), original_data_idx]
        Moskalev_data = data[seq(1, n_data, n_AMP), Moskalev_data_idx]
        boxplot_out_Moskalev_Control = graphics::boxplot(Moskalev_data[, Control_ix], data = data, plot = FALSE)
        m = boxplot_out_Moskalev_Control$stats[3, ]
        min = boxplot_out_Moskalev_Control$stats[1, ]
        max = boxplot_out_Moskalev_Control$stats[5, ]
        
        if (opt_BoxPlot == 1)
        {
            m = colMeans(Moskalev_data[, Control_ix], na.rm = TRUE)
            names(m) = NULL 
        }
        
        Moskalev_mean = mean(m, na.rm = TRUE)
        Moskalev_sd = stats::sd(m, na.rm = TRUE)
        
        if (Sys_info[[1]] == "Windows")
        {
            grDevices::dev.new()
        } else {
            name_fig = paste(paste(dir, "/", "MoskalevCorrection_", levels(droplevels(data$Target[1])), "_Fig1", sep = ""), "pdf", sep = ".")
            grDevices::pdf(file = name_fig)
        }
        
        graphics::plot(c(1 : n_Control), m, type = "l", lty = 2, lwd = 2, col = "white", yaxt = "n", xaxt = "n",
        xlim = c(1 - 0.25, n_Control + 0.25), ylim = c(-1, 101),
        cex.axis = cex_par, cex.lab = cex_par, xlab = "", ylab = "Corrected methylation degree", main = title_name)
        llim = max(0, Moskalev_mean - q_Control * Moskalev_sd)
        ulim = min(100, Moskalev_mean + q_Control * Moskalev_sd)
        graphics::lines(c(1 - 0.25, n_Control + 0.25), c(llim, llim), type = "l", lty = 6, lwd = 2, col = "grey", cex.lab = cex_par)
        graphics::lines(c(1 - 0.25, n_Control + 0.25), c(ulim, ulim), type = "l", lty = 6, lwd = 2, col = "grey", cex.lab = cex_par)
        graphics::axis(1, at = c(1 : n_Control), labels = c(Control[Control_ix]), las = 2, cex.axis = cex_par)
        if (n_AMP <= 5)
        {
            graphics::axis(2, at = unique(data$x), cex.axis = cex_par)
        } else {
            graphics::axis(2, at = unique(data$x), labels = FALSE, cex.axis = cex_par)
            x_pos = graphics::par()$usr[1] - 0.05 * (graphics::par()$usr[2] - graphics::par()$usr[1])
            graphics::text(x = x_pos, y = unique(data$x), labels = unique(data$x), srt = 0, adj = 1, xpd = TRUE, cex = cex_par)
        }
        plot.errbars = plotsegraph(c(1 : n_Control), m, min, max, 0.25, col = "black")
        graphics::points(c(1 : n_Control), m, type = "p", pch = 21, lwd = 1.5, cex = 1.05, bg = "white", col = "black")
        
        if (Sys_info[[1]] != "Windows")
        {
            grDevices::dev.off()
        }
        
        if (n_Case > 0)
        {
            title_name = paste(levels(droplevels(data$Target[1])), "Moskalev etal. 2011b: Patients", sep = " - ")
            
            boxplot_out_Moskalev_Case = graphics::boxplot(Moskalev_data[, Case_ix + n_Control], data = data, plot = FALSE)
            m = boxplot_out_Moskalev_Case$stats[3, ]
            min = boxplot_out_Moskalev_Case$stats[1, ]
            max = boxplot_out_Moskalev_Case$stats[5, ]
            
            if (opt_BoxPlot == 1)
            {
                m = colMeans(Moskalev_data[, Case_ix + n_Control], na.rm = TRUE)
                names(m) = NULL  
            }
            
            if (Sys_info[[1]] == "Windows")
            {
                grDevices::dev.new()
            } else {
                name_fig = paste(paste(dir, "/", "MoskalevCorrection_", levels(droplevels(data$Target[1])), "_Fig2", sep = ""), "pdf", sep = ".")
                grDevices::pdf(file = name_fig)
            }
            
            graphics::plot(c(1 : n_Case), m, type = "l", lty = 2, lwd = 2, col = "white", yaxt = "n", xaxt = "n",
            xlim = c(1 - 0.25, n_Case + 0.25), ylim = c(-1, 101),
            cex.axis = cex_par, cex.lab = cex_par, xlab = "", ylab = "Corrected methylation degree", main = title_name)
            llim = max(0, Moskalev_mean - q_Control * Moskalev_sd)
            ulim = min(100, Moskalev_mean + q_Control * Moskalev_sd)
            graphics::lines(c(1 - 0.25, n_Case + 0.25), c(llim, llim), type = "l", lty = 6, lwd = 2, col = "grey", cex.lab = cex_par)
            graphics::lines(c(1 - 0.25, n_Case + 0.25), c(ulim, ulim), type = "l", lty = 6, lwd = 2, col = "grey", cex.lab = cex_par)
            graphics::axis(1, at = c(1 : n_Case), labels = c(Case[Case_ix]), las = 2, cex.axis = cex_par * (18 / n_Case))
            if (n_AMP <= 5)
            {
                graphics::axis(2, at = unique(data$x), cex.axis = cex_par)
            } else {
                graphics::axis(2, at = unique(data$x), labels = FALSE, cex.axis = cex_par)
                x_pos = graphics::par()$usr[1] - 0.05 * (graphics::par()$usr[2] - graphics::par()$usr[1])
                graphics::text(x = x_pos, y = unique(data$x), labels = unique(data$x), srt = 0, adj = 1, xpd = TRUE, cex = cex_par)
            }
            plot.errbars = plotsegraph(c(1 : n_Case), m, min, max, 0.25, col = "black")
            graphics::points(c(1 : n_Case), m, type = "p", pch = 21, lwd = 1.5, cex = 1.05, bg = "white", col = "black")
            
            Moskalev_idx = (llim > m | m > ulim)
            n_Moskalev_idx = sum(Moskalev_idx, na.rm = TRUE)
            if (n_Moskalev_idx > 0)
            {
                graphics::points(which(Moskalev_idx), rep(100, n_Moskalev_idx), type = "p", pch = 25, lwd = 1.25, cex = 1.25, bg = "red", col = "red")
            }
            
            if (Sys_info[[1]] != "Windows")
            {
                grDevices::dev.off()
            }
            
            Moskalev_idx = which(Moskalev_idx)
            if ((printing == TRUE) & length(Moskalev_idx > 0))
            {
                print("Moskalev hyper/hypomethylated samples:")
                print(Case[Case_ix][Moskalev_idx])
                print(round(m[Moskalev_idx], 2))
                cat(rep("\n", 1))
            }
        }
    }
    
    return(data)
}
