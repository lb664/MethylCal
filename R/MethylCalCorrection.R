#' MethylCal correction of case/controls samples
#'
#' Correction of case/controls samples using MethylCal calibration
#'
#' @param data Formatted input data frame obtained from the function
#' \code{\link{Formatting}}.
#' @param Target Name of the target DMR/CpG island/gene to be visualised.
#' @param prior Prior distribution set-up for the random effects and 
#' the Latent Gaussian Field (Rue et al., 2009). Three different priors
#' are implemented:
#' \itemize{
#'   \item \code{LG}: log-Gamma prior is the default prior with 
#'   \code{a = 1} and \code{b = 0.1} parametrization;
#'   \item \code{PC}: Penalized Complexity prior (Simpson et al., 
#'   2017);
#'   \item \code{HC}: Half-Cauchy prior (Wang X et al., 2018).
#' }
#' @param n_Control Number of controls samples.
#' @param n_Case Number of case samples.
#' @param level_Control Level of significance of the differential
#' methylation test between case and control corrected samples. Default
#' value is \code{level_Control = 0.9986501} which correspond to z = 3
#' quantile for normally distributed data.
#' @param opt_BoxPlot Boxplot option: If \code{opt_BoxPlot = 0} (default), 
#' BoxPlot is centered around the median, whereas if \code{opt_BoxPlot = 1} 
#' it is centered around the mean.
#' @param dir In Unix-specific OS, user-specified directory where the
#' plots in \code{\link[grDevices]{pdf}} format are saved. If the directory
#' is not specified, figures are saved in the current working directory.
#' @param printing If \code{printing = TRUE} (default), the corrected 
#' methylation levels of the case/control samples using MethylCal calibration
#' (Ochoa et al., 2019) are printed on the screen.
#' @param plotting If \code{plotting = TRUE} (default), the corrected 
#' methylation levels for the control samples as well as its (1-\code{level_Control)}\%
#' confidence interval are depicted. BoxPlot of the corrected methylation
#' levels are also shown for each control sample. If cases are
#' included, a second figure presents a BoxPlot for each case sample
#' as well as the controls' (1-\code{level_Control)}\% confidence 
#' interval and the hyper- and hypo-methylated cases (top red triangles).
#' @param cex_par Number indicating the amount by which plotting text
#' and symbols should be scaled relative to the default (\code{cex_par = 1}).
#'
#' @keywords Case and control samples, MethylCal calibration
#'
#' @import INLA
#' @import lattice
#' @import latticeExtra
#'
#' @export
#'
#' @return This function returns the corrected methylation level for
#' the control and (if selected) case samples using MethylCal calibration.
#' Based a parametric t-test at (1-\code{level_Control)}\%, hyper- and
#' hypo-methylated cases are also flagged.
#'
#' In Unix-specific OS, figures are saved in the current directory,
#' unless otherwise specified by the user, in \code{\link{pdf}} format.
#' In Windows OS, figures are printed on the screen.
#'
#' @references Ochoa E, Zuber V, Fernandez-Jimenez N, Bilbao JR, Clark 
#' GR, Maher ER and Bottolo L. MethylCal: Bayesian calibration of methylation
#' levels. Submitted. 2019.
#' @references Wang X, Ryan YY, Faraway JJ. Bayesian Regression Modeling
#' with INLA. 2018, 1st edition. Chapman and Hall/CRC.
#' @references Simpson S, Rue H, Riebler A, Martins TG, Sorbye SH. Penalising
#' model component complexity: A principled, practical approach to constructing
#' priors. Statist Sci. 2017; 1:1-28.  (\href{https://doi.org/10.1214/16-STS576}{doi})
#' @references Rue H, Martino S, Chopin N. Approximate Bayesian inference
#' for latent Gaussian models by using integrated nested Laplace approximations.
#' J Roy Stat Soc B Met. 2009; 71(2):319-392. (\href{https://doi.org/10.1111/j.1467-9868.2008.00700.x}{doi})
#'
#' @examples
#' data(BWS_data)
#' AMP = c(0, 25, 50, 75, 100)
#' data = Formatting(BWS_data, AMP = AMP, n_Control = 15)
#' corr_data = MethylCalCorrection(data, Target = "KCNQ1OT1", n_Control = 15)
#'
#' data(Celiac_data)
#' AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#' data = Formatting(Celiac_data, AMP = AMP, n_Control = 13, n_Case = (2 * 17))
#' corr_data = MethylCalCorrection(data, Target = "NFKBIA", n_Control = 13, n_Case = (2 * 17), 
#' opt_BoxPlot = 1)


MethylCalCorrection = function(data, Target = NULL, prior = "LG", n_Control = 0, n_Case = 0, level_Control = 0.9986501, opt_BoxPlot = 0, dir = NULL, printing = TRUE, plotting = TRUE, cex_par = 1.25)
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

    f_MethylCal = function(data, CpG_group, fit_MethylCal, Control_Case_name, prior, n_predict = 10 ^4)
    {
        n_data = nrow(data)
        p_data = ncol(data)
        n_CpG = max(data$CpG_group)
        n_AMP = max(data$AMP_group)
        
        CpG_group_idx = which(data$CpG_group == CpG_group)
        CpG_idx = which(data$CpG_group == CpG_group)[1]
        Control_Case_idx = which(colnames(data) == Control_Case_name)
        
        predict_MethylCal = rep(NA, n_predict)
        fit_spline = stats::spline(data$x[CpG_group_idx], fit_MethylCal$summary.fitted$mean[CpG_group_idx], method = "natural", n = n_predict)
        predict_MethylCal = fit_spline$y
        m = fit_spline$y - data[CpG_idx, Control_Case_idx]
        
        if (all(diff(m) >= 0))
        {
            m = m ^2
            idx_min = which.min(m)
            x = seq(0, 100, length.out = n_predict)[idx_min]
            m = (fit_spline$y - data[CpG_idx, Control_Case_idx])[idx_min]
            m = (fit_spline$y)[idx_min]
        } else {
            m = m ^2
            idx_min = sort(m, index.return = TRUE)$ix[1 : (10 ^2)]
            idx_min = idx_min[which.min((seq(0, 100, length.out = n_predict)[idx_min] - data[CpG_idx, Control_Case_idx]) ^2)]
            x = seq(0, 100, length.out = n_predict)[idx_min]
            m = (fit_spline$y - data[CpG_idx, Control_Case_idx])[idx_min]
            print("Non-decreasing function")
        }
         
        return(x)
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
    
    DIC_remove = 1
    level = 0.95
    n_predict = 200
    
    idx = which(data$Target == Target)
    data = data[idx, ]
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
    data$CpG_group = rep(1 : n_CpG, each = n_AMP)
    n_predict = n_predict + 1
    x_pred = seq(from = 0, to = 100, length.out = n_predict)
    CpG_group = rep(NA, n_CpG * n_predict)
    x_pred_plot = rep(NA, n_CpG * n_predict)
    q = stats::qnorm(level + (1 - level) /2)
    q_Control = stats::qnorm(level_Control + (1 - level_Control) /2)
    
    for (s in 1 : n_CpG) 
    {
        CpG_idx = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        x_pred_plot[(1 : n_predict) + n_predict * (s - 1)] = x_pred
    }
    
    data$xCpG_group = data$CpG_group + max(data$CpG_group)
    
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
    
    Control_Case = c(Control, Case)
    n_Control = length(Control)
    n_Case = length(Case)
    n_Control_Case = length(Control_Case)
    
    if (prior == "LG")
    {
        prior = list(prec = list(prior = "loggamma", param = c(1, 0.1)))
    } else if (prior == "LGEB") {
        stop("Empirical Bayes LogGamma prior to be implemented")
    } else if (prior == "PC") {
        sdy = stats::sd(data$y, na.rm = TRUE)
        prior = list(prec = list(prior = "pc.prec", param = c(3 * sdy, 0.01)))
    } else if (prior == "HC") {
        halfcauchy = "expression:
                     lambda = 0.022;
                     precision = exp(log_precision);
                     logdens = -1.5 * log_precision - log(pi * lambda) - log(1 + 1 /(precision * lambda ^2));
                     log_jacobian = log_precision;
                     return(logdens + log_jacobian);"
        prior = list(prec = list(prior = halfcauchy))
    }
    
    formula_MethylCal_1 = y ~ x1 + x2 + x3
    formula_MethylCal_2 = y ~ x1 + x2 + x3 + f(AMP_group, model = "iid", hyper = prior, constr = TRUE)
    formula_MethylCal_3 = y ~ x1 + x2 + x3 + f(AMP_group, model = "iid", hyper = prior, constr = TRUE) + f(CpG_group, model = "iid", hyper = prior, constr = TRUE)
    formula_MethylCal_4 = y ~ x1 + x2 + x3 + f(AMP_group, model = "iid", hyper = prior, constr = TRUE) + f(CpG_pos, model = "rw1", hyper = prior, scale.model = TRUE)
    formula_MethylCal_5 = y ~ x1 + x2 + x3 + f(AMP_group, model = "iid", hyper = prior, constr = TRUE) + f(CpG_group, model = "iid2d", n = 2 * n_CpG, constr = TRUE) + f(xCpG_group, AMP_group, copy = "CpG_group")
    
    fit_MethylCal_1 = inla(formula_MethylCal_1, family = "gaussian", control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE), data = data)
    if (printing == TRUE)
    {
        print("Model 1 ...")
    }
    fit_MethylCal_2 = inla(formula_MethylCal_2, family = "gaussian", control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE), data = data)
    if (printing == TRUE)
    {
        print("Model 1 estimated")
        print("Model 2 ...")
    }
    fit_MethylCal_3 = inla(formula_MethylCal_3, family = "gaussian", control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE), data = data)
    if (printing == TRUE)
    {
        print("Model 2 estimated")
        print("Model 3 ...")
    }
    fit_MethylCal_4 = inla(formula_MethylCal_4, family = "gaussian", control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE), data = data)
    if (printing == TRUE)
    {
        print("Model 3 estimated")
        print("Model 4 ...")
    }
    fit_MethylCal_5 = inla(formula_MethylCal_5, family = "gaussian", control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE), data = data)
    if (printing == TRUE)
    {
        print("Model 4 estimated")
        cat(rep("\n", 1))
    }
    
    mlik_MethylCal = rep(NA, 5)
    mlik_MethylCal[1] = fit_MethylCal_1$mlik[2]
    mlik_MethylCal[2] = fit_MethylCal_2$mlik[2]
    mlik_MethylCal[3] = fit_MethylCal_3$mlik[2]
    mlik_MethylCal[4] = fit_MethylCal_4$mlik[2]
    mlik_MethylCal[5] = fit_MethylCal_5$mlik[2]
    
    DIC_MethylCal = rep(NA, 5)
    DIC_MethylCal[1] = fit_MethylCal_1$dic[[1]]
    DIC_MethylCal[2] = fit_MethylCal_2$dic[[1]]
    DIC_MethylCal[3] = fit_MethylCal_3$dic[[1]]
    DIC_MethylCal[4] = fit_MethylCal_4$dic[[1]]
    DIC_MethylCal[5] = fit_MethylCal_5$dic[[1]]
    
    fitted_MethylCal = matrix(data = NA, nrow = n_CpG * n_AMP, ncol = 5)
    fitted_MethylCal[, 1] = fit_MethylCal_1$summary.fitted.values$mean
    fitted_MethylCal[, 2] = fit_MethylCal_2$summary.fitted.values$mean
    fitted_MethylCal[, 3] = fit_MethylCal_3$summary.fitted.values$mean
    fitted_MethylCal[, 4] = fit_MethylCal_4$summary.fitted.values$mean
    fitted_MethylCal[, 5] = fit_MethylCal_5$summary.fitted.values$mean
    
    fitted_MethylCal[fitted_MethylCal[, 1] < 0, 1] = 0 
    fitted_MethylCal[fitted_MethylCal[, 1] > 100, 1] = 100
    fitted_MethylCal[fitted_MethylCal[, 2] < 0, 2] = 0 
    fitted_MethylCal[fitted_MethylCal[, 2] > 100, 2] = 100
    fitted_MethylCal[fitted_MethylCal[, 3] < 0, 3] = 0 
    fitted_MethylCal[fitted_MethylCal[, 3] > 100, 3] = 100
    fitted_MethylCal[fitted_MethylCal[, 4] < 0, 4] = 0 
    fitted_MethylCal[fitted_MethylCal[, 4] > 100, 4] = 100
    fitted_MethylCal[fitted_MethylCal[, 5] < 0, 5] = 0 
    fitted_MethylCal[fitted_MethylCal[, 5] > 100, 5] = 100
    
    RSS_MethylCal = rep(NA, 5)
    RSS_MethylCal[1] = sum((data$y - fitted_MethylCal[, 1]) ^2, na.rm = TRUE)
    RSS_MethylCal[2] = sum((data$y - fitted_MethylCal[, 2]) ^2, na.rm = TRUE)
    RSS_MethylCal[3] = sum((data$y - fitted_MethylCal[, 3]) ^2, na.rm = TRUE)
    RSS_MethylCal[4] = sum((data$y - fitted_MethylCal[, 4]) ^2, na.rm = TRUE)
    RSS_MethylCal[5] = sum((data$y - fitted_MethylCal[, 5]) ^2, na.rm = TRUE)
    
    DIC_MethylCal_tmp = DIC_MethylCal
    DIC_MethylCal_tmp[DIC_remove] = NA
    MethylCal_idx = which(DIC_MethylCal_tmp == min(DIC_MethylCal_tmp, na.rm = TRUE))
    
    if (printing == TRUE)
    {
        print("Estimated models:")
        print(c("", "Model 1", "Model 2", "Model 3", "Model 4"))
        print(c("mlik", round(mlik_MethylCal, 2)))
        print(c("DIC", round(DIC_MethylCal, 2)))
        print(c("RSS", round(RSS_MethylCal, 2)))
        cat(rep("\n", 1))
    }
    
    if (MethylCal_idx == 1) 
    {
        fit_MethylCal = fit_MethylCal_1
    } else if (MethylCal_idx == 2)
    {
        fit_MethylCal = fit_MethylCal_2
    } else if (MethylCal_idx == 3)
    {
        fit_MethylCal = fit_MethylCal_3
    } else if (MethylCal_idx == 4)
    {
        fit_MethylCal = fit_MethylCal_4
    } else if (MethylCal_idx == 5)
    {
        fit_MethylCal = fit_MethylCal_5
    }
    
    MethylCal_colnames = rep(NA, n_Control_Case)
    
    for (Control_Case_idx in 1 : n_Control_Case)
    {
        MethylCal_colnames[Control_Case_idx] = paste(Control_Case[Control_Case_idx], "MethylCal_corrected", sep = "_")
    }
    
    data_Control_Case = matrix(data = NA, nrow = n_data, ncol = n_Control_Case)
    colnames(data_Control_Case) = MethylCal_colnames
    data = cbind(data, data_Control_Case)
    
    if (printing == TRUE)
    {
        print("Corrected case/control samples:")
    }
    
    for (i in 1 : n_Control_Case)
    {
        Control_Case_name = Control_Case[i]
        Control_Case_idx = which(colnames(data) == Control_Case_name)
        data_colnames_idx = which(colnames(data) == MethylCal_colnames[i])
        
        if (printing == TRUE)
        {
            print(paste("Corrected Target:", Target, "-", "Case Control:", Control_Case_name, sep = " "))
        }
        
        for (s in 1 : n_CpG)
        {
            CpG_idx = which(data$CpG_group == s)
            CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
            data_CpG = data[CpG_idx, ]
            y_correct = NA
            
            if (!all(is.na(data_CpG$y)))
            {
                if (is.na(data[CpG_idx, Control_Case_idx][1]))
                {
                    y_correct = NA
                    data[CpG_idx[1], data_colnames_idx] = y_correct
                } else {
                    y_correct = f_MethylCal(data = data, CpG_group = s, fit_MethylCal = fit_MethylCal, Control_Case_name = Control_Case_name, prior = prior)
                    y_correct[y_correct < 0] = 0
                    y_correct[y_correct > 100] = 100
                    data[CpG_idx[1], data_colnames_idx] = y_correct
                }
            }
            
            if (printing == TRUE)
            {
                print(paste("Original:", round(data[CpG_idx[1], Control_Case_idx], 2), "-", "Methylcal corrected:", round(y_correct, 2), sep = " "))
            }
        
        }
    
    }
    
    if (printing == TRUE)
    {
        cat(rep("\n", 1))
    }
    
    if (plotting == TRUE)
    {
        title_name = paste(levels(droplevels(data$Target[1])), "MethylCal: Controls", sep = " - ")
        
        original_data_idx = match(Control_Case_name, colnames(data))
        MethylCal_data_idx = match(MethylCal_colnames, colnames(data))
        original_data = data[seq(1, n_data, n_AMP), original_data_idx]
        MethylCal_data = data[seq(1, n_data, n_AMP), MethylCal_data_idx]
        boxplot_out_MethylCal_Control = graphics::boxplot(MethylCal_data[, Control_ix], data = data, plot = FALSE)
        m = boxplot_out_MethylCal_Control$stats[3, ]
        min = boxplot_out_MethylCal_Control$stats[1, ]
        max = boxplot_out_MethylCal_Control$stats[5, ]
        
        if (opt_BoxPlot == 1)
        {
            m = colMeans(MethylCal_data[, Control_ix], na.rm = TRUE)
            names(m) = NULL 
        }
        
        MethylCal_mean = mean(m, na.rm = TRUE)
        MethylCal_sd = stats::sd(m, na.rm = TRUE)
        
        if (Sys_info[[1]] == "Windows")
        {
            grDevices::dev.new()
        } else {
            name_fig = paste(paste(dir, "/", "MethylCalCorrection_", levels(droplevels(data$Target[1])), "_Fig1", sep = ""), "pdf", sep = ".")
            grDevices::pdf(file = name_fig)
        }
        
        graphics::plot(c(1 : n_Control), m, type = "l", lty = 2, lwd = 2, col = "white", yaxt = "n", xaxt = "n",
        xlim = c(1 - 0.25, n_Control + 0.25), ylim = c(-1, 101),
        cex.axis = cex_par, cex.lab = cex_par, xlab = "", ylab = "Corrected methylation degree", main = title_name)
        llim = max(0, MethylCal_mean - q_Control * MethylCal_sd)
        ulim = min(100, MethylCal_mean + q_Control * MethylCal_sd)
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
            title_name = paste(levels(droplevels(data$Target[1])), "MethylCal: Patients", sep = " - ")
            
            boxplot_out_MethylCal_Case = graphics::boxplot(MethylCal_data[, Case_ix + n_Control], data = data, plot = FALSE)
            m = boxplot_out_MethylCal_Case$stats[3, ]
            min = boxplot_out_MethylCal_Case$stats[1, ]
            max = boxplot_out_MethylCal_Case$stats[5, ]
            
            if (opt_BoxPlot == 1)
            {
                m = colMeans(MethylCal_data[, Case_ix + n_Control], na.rm = TRUE)
                names(m) = NULL  
            }
            
            if (Sys_info[[1]] == "Windows")
            {
                grDevices::dev.new()
            } else {
                name_fig = paste(paste(dir, "/", "MethylCalCorrection_", levels(droplevels(data$Target[1])), "_Fig2", sep = ""), "pdf", sep = ".")
                grDevices::pdf(file = name_fig)
            }
            
            graphics::plot(c(1 : n_Case), m, type = "l", lty = 2, lwd = 2, col = "white", yaxt = "n", xaxt = "n",
            xlim = c(1 - 0.25, n_Case + 0.25), ylim = c(-1, 101),
            cex.axis = cex_par, cex.lab = cex_par, xlab = "", ylab = "Corrected methylation degree", main = title_name)
            llim = max(0, MethylCal_mean - q_Control * MethylCal_sd)
            ulim = min(100, MethylCal_mean + q_Control * MethylCal_sd)
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
            
            MethylCal_idx = (llim > m | m > ulim)
            n_MethylCal_idx = sum(MethylCal_idx, na.rm = TRUE)
            if (n_MethylCal_idx > 0)
            {
                graphics::points(which(MethylCal_idx), rep(100, n_MethylCal_idx), type = "p", pch = 25, lwd = 1.25, cex = 1.25, bg = "red", col = "red")
            }
            
            if (Sys_info[[1]] != "Windows")
            {
                grDevices::dev.off()
            }
            
            MethylCal_idx = which(MethylCal_idx)
            if ((printing == TRUE) & length(MethylCal_idx > 0))
            {
                print("MethylCal hyper/hypomethylated samples:")
                print(Case[Case_ix][MethylCal_idx])
                print(round(m[MethylCal_idx], 2))
                cat(rep("\n", 1))
            }
        }
    }
    
    return(data)
}
