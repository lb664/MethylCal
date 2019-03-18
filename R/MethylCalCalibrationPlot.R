#' MethylCal calibration plots
#'
#' Visualisation of MethylCal calibration of the standard control experiment.
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
#' @param level Level of the posterior predictive region.
#' @param dir In Unix-specific OS, user-specified directory where the
#' plots in \code{\link[grDevices]{pdf}} format are saved. If the directory
#' is not specified, figures are saved in the current working directory.
#' @param printing If \code{printing = TRUE} (default), messages are 
#' printed on the screen regarding the estimated models (mlik = marginal
#' likelihood, DIC = Deviance Information Criteria, RSS = Residual
#' Sum of Squares) and the correction of the apparent methylation levels
#' (Ochoa et al., 2019).
#' @param cex_par Number indicating the amount by which plotting text
#' and symbols should be scaled relative to the default (\code{cex_par = 1}).
#'
#' @keywords Data visualisation, MethylCal calibration, standard control
#' experiment
#'
#' @import INLA
#' @import lattice
#' @import latticeExtra
#'
#' @export
#'
#' @return This function returns three plots. The first scatterplot 
#' depicts the values of the recorded apparent methylation levels at
#' each Actual Methylation Percentage (AMP) with superimposed MethylCal's
#' calibration curve (Ochoa et al., 2019) for each CpG (red dashed
#' line). The second plot presents the apparent methylation levels at 
#' consecutive CpGs stratified by AMPs with superimposed the predicted
#' values (red dashed line) as well as the (1-level)\% posterior predictive
#' region (dashed-dotted red lines). Finally, the third plot is the 
#' scatterplot of the corrected actual methylation percentage at each
#' AMP for all CpGs within a DRM/CpG island/gene.
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
#' data = Formatting(BWS_data, AMP = AMP)
#' MethylCalCalibrationPlot(data, Target = "KCNQ1OT1", prior = "HC")
#'
#' data(Celiac_data)
#' AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#' data = Formatting(Celiac_data, AMP = AMP)
#' MethylCalCalibrationPlot(data, Target = "NFKBIA", level = 0.99, printing = FALSE)

MethylCalCalibrationPlot = function(data, Target = NULL, prior = "LG", level = 0.95, dir = NULL, printing = TRUE, cex_par = 1.25)
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
    
    f_MethylCal = function(data, CpG_group, AMP_group, MethylCal_idx, prior, n_predict = 10 ^4)
    {
        n_data = nrow(data)
        p_data = ncol(data)
        n_CpG = max(data$CpG_group)
        n_AMP = max(data$AMP_group)
        CpG_group_idx = which(data$CpG_group == CpG_group)
        CpG_AMP_idx = which(data$CpG_group == CpG_group & data$AMP_group == AMP_group)
        predict_MethylCal = rep(NA, n_predict)
        
        if (MethylCal_idx == 1)
        {
            formula_MethylCal = y ~ x1 + x2 + x3
        }
        if (MethylCal_idx == 2)
        {
            formula_MethylCal = y ~ x1 + x2 + x3 + f(AMP_group, model = "iid", hyper = prior, constr = TRUE)
        }
        if (MethylCal_idx == 3)
        {
            formula_MethylCal = y ~ x1 + x2 + x3 + f(AMP_group, model = "iid", hyper = prior, constr = TRUE) + f(CpG_group, model = "iid", hyper = prior, constr = TRUE)
        }
        if (MethylCal_idx == 4)
        {
            formula_MethylCal = y ~ x1 + x2 + x3 + f(AMP_group, model = "iid", hyper = prior, constr = TRUE) + f(CpG_pos, model = "rw1", hyper = prior, scale.model = TRUE)
        }
        if (MethylCal_idx == 5)
        {
            formula_MethylCal = y ~ x1 + x2 + x3 + f(AMP_group, model = "iid", hyper = prior, constr = TRUE) + f(CpG_group, model = "iid2d", n = 2 * n_CpG, constr = TRUE) + f(xCpG_group, AMP_group, copy = "CpG_group")
        }
        
        fit_MethylCal = inla(formula_MethylCal, family = "gaussian", data = data, control.predictor = list(compute = TRUE))
        fit_spline = stats::spline(data$x[CpG_group_idx], fit_MethylCal$summary.fitted$mean[CpG_group_idx], method = "natural", n = n_predict)
        predict_MethylCal = fit_spline$y
        m = fit_spline$y - data$y[CpG_AMP_idx]
        
        if (all(diff(m) >= 0))
        {
            m = m ^2
            idx_min = which.min(m)
            x = seq(0, 100, length.out = n_predict)[idx_min]
            m = (fit_spline$y - data$y[CpG_AMP_idx])[idx_min]
        } else {
            m = m ^2
            idx_min = sort(m, index.return = TRUE)$ix[1 : (10 ^2)]
            idx_min = idx_min[which.min((seq(0, 100, length.out = n_predict)[idx_min] - data$x[CpG_AMP_idx]) ^2)]
            x = seq(0, 100, length.out = n_predict)[idx_min]
            m = (fit_spline$y - data$y[CpG_AMP_idx])[idx_min]
            
            print("Non-decreasing function")
        }
        
       return(x)
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
    
    for (s in 1 : n_CpG) 
    {
        CpG_idx = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        x_pred_plot[(1 : n_predict) + n_predict * (s - 1)] = x_pred
    }
    
    data$xCpG_group = data$CpG_group + max(data$CpG_group)
    
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
    
    if (printing == TRUE)
    {
        print("Estimated models:")
        print(c("", "Model 1", "Model 2", "Model 3", "Model 4"))
        print(c("mlik", round(mlik_MethylCal[-DIC_remove], 2)))
        print(c("DIC", round(DIC_MethylCal[-DIC_remove], 2)))
        print(c("RSS", round(RSS_MethylCal[-DIC_remove], 2)))
        cat(rep("\n", 1))
    }
    
    DIC_MethylCal_tmp = DIC_MethylCal
    DIC_MethylCal_tmp[DIC_remove] = NA
    MethylCal_idx = which(DIC_MethylCal_tmp == min(DIC_MethylCal_tmp, na.rm = TRUE))
    predict_MethylCal = rep(NA, n_CpG * n_predict)
    
    if (MethylCal_idx == 1)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(data$x[CpG_group_idx], fit_MethylCal_1$summary.fitted$mean[CpG_group_idx], method = "natural", n = n_predict)
            predict_MethylCal[(1 : n_predict) + n_predict * (s - 1)] = fit_spline$y
        }
        predict_MethylCal[predict_MethylCal < 0] = 0
        predict_MethylCal[predict_MethylCal > 100] = 100
        predict_MethylCal_CI_low = fit_MethylCal_1$summary.fitted.values$mean - q * fit_MethylCal_1$summary.fitted.values$sd
        predict_MethylCal_CI_high = fit_MethylCal_1$summary.fitted.values$mean + q * fit_MethylCal_1$summary.fitted.values$sd
        predict_MethylCal_CI_low[predict_MethylCal_CI_low < 0] = 0
        predict_MethylCal_CI_low[predict_MethylCal_CI_low > 100] = 100
        predict_MethylCal_CI_high[predict_MethylCal_CI_high < 0 ] = 0
        predict_MethylCal_CI_high[predict_MethylCal_CI_high > 100] = 100
    }
    if (MethylCal_idx == 2)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(data$x[CpG_group_idx], fit_MethylCal_2$summary.fitted$mean[CpG_group_idx], method = "natural", n = n_predict)
            predict_MethylCal[(1 : n_predict) + n_predict * (s - 1)] = fit_spline$y
        }
        predict_MethylCal[predict_MethylCal < 0] = 0
        predict_MethylCal[predict_MethylCal > 100] = 100
        predict_MethylCal_CI_low = fit_MethylCal_2$summary.fitted.values$mean - q * fit_MethylCal_2$summary.fitted.values$sd
        predict_MethylCal_CI_high = fit_MethylCal_2$summary.fitted.values$mean + q * fit_MethylCal_2$summary.fitted.values$sd
        predict_MethylCal_CI_low[predict_MethylCal_CI_low < 0] = 0
        predict_MethylCal_CI_low[predict_MethylCal_CI_low > 100] = 100
        predict_MethylCal_CI_high[predict_MethylCal_CI_high < 0 ] = 0
        predict_MethylCal_CI_high[predict_MethylCal_CI_high > 100] = 100
    }
    if (MethylCal_idx == 3)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(data$x[CpG_group_idx], fit_MethylCal_3$summary.fitted$mean[CpG_group_idx], method = "natural", n = n_predict)
            predict_MethylCal[(1 : n_predict) + n_predict * (s - 1)] = fit_spline$y
        }
        predict_MethylCal[predict_MethylCal < 0] = 0
        predict_MethylCal[predict_MethylCal > 100] = 100
        predict_MethylCal_CI_low = fit_MethylCal_3$summary.fitted.values$mean - q * fit_MethylCal_3$summary.fitted.values$sd
        predict_MethylCal_CI_high = fit_MethylCal_3$summary.fitted.values$mean + q * fit_MethylCal_3$summary.fitted.values$sd
        predict_MethylCal_CI_low[predict_MethylCal_CI_low < 0] = 0
        predict_MethylCal_CI_low[predict_MethylCal_CI_low > 100] = 100
        predict_MethylCal_CI_high[predict_MethylCal_CI_high < 0 ] = 0
        predict_MethylCal_CI_high[predict_MethylCal_CI_high > 100] = 100
    }
    if (MethylCal_idx == 4)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(data$x[CpG_group_idx], fit_MethylCal_4$summary.fitted$mean[CpG_group_idx], method = "natural", n = n_predict)
            predict_MethylCal[(1 : n_predict) + n_predict * (s - 1)] = fit_spline$y
        }
        predict_MethylCal[predict_MethylCal < 0] = 0
        predict_MethylCal[predict_MethylCal > 100] = 100
        predict_MethylCal_CI_low = fit_MethylCal_4$summary.fitted.values$mean - q * fit_MethylCal_4$summary.fitted.values$sd
        predict_MethylCal_CI_high = fit_MethylCal_4$summary.fitted.values$mean + q * fit_MethylCal_4$summary.fitted.values$sd
        predict_MethylCal_CI_low[predict_MethylCal_CI_low < 0] = 0
        predict_MethylCal_CI_low[predict_MethylCal_CI_low > 100] = 100
        predict_MethylCal_CI_high[predict_MethylCal_CI_high < 0 ] = 0
        predict_MethylCal_CI_high[predict_MethylCal_CI_high > 100] = 100
    }
    if (MethylCal_idx == 5)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(data$x[CpG_group_idx], fit_MethylCal_5$summary.fitted$mean[CpG_group_idx], method = "natural", n = n_predict)
            predict_MethylCal[(1 : n_predict) + n_predict * (s - 1)] = fit_spline$y
        }
        predict_MethylCal[predict_MethylCal < 0] = 0
        predict_MethylCal[predict_MethylCal > 100] = 100
        predict_MethylCal_CI_low = fit_MethylCal_5$summary.fitted.values$mean - q * fit_MethylCal_5$summary.fitted.values$sd
        predict_MethylCal_CI_high = fit_MethylCal_5$summary.fitted.values$mean + q * fit_MethylCal_5$summary.fitted.values$sd
        predict_MethylCal_CI_low[predict_MethylCal_CI_low < 0] = 0
        predict_MethylCal_CI_low[predict_MethylCal_CI_low > 100] = 100
        predict_MethylCal_CI_high[predict_MethylCal_CI_high < 0 ] = 0
        predict_MethylCal_CI_high[predict_MethylCal_CI_high > 100] = 100
    }
    
    if (printing != TRUE)
    {
        print("Correcting apparent methylation levels ...")
    }
    
    data$y_correct_MethylCal = NA
    counter = 0
    
    for (s in 1 : n_CpG)
    {
        CpG_idx = which(data$CpG_group == s)
        data_CpG = data[CpG_idx, ]
            
        if (!all(is.na(data_CpG$y)))
        {
        
            for (l in 1 : n_AMP)
            {
                y_correct = f_MethylCal(data = data, CpG_group = s, AMP_group = l, MethylCal_idx = MethylCal_idx, prior = prior)
                y_correct[y_correct < 0] = 0
                y_correct[y_correct > 100] = 100
                data$y_correct_MethylCal[l + (s - 1) * n_AMP] = y_correct
                counter = counter + 1
                
                if (printing == TRUE)
                {
                    print(paste(paste("Correcting apparent methylation levels:", round(counter / (n_AMP * n_CpG) * 100, 1), sep = " "), "%", sep = ""))
                }
            
            }
         
        }
    
    }
    print("Correcting apparent methylation levels DONE")
    cat(rep("\n", 1))
    
    title_name = paste(levels(droplevels(data$Target[1])), "MethylCal", sep = " - ")
    
    a = lattice::xyplot(c(0, 100) ~ c(0, 100), type = "l", lty = 6, lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white"), data = data, 
    scales = list(cex = cex_par, x = list(at = unique(data$x), limits = c(-5, 105), rot = 45), y = list(at = AMP, limits = c(-5, 105))), 
    xlab = list(cex = cex_par, label = "% Actual methylation"), ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"), main = title_name)
    b = lattice::xyplot(y ~ x, groups = factor(CpG_group), type = "l", lty = 2, lwd = 2, col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    c = lattice::xyplot(predict_MethylCal ~ x_pred_plot, groups = factor(CpG_group), type = "l", lty = 2, lwd = 2, col = "red", strip = lattice::strip.custom(bg = "white"))
    d = lattice::xyplot(y ~ x, groups = factor(CpG_group), type = "p", pch = 21, lwd = 1.5, cex = 1.05, fill = "white", col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    p = a + latticeExtra::as.layer(c) + latticeExtra::as.layer(d)
    if (Sys_info[[1]] == "Windows")
    {
        grDevices::dev.new()
        print(p)
    } else {
        name_fig = paste(paste(dir, "/", "MethylCal[CalibrationPlot_", levels(droplevels(data$Target[1])), "_Fig1", sep = ""), "pdf", sep = ".")
        grDevices::pdf(file = name_fig)
        print(p)
        grDevices::dev.off()
    }
    
    a = lattice::xyplot(x ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "l", lty = 6, lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white", par.strip.text = list(cex = .70)), data = data, 
    scales = list(cex = cex_par, x = list(tick.number = n_CpG_plot), y = list(at = unique(data$x), limits = c(-5, 105))), layout = c(n_AMP, 1),
    xlab = list(cex = cex_par, label = "CpG"), ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"), main = title_name)
    b = lattice::xyplot(y ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "l", lty = 2, lwd = 2, col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    c = lattice::xyplot(y ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "p", pch = 21, lwd = 1.5, cex = 0.85, fill = "white", col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    d = lattice::xyplot(fitted_MethylCal[, MethylCal_idx] ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "l", lty = 2, lwd = 2, col = "red", strip = lattice::strip.custom(bg = "white"), data = data)
    e = lattice::xyplot(predict_MethylCal_CI_low ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "l", lty = 6, lwd = 2, col = "red", strip = lattice::strip.custom(bg = "white"), data = data)
    f = lattice::xyplot(predict_MethylCal_CI_high ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "l", lty = 6, lwd = 2, col = "red", strip = lattice::strip.custom(bg = "white"), data = data)
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c) + latticeExtra::as.layer(d) + latticeExtra::as.layer(e) + latticeExtra::as.layer(f)
    if (Sys_info[[1]] == "Windows")
    {
        grDevices::dev.new()
        print(p)
    } else {
        name_fig = paste(paste(dir, "/", "MethylCal[CalibrationPlot_", levels(droplevels(data$Target[1])), "_Fig1", sep = ""), "pdf", sep = ".")
        grDevices::pdf(file = name_fig)
        print(p)
        grDevices::dev.off()
    }
    
    a = lattice::xyplot(c(0, 100) ~ c(0, 100), type = "l", lty = 6, lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white"), data = data, 
    scales = list(cex = cex_par, x = list(at = unique(data$x), limits = c(-5, 105), rot = 45), y = list(at = AMP, limits = c(-5, 105))), 
    xlab = list(cex = cex_par, label = "% Actual methylation"), ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"), main = title_name)
    b = lattice::xyplot(y_correct_MethylCal ~ x, groups = factor(CpG_group), type = "p", pch = 21, lwd = 1.5, cex = 1.05, fill = "white", col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    p = a + latticeExtra::as.layer(b)
    if (Sys_info[[1]] == "Windows")
    {
        grDevices::dev.new()
        print(p)
    } else {
        name_fig = paste(paste(dir, "/", "MethylCal[CalibrationPlot_", levels(droplevels(data$Target[1])), "_Fig1", sep = ""), "pdf", sep = ".")
        grDevices::pdf(file = name_fig)
        print(p)
        grDevices::dev.off()
    }
}
