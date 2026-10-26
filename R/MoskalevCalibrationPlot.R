#' Moskalev calibration plots
#'
#' Visualisation of Moskalev calibration of the standard control experiment.
#'
#' @param data Formatted input data frame obtained from the function
#'   \code{\link{Formatting}}.
#' @param Target Name of the target DMR/CpG island/gene to be visualised.
#' @param ModA Extension of Warnecke's calibration curve
#'   \insertCite{Warnecke1997}{MethylCal} by
#'   \insertCite{Moskalev2011}{MethylCal}. Since the default value is
#'   \code{ModA = FALSE}, results of this extension do not appear. To show
#'   results, select \code{ModA = TRUE}.
#' @param ModB Moskalev's cubic polynomial calibration curve
#'   \insertCite{Moskalev2011}{MethylCal}. Default value is \code{ModB = TRUE}.
#' @param level Level of significance of Moskalev's cubic polynomial
#'   prediction interval.
#' @param dir A user-specified directory where plots are saved in
#'   \code{\link[grDevices]{pdf}} format. If not specified, figures are saved
#'   in the current working directory.
#' @param cex_par Number indicating the amount by which plotting text and
#'   symbols should be scaled relative to the default (\code{cex_par = 1}).
#' @param printing If \code{printing = TRUE} (default), a message with the
#'   Target name is printed to the screen.
#' @keywords Data visualisation, Moskalev calibration, standard control experiment
#'
#' @export
#'
#' @return This function returns four plots. The first scatterplot depicts
#'   the values of the recorded apparent methylation levels at each Actual
#'   Methylation Percentage (AMP) with superimposed Moskalev's cubic polynomial
#'   calibration curve \insertCite{Moskalev2011}{MethylCal} for each CpG (red
#'   dashed line). The second plot presents the apparent methylation levels at
#'   consecutive CpGs stratified by AMPs with the predicted values superimposed
#'   (red dashed line) as well as the (1-level)\% prediction interval
#'   (dashed-dotted red lines). The third plot is the scatterplot of the
#'   corrected methylation levels at each AMP for all CpGs within a
#'   DMR/CpG island/gene. Finally, the fourth plot is the hyperbolic
#'   coefficient \insertCite{Moskalev2011}{MethylCal} of the corrected
#'   methylation levels regressed on the AMPs for each CpG.
#'
#'   Figures are saved in the specified directory (or the current working
#'   directory if none is given) in \code{\link[grDevices]{pdf}} format.
#'
#' @references \insertRef{Ochoa2019}{MethylCal}
#' @references \insertRef{Moskalev2011}{MethylCal}
#'
#' @examples
#' \dontrun{
#'    # Generate the Moskalev cubic polynomial calibration plots for a selected DMR
#'    # target. Model B (default) fits a CpG-specific cubic polynomial regression.
#'  
#'    data(BWS_data)
#'    AMP = c(0, 25, 50, 75, 100)
#'    data = Formatting(BWS_data, AMP = AMP)
#'    MoskalevCalibrationPlot(data, Target = "KCNQ1OT1")
#'  
#'    data(Celiac_data)
#'    AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#'    data = Formatting(Celiac_data, AMP = AMP)
#'    MoskalevCalibrationPlot(data, Target = "NFKBIA", level = 0.99)
#'    
#' }


MoskalevCalibrationPlot = function(data, Target = NULL, ModA = FALSE,
                                   ModB = TRUE, level = 0.95, dir = NULL,
                                   cex_par = 1.25, printing = TRUE)
{
    
    if (!is.null(Target)) Target = gsub("\\\\", "/", Target)   # normalise \\ -> / for Windows compatibility
    if (length(unique(data$Target)) > 1)
    {
        if (is.null(Target))
        {
            stop("Please provide the name of the Target for the analysis")
        } else if (sum(data$Target == Target) == 0) {
            stop("Target name not recognised. Available targets: ",
             paste(sort(unique(data$Target)), collapse = ", "))
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
    
    
    if (is.null(dir))
    {
        dir = getwd()
    }
    dir = sub("/+$", "", dir)   # strip trailing slashes to avoid double // in filenames
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    target_safe = gsub("[/\\:*?\"<>|]", "_", Target)
    
    if (printing) message("Target: ", Target)
    
    n_predict = 200
    optimize_interval = c(-10, 110)
    
    idx = which(data$Target == Target)
    data = data[idx, ]
    data$x1  = data$x
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
    
    CpG_group = rep(NA, n_CpG * n_predict)
    b2_lm = rep(NA, n_CpG)
    std_b2_lm = rep(NA, n_CpG)
    x_pred_plot = rep(NA, n_CpG)
    predict_b2_lm = rep(NA, n_CpG * n_predict)
    predict_b2_lm_CI_low = rep(NA, n_CpG * n_AMP)
    predict_b2_lm_CI_high = rep(NA, n_CpG * n_AMP)
    
    for (s in 1 : n_CpG)
    {
        CpG_idx = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        data_CpG$y_tmp = data_CpG$x * (min(data_CpG$y) - data_CpG$y) - 100 * (min(data_CpG$y) - data_CpG$y)
        data_CpG$x_tmp = data_CpG$x * (max(data_CpG$y) - data_CpG$y)
        
        if (!all(is.na(data_CpG$y)))
        {
            fit_lm = stats::lm(y_tmp ~ -1 + x_tmp, data = data_CpG)
            b2_lm[s] = stats::coef(summary(fit_lm))[1]
            std_b2_lm[s] = stats::coef(summary(fit_lm))[2]
            
            x_pred_plot[(1 : n_predict) + n_predict * (s - 1)] = x_pred
            predict_b2_lm[(1 : n_predict) + n_predict * (s - 1)] = ((b2_lm[s] * max(data_CpG$y) - min(data_CpG$y)) * x_pred +
                100 * min(data_CpG$y)) /
                (b2_lm[s] * x_pred - x_pred + 100)
            predict_b2_lm_CI_low[(1 : n_AMP) + n_AMP * (s - 1)] = (((b2_lm[s] - q * std_b2_lm[s]) * max(data_CpG$y) -
                  min(data_CpG$y)) * data_CpG$x + 100 * min(data_CpG$y)) /
                ((b2_lm[s] - q * std_b2_lm[s]) * data_CpG$x -
                  data_CpG$x + 100)
            predict_b2_lm_CI_high[(1 : n_AMP) + n_AMP * (s - 1)] = (((b2_lm[s] + q * std_b2_lm[s]) * max(data_CpG$y) -
                  min(data_CpG$y)) * data_CpG$x + 100 * min(data_CpG$y)) /
                ((b2_lm[s] + q * std_b2_lm[s]) * data_CpG$x -
                  data_CpG$x + 100)
        }
    
    }
    
    predict_b2_lm[predict_b2_lm < 0] = 0
    predict_b2_lm[predict_b2_lm > 100] = 100
    predict_b2_lm_CI_low[predict_b2_lm_CI_low < 0] = 0
    predict_b2_lm_CI_low[predict_b2_lm_CI_low > 100] = 100
    predict_b2_lm_CI_high[predict_b2_lm_CI_high < 0] = 0
    predict_b2_lm_CI_high[predict_b2_lm_CI_high > 100] = 100
    
    data$y_correct_b2_lm = NA
    
    for (s in 1 : n_CpG)
    {
        CpG_idx = which(data$CpG_group == s)
        data_CpG = data[CpG_idx, ]
        y_min = min(data_CpG$y)
        y_max = max(data_CpG$y)
        
        if (!all(is.na(data_CpG$y)))
        {
            y_correct = (100 * y_min - 100 * data_CpG$y) / (b2_lm[s] * data_CpG$y - b2_lm[s] * y_max + y_min - data_CpG$y)
            y_correct[y_correct < 0] = 0
            y_correct[y_correct > 100] = 100
            data$y_correct_b2_lm[CpG_idx] = y_correct
        }
        
    }
    
    correct_b2_lm = rep(NA, n_CpG)
    correct_predict_b2_lm = rep(NA, n_CpG * n_predict)
    
    for (s in 1 : n_CpG)
    {
        CpG_idx = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        data_CpG$y_tmp = (100 - data_CpG$x) * data$y_correct_b2_lm[CpG_idx]
        data_CpG$x_tmp = (100 - data$y_correct_b2_lm[CpG_idx]) * data_CpG$x
        
        if (!all(is.na(data_CpG$y)))
        {
            fit_lm = stats::lm(y_tmp ~ -1 + x_tmp, data = data_CpG)
            correct_b2_lm[s] = stats::coef(summary(fit_lm))[1]
            correct_predict_b2_lm[(1 : n_predict) + n_predict * (s - 1)] = (100 * correct_b2_lm[s] * x_pred) /
                (correct_b2_lm[s] * x_pred - x_pred + 100)
        }
    
    }
    
    if (ModA == TRUE)
    {
        title_name = paste(levels(droplevels(data$Target[1])), "Moskalev etal. 2011a", sep = " - ")
        
        a = lattice::xyplot(
                c(0, 100) ~ c(0, 100),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white"),
                data = data,
                scales = list(cex = cex_par,
                      x = list(at = unique(data$x), limits = c(-5, 105), rot = 45),
                      y = list(at = AMP, limits = c(-5, 105))),
        xlab = list(cex = cex_par, label = "% Actual methylation"),
        ylab = list(cex = cex_par,
                    label = "% Apparent methylation after PCR"),
        main = title_name
    )
        b = lattice::xyplot(
                y ~ x,
                groups = factor(CpG_group),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        c = lattice::xyplot(
                predict_b2_lm ~ x_pred_plot,
                groups = factor(CpG_group),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white")
        )
        d = lattice::xyplot(
                y ~ x,
                groups = factor(CpG_group),
                type = "p",
                pch = 21,
                lwd = 1.5,
                cex = 1.05,
                fill = "white",
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        p = a + latticeExtra::as.layer(c) + latticeExtra::as.layer(d)
                    name_fig = paste0(dir, "/MoskalevCalibrationPlot_ModA_",
                target_safe,
                "_Fig1.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
        
        a = lattice::xyplot(
                x ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white", par.strip.text = list(cex = .70)),
                data = data,
                scales = list(cex = cex_par,
                      x = list(tick.number = n_CpG_plot),
                      y = list(at = unique(data$x), limits = c(-5, 105))),
        layout = c(n_AMP, 1),
        xlab = list(cex = cex_par, label = "CpG"),
        ylab = list(cex = cex_par,
                    label = "% Apparent methylation after PCR"),
        main = title_name
    )
        b = lattice::xyplot(
                y ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        c = lattice::xyplot(
                y ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "p",
                pch = 21,
                lwd = 1.5,
                cex = 0.85,
                fill = "white",
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        d = lattice::xyplot(
                predict_b2_lm_CI_low ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        e = lattice::xyplot(
                predict_b2_lm_CI_high ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c) + latticeExtra::as.layer(d) + latticeExtra::as.layer(e)
                    name_fig = paste0(dir, "/MoskalevCalibrationPlot_ModA_",
                target_safe,
                "_Fig2.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
        
        a = lattice::xyplot(
                c(0, 100) ~ c(0, 100),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white"),
                data = data,
                scales = list(cex = cex_par,
                      x = list(at = unique(data$x), limits = c(-5, 105), rot = 45),
                      y = list(at = AMP, limits = c(-5, 105))),
        xlab = list(cex = cex_par, label = "% Actual methylation"),
        ylab = list(cex = cex_par, label = "% Corrected methylation"),
        main = title_name
    )
        b = lattice::xyplot(
                correct_predict_b2_lm ~ x_pred_plot,
                groups = factor(CpG_group),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white")
        )
        c = lattice::xyplot(
                y_correct_b2_lm ~ x,
                groups = factor(CpG_group),
                type = "p",
                pch = 21,
                lwd = 1.5,
                cex = 1.05,
                fill = "white",
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        p = a + latticeExtra::as.layer(c)
                    name_fig = paste0(dir, "/MoskalevCalibrationPlot_ModA_",
                target_safe,
                "_Fig3.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
        
        a = lattice::xyplot(
                c(1, 1) ~ c(1, n_CpG),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white"),
                data = data,
                scales = list(cex = cex_par,
                              y = list(limits = c(min(0.75, min(correct_b2_lm) * 0.95),
                                                  max(1.25, max(correct_b2_lm) * 1.05)))),
                xlab = list(cex = cex_par, label = "CpG"),
                ylab = list(cex = cex_par, label = "Hyperbolic coefficient corrected"),
                main = title_name
        )
        b = lattice::xyplot(
                correct_b2_lm ~ unique(CpG_group),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        c = lattice::xyplot(
                correct_b2_lm ~ unique(CpG_group),
                type = "p",
                pch = 21,
                lwd = 1.5,
                cex = 1.05,
                fill = "white",
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
                    name_fig = paste0(dir, "/MoskalevCalibrationPlot_ModA_",
                target_safe,
                "_Fig4.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
    }
    
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
            predict_c_lm[(1 : n_predict) + n_predict * (s - 1)] = c_lm[1] + c_lm[2] * x_pred +
                c_lm[3] * x_pred^2 + c_lm[4] * x_pred^3
            predict_c_lm_CI_low[(1 : n_AMP) + n_AMP * (s - 1)] = stats::predict(fit_lm, data_CpG,
                interval = "confidence", level = level)[, 2]
            predict_c_lm_CI_high[(1 : n_AMP) + n_AMP * (s - 1)] = stats::predict(fit_lm, data_CpG,
                interval = "confidence", level = level)[, 3]
        }
    
    }
    
    fitted_c_lm[fitted_c_lm < 0] = 0
    fitted_c_lm[fitted_c_lm > 100] = 100
    predict_c_lm[predict_c_lm < 0] = 0
    predict_c_lm[predict_c_lm > 100] = 100
    predict_c_lm_CI_low[predict_c_lm_CI_low < 0] = 0
    predict_c_lm_CI_low[predict_c_lm_CI_low > 100] = 100
    predict_c_lm_CI_high[predict_c_lm_CI_high < 0] = 0
    predict_c_lm_CI_high[predict_c_lm_CI_high > 100] = 100
    
    data$y_correct_c_lm = NA
    
    for (s in 1 : n_CpG)
    {
        CpG_idx = which(data$CpG_group == s)
        data_CpG = data[CpG_idx, ]
        
        if (!all(is.na(data_CpG$y)))
        {
            fit_lm = stats::lm(y ~ x1 + x2 + x3, data = data_CpG)
            par = stats::coef(summary(fit_lm))[, 1]
            
            for (l in 1 : n_AMP) 
            {
                y_correct = stats::optimize(f_Moskalev, y = data_CpG$y[l], par = par, interval = optimize_interval)$minimum
                y_correct[y_correct < 0] = 0
                y_correct[y_correct > 100] = 100
                data$y_correct_c_lm[l + (s - 1) * n_AMP] = y_correct
            }
        
        }
    
    }
    
    correct_c_lm = rep(NA, n_CpG)
    correct_predict_c_lm = rep(NA, n_CpG * n_predict)
    
    for (s in 1 : n_CpG)
    {
        CpG_idx = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        data_CpG$y_tmp = (100 - data_CpG$x) * data$y_correct_c_lm[CpG_idx]
        data_CpG$x_tmp = (100 - data$y_correct_c_lm[CpG_idx]) * data_CpG$x
        
        if (!all(is.na(data_CpG$y))) {
            fit_lm = stats::lm(y_tmp ~ -1 + x_tmp, data = data_CpG)
            correct_c_lm[s] = stats::coef(summary(fit_lm))[1]
            correct_predict_c_lm[(1 : n_predict) + n_predict * (s - 1)] = (100 * correct_c_lm[s] * x_pred) /
                (correct_c_lm[s] * x_pred - x_pred + 100)
        }
    
    }
    
    if (ModB == TRUE)
    {
        title_name = paste(levels(droplevels(data$Target[1])), "Moskalev etal. 2011b", sep = " - ")
        
        a = lattice::xyplot(
                c(0, 100) ~ c(0, 100),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white"),
                data = data,
                scales = list(cex = cex_par,
                      x = list(at = unique(data$x), limits = c(-5, 105), rot = 45),
                      y = list(at = AMP, limits = c(-5, 105))),
        xlab = list(cex = cex_par, label = "% Actual methylation"),
        ylab = list(cex = cex_par,
                    label = "% Apparent methylation after PCR"),
        main = title_name
    )
        b = lattice::xyplot(
                y ~ x,
                groups = factor(CpG_group),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        c = lattice::xyplot(
                predict_c_lm ~ x_pred_plot,
                groups = factor(CpG_group),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white")
        )
        d = lattice::xyplot(
                y ~ x,
                groups = factor(CpG_group),
                type = "p",
                pch = 21,
                lwd = 1.5,
                cex = 1.05,
                fill = "white",
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        p = a + latticeExtra::as.layer(c) + latticeExtra::as.layer(d)
                    name_fig = paste0(dir, "/MoskalevCalibrationPlot_ModB_",
                target_safe,
                "_Fig1.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
        
        a = lattice::xyplot(
                x ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white", par.strip.text = list(cex = .70)),
                data = data,
                scales = list(cex = cex_par,
                      x = list(tick.number = n_CpG_plot),
                      y = list(at = unique(data$x), limits = c(-5, 105))),
        layout = c(n_AMP, 1),
        xlab = list(cex = cex_par, label = "CpG"),
        ylab = list(cex = cex_par,
                    label = "% Apparent methylation after PCR"),
        main = title_name
    )
        b = lattice::xyplot(
                y ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        c = lattice::xyplot(
                y ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "p",
                pch = 21,
                lwd = 1.5,
                cex = 0.85,
                fill = "white",
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        d = lattice::xyplot(
                fitted_c_lm ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        e = lattice::xyplot(
                predict_c_lm_CI_low ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        f = lattice::xyplot(
                predict_c_lm_CI_high ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c) +
            latticeExtra::as.layer(d) + latticeExtra::as.layer(e) +
            latticeExtra::as.layer(f)
                    name_fig = paste0(dir, "/MoskalevCalibrationPlot_ModB_",
                target_safe,
                "_Fig2.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
        
        a = lattice::xyplot(
                c(0, 100) ~ c(0, 100),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white"),
                data = data,
                scales = list(cex = cex_par,
                      x = list(at = unique(data$x), limits = c(-5, 105), rot = 45),
                      y = list(at = AMP, limits = c(-5, 105))),
        xlab = list(cex = cex_par, label = "% Actual methylation"),
        ylab = list(cex = cex_par, label = "% Corrected methylation"),
        main = title_name
    )
        b = lattice::xyplot(
                correct_predict_c_lm ~ x_pred_plot,
                groups = factor(CpG_group),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white")
        )
        c = lattice::xyplot(
                y_correct_c_lm ~ x,
                groups = factor(CpG_group),
                type = "p",
                pch = 21,
                lwd = 1.5,
                cex = 1.05,
                fill = "white",
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        p = a + latticeExtra::as.layer(c)
                    name_fig = paste0(dir, "/MoskalevCalibrationPlot_ModB_",
                target_safe,
                "_Fig3.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
        
        a = lattice::xyplot(
                c(1, 1) ~ c(1, n_CpG),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white"),
                data = data,
                scales = list(cex = cex_par,
                              y = list(limits = c(min(0.75, min(correct_c_lm) * 0.95),
                                                  max(1.25, max(correct_c_lm) * 1.05)))),
                xlab = list(cex = cex_par, label = "CpG"),
                ylab = list(cex = cex_par, label = "Hyperbolic coefficient corrected"),
                main = title_name
        )
        b = lattice::xyplot(
                correct_c_lm ~ unique(CpG_group),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        c = lattice::xyplot(
                correct_c_lm ~ unique(CpG_group),
                type = "p",
                pch = 21,
                lwd = 1.5,
                cex = 1.05,
                fill = "white",
                col = "black",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
                    name_fig = paste0(dir, "/MoskalevCalibrationPlot_ModB_",
                target_safe,
                "_Fig4.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
    
    }
 
}
