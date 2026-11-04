#' Warnecke calibration plots
#'
#' Visualisation of Warnecke calibration of the standard control experiment.
#'
#' @param data Formatted input data frame obtained from the function
#'   \code{\link{Formatting}}.
#' @param Target Name of the target DMR/CpG island/gene to be visualised.
#' @param level Level of significance of the prediction interval obtained with
#'   Warnecke's calibration curve.
#' @param dir A user-specified directory where plots are saved in
#'   \code{\link[grDevices]{pdf}} format. If not specified, figures are saved
#'   in the current working directory.
#' @param cex_par Number indicating the amount by which plotting text and
#'   symbols should be scaled relative to the default (\code{cex_par = 1}).
#' @param printing If \code{printing = TRUE} (default), a message with the
#'   Target name is printed to the screen.
#' @keywords Data visualisation, Warnecke calibration, standard control experiment
#'
#' @export
#'
#' @return This function returns four plots. The first scatterplot depicts the
#'   values of the recorded apparent methylation levels at each Actual
#'   Methylation Percentage (AMP) with superimposed Warnecke's calibration
#'   curve \insertCite{Warnecke1997}{MethylCal} for each CpG (red dashed line).
#'   The second plot presents the apparent methylation levels at consecutive
#'   CpGs stratified by AMPs with the (1-level)\% prediction interval
#'   superimposed (dashed-dotted red lines). The third plot is the scatterplot
#'   of the corrected methylation levels at each AMP for all CpGs within a
#'   DMR/CpG island/gene. Finally, the fourth plot is the hyperbolic
#'   coefficient \insertCite{Moskalev2011}{MethylCal} of the corrected
#'   methylation levels regressed on the AMPs for each CpG.
#'
#'   Figures are saved in the specified directory (or the current working
#'   directory if none is given) in \code{\link[grDevices]{pdf}} format.
#'
#' @references \insertRef{Ochoa2019}{MethylCal}
#' @references \insertRef{Moskalev2011}{MethylCal}
#' @references \insertRef{Warnecke1997}{MethylCal}
#'
#' @examples
#' \dontrun{
#'    # Generate the four Warnecke hyperbolic calibration plots for a selected
#'    # DMR target. The Celiac example requests the 99% prediction interval.
#'  
#'    data(BWS_data)
#'    AMP = c(0, 25, 50, 75, 100)
#'    data = Formatting(BWS_data, AMP = AMP)
#'    WarneckeCalibrationPlot(data, Target = "KCNQ1OT1")
#'  
#'    data(Celiac_data)
#'    AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#'    data = Formatting(Celiac_data, AMP = AMP)
#'    WarneckeCalibrationPlot(data, Target = "NFKBIA", level = 0.99)
#'    
#' }


WarneckeCalibrationPlot = function(data, Target = NULL, level = 0.95,
                                   dir = NULL, cex_par = 1.25, printing = TRUE)
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

    if (is.null(dir)) dir = getwd()
    dir = sub("/+$", "", dir)   # strip trailing slashes to avoid double // in filenames
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    target_safe = gsub("[/\\:*?\"<>|]", "_", Target)
    if (printing) message("Target: ", Target)

    n_predict = 200
    idx       = which(data$Target == Target)
    data      = data[idx, ]
    data$x1   = data$x
    data$x2   = data$x ^ 2
    data$x3   = data$x ^ 3

    data       = data[!is.na(data$CpG_pos), ]
    n_data     = nrow(data)
    p_data     = ncol(data)
    n_CpG      = length(unique(data$CpG_pos))
    n_CpG_plot = min(6, n_CpG)
    n_AMP      = length(unique(data$AMP_group_label))
    AMP        = sort(unique(data$x))
    AMP_level  = as.character(data$AMP_group_label)[1 : n_AMP]
    n_predict  = n_predict + 1
    x_pred     = seq(from = 0, to = 100, length.out = n_predict)
    q          = stats::qnorm(level + (1 - level) / 2)

    CpG_group          = rep(NA, n_CpG * n_predict)
    b1_lm              = rep(NA, n_CpG)
    std_b1_lm          = rep(NA, n_CpG)
    x_pred_plot        = rep(NA, n_CpG * n_predict)
    predict_b1_lm      = rep(NA, n_CpG * n_predict)
    predict_b1_lm_CI_low  = rep(NA, n_CpG * n_AMP)
    predict_b1_lm_CI_high = rep(NA, n_CpG * n_AMP)

    for (s in 1 : n_CpG)
    {
        CpG_idx  = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        data_CpG$y_tmp = (100 - data_CpG$x) * data_CpG$y
        data_CpG$x_tmp = (100 - data_CpG$y) * data_CpG$x

        if (!all(is.na(data_CpG$y)))
        {
            fit_lm    = stats::lm(y_tmp ~ -1 + x_tmp, data = data_CpG)
            b1_lm[s]     = stats::coef(summary(fit_lm))[1]
            std_b1_lm[s] = stats::coef(summary(fit_lm))[2]
            x_pred_plot[(1 : n_predict) + n_predict * (s - 1)] = x_pred
            predict_b1_lm[(1 : n_predict) + n_predict * (s - 1)] =
                (100 * b1_lm[s] * x_pred) /
                (b1_lm[s] * x_pred - x_pred + 100)
            predict_b1_lm_CI_low[(1 : n_AMP) + n_AMP * (s - 1)] =
                (100 * (b1_lm[s] - q * std_b1_lm[s]) * data_CpG$x) /
                ((b1_lm[s] - q * std_b1_lm[s]) * data_CpG$x - data_CpG$x + 100)
            predict_b1_lm_CI_high[(1 : n_AMP) + n_AMP * (s - 1)] =
                (100 * (b1_lm[s] + q * std_b1_lm[s]) * data_CpG$x) /
                ((b1_lm[s] + q * std_b1_lm[s]) * data_CpG$x - data_CpG$x + 100)
        }
    }

    predict_b1_lm           = clamp_vector_cpp(predict_b1_lm,           0, 100)
    predict_b1_lm_CI_low    = clamp_vector_cpp(predict_b1_lm_CI_low,    0, 100)
    predict_b1_lm_CI_high   = clamp_vector_cpp(predict_b1_lm_CI_high,   0, 100)

    y_min = 0; y_max = 100
    data$y_correct_b1_lm = NA

    for (s in 1 : n_CpG)
    {
        CpG_idx  = which(data$CpG_group == s)
        data_CpG = data[CpG_idx, ]
        if (!all(is.na(data_CpG$y)))
        {
            y_correct =
                (100 * y_min - 100 * data_CpG$y) /
                (b1_lm[s] * data_CpG$y - b1_lm[s] * y_max +
                 y_min - data_CpG$y)
            y_correct = clamp_vector_cpp(y_correct, 0, 100)
            data$y_correct_b1_lm[CpG_idx] = y_correct
        }
    }

    correct_b1_lm         = rep(NA, n_CpG)
    correct_predict_b1_lm = rep(NA, n_CpG * n_predict)

    for (s in 1 : n_CpG)
    {
        CpG_idx  = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        data_CpG$y_tmp = (100 - data_CpG$x) * data$y_correct_b1_lm[CpG_idx]
        data_CpG$x_tmp = (100 - data$y_correct_b1_lm[CpG_idx]) * data_CpG$x
        if (!all(is.na(data_CpG$y)))
        {
            fit_lm = stats::lm(y_tmp ~ -1 + x_tmp, data = data_CpG)
            correct_b1_lm[s] = stats::coef(summary(fit_lm))[1]
            correct_predict_b1_lm[(1 : n_predict) + n_predict * (s - 1)] =
                (100 * correct_b1_lm[s] * x_pred) /
                (correct_b1_lm[s] * x_pred - x_pred + 100)
        }
    }

    title_name = paste(levels(droplevels(data$Target[1])),
                       "Warnecke etal. 1997", sep = " - ")

    # ---- Figure 1: calibration curve ----
    a = lattice::xyplot(
        c(0, 100) ~ c(0, 100), type = "l", lty = 6, lwd = 2, col = "grey",
        strip = lattice::strip.custom(bg = "white"), data = data,
        scales = list(cex = cex_par,
                      x = list(at = unique(data$x), limits = c(-5, 105),
                                rot = 45),
                      y = list(at = AMP, limits = c(-5, 105))),
        xlab = list(cex = cex_par, label = "% Actual methylation"),
        ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"),
        main = title_name
    )
    b = lattice::xyplot(
        y ~ x, groups = factor(CpG_group), type = "l", lty = 2, lwd = 2,
        col = "black", strip = lattice::strip.custom(bg = "white"), data = data
    )
    c = lattice::xyplot(
        predict_b1_lm ~ x_pred_plot, groups = factor(CpG_group), type = "l",
        lty = 2, lwd = 2, col = "red",
        strip = lattice::strip.custom(bg = "white")
    )
    d = lattice::xyplot(
        y ~ x, groups = factor(CpG_group), type = "p", pch = 21, lwd = 1.5,
        cex = 1.05, fill = "white", col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    p = a + latticeExtra::as.layer(c) + latticeExtra::as.layer(d)
            name_fig = paste0(dir, "/WarneckeCalibrationPlot_",
                          target_safe, "_Fig1.pdf")
        grDevices::pdf(file = name_fig); print(p); grDevices::dev.off()

    # ---- Figure 2: CpG-stratified apparent methylation with CI ----
    a = lattice::xyplot(
        x ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
        type = "l", lty = 6, lwd = 2, col = "grey",
        strip = lattice::strip.custom(bg = "white",
                                      par.strip.text = list(cex = .70)),
        data = data,
        scales = list(cex = cex_par,
                      x = list(tick.number = n_CpG_plot),
                      y = list(at = unique(data$x), limits = c(-5, 105))),
        layout = c(n_AMP, 1),
        xlab = list(cex = cex_par, label = "CpG"),
        ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"),
        main = title_name
    )
    b = lattice::xyplot(
        y ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
        type = "l", lty = 2, lwd = 2, col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    c = lattice::xyplot(
        y ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
        type = "p", pch = 21, lwd = 1.5, cex = 0.85, fill = "white",
        col = "black", strip = lattice::strip.custom(bg = "white"), data = data
    )
    d = lattice::xyplot(
        predict_b1_lm_CI_low ~ CpG_group |
            factor(AMP_group_label, levels = AMP_level),
        type = "l", lty = 6, lwd = 2, col = "red",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    e = lattice::xyplot(
        predict_b1_lm_CI_high ~ CpG_group |
            factor(AMP_group_label, levels = AMP_level),
        type = "l", lty = 6, lwd = 2, col = "red",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c) +
        latticeExtra::as.layer(d) + latticeExtra::as.layer(e)
            name_fig = paste0(dir, "/WarneckeCalibrationPlot_",
                          target_safe, "_Fig2.pdf")
        grDevices::pdf(file = name_fig); print(p); grDevices::dev.off()

    # ---- Figure 3: corrected methylation scatter (bug fix: add layer b) ----
    a = lattice::xyplot(
        c(0, 100) ~ c(0, 100), type = "l", lty = 6, lwd = 2, col = "grey",
        strip = lattice::strip.custom(bg = "white"), data = data,
        scales = list(cex = cex_par,
                      x = list(at = unique(data$x), limits = c(-5, 105),
                                rot = 45),
                      y = list(at = AMP, limits = c(-5, 105))),
        xlab = list(cex = cex_par, label = "% Actual methylation"),
        ylab = list(cex = cex_par, label = "% Corrected methylation"),
        main = title_name
    )
    b = lattice::xyplot(
        correct_predict_b1_lm ~ x_pred_plot, groups = factor(CpG_group),
        type = "l", lty = 2, lwd = 2, col = "red",
        strip = lattice::strip.custom(bg = "white")
    )
    c = lattice::xyplot(
        y_correct_b1_lm ~ x, groups = factor(CpG_group), type = "p",
        pch = 21, lwd = 1.5, cex = 1.05, fill = "white", col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
            name_fig = paste0(dir, "/WarneckeCalibrationPlot_",
                          target_safe, "_Fig3.pdf")
        grDevices::pdf(file = name_fig); print(p); grDevices::dev.off()

    # ---- Figure 4: hyperbolic coefficient ----
    a = lattice::xyplot(
        c(1, 1) ~ c(1, n_CpG), type = "l", lty = 6, lwd = 2, col = "grey",
        strip = lattice::strip.custom(bg = "white"), data = data,
        scales = list(cex = cex_par,
                      y = list(limits = c(min(0.75, min(correct_b1_lm) * 0.95),
                                          max(1.25, max(correct_b1_lm) * 1.05)))),
        xlab = list(cex = cex_par, label = "CpG"),
        ylab = list(cex = cex_par, label = "Hyperbolic coefficient corrected"),
        main = title_name
    )
    b = lattice::xyplot(
        correct_b1_lm ~ unique(CpG_group), type = "l", lty = 2, lwd = 2,
        col = "black", strip = lattice::strip.custom(bg = "white"), data = data
    )
    c = lattice::xyplot(
        correct_b1_lm ~ unique(CpG_group), type = "p", pch = 21, lwd = 1.5,
        cex = 1.05, fill = "white", col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
            name_fig = paste0(dir, "/WarneckeCalibrationPlot_",
                          target_safe, "_Fig4.pdf")
        grDevices::pdf(file = name_fig); print(p); grDevices::dev.off()
}
