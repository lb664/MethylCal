#' Preliminary data exploration
#'
#' Data visualisation of the standard control experiment.
#'
#' @param data Formatted input data frame obtained from the function
#'   \code{\link{Formatting}}.
#' @param Target Name of the target DMR/CpG island/gene to be visualised.
#' @param opt_BoxPlot Boxplot option: If \code{opt_BoxPlot = 0} (default), the
#'   boxplot is centred on the median, whereas if \code{opt_BoxPlot = 1} it
#'   is centred on the mean.
#' @param dir A user-specified directory where plots are saved in
#'   \code{\link[grDevices]{pdf}} format. If not specified, figures are saved
#'   in the current working directory.
#' @param cex_par Number indicating the amount by which plotting text and
#'   symbols should be scaled relative to the default (\code{cex_par = 1}).
#' @param printing If \code{printing = TRUE} (default), a message with the
#'   target name and basic run information is printed to the screen.
#' @keywords Data exploration, data visualisation, standard control experiment
#'
#' @export
#'
#' @return This function returns five plots. The first scatterplot depicts
#' the values of the recorded apparent methylation levels at each Actual
#' Methylation Percentage (AMP). The second one summarises the recorded
#' apparent methylation levels at each AMP by a boxplot, fitting a cubic
#' polynomial regression \insertCite{Moskalev2011}{MethylCal} on the median
#' (or mean) of the data. The third and fourth plots show the apparent
#' methylation percentage after PCR at consecutive CpGs and at each CpG
#' physical position on the chromosome. Finally, the fifth plot presents
#' the apparent methylation levels at consecutive CpGs stratified by AMPs.
#'
#' Figures are saved in the specified directory (or the current working
#' directory if none is given) in \code{\link[grDevices]{pdf}} format.
#'
#' @references \insertRef{Ochoa2019}{MethylCal}
#' @references \insertRef{Moskalev2011}{MethylCal}
#'
#' @examples
#' \dontrun{
#'    # BWS dataset: five calibration AMPs, boxplot centred on the median (default).
#'    data(BWS_data)
#'    AMP = c(0, 25, 50, 75, 100)
#'    data = Formatting(BWS_data, AMP = AMP)
#'    ExploratoryPlot(data, Target = "KCNQ1OT1")
#'  
#'    # Celiac dataset: eight calibration AMPs, boxplot centred on the mean.
#'    data(Celiac_data)
#'    AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#'    data = Formatting(Celiac_data, AMP = AMP)
#'    ExploratoryPlot(data, Target = "NFKBIA", opt_BoxPlot = 1)
#'    
#' }


ExploratoryPlot = function(data, Target = NULL, opt_BoxPlot = 0,
                           dir = NULL, cex_par = 1.25, printing = TRUE)
{
    plotsegraph = function(loc, mean, min, max, wiskwidth, lwd = 1,
                           col = "black") {
        w = wiskwidth / 2
        graphics::segments(x0 = loc, x1 = loc, y0 = min, y1 = max,
                           lwd = lwd, col = col)
        graphics::segments(x0 = loc - w, x1 = loc + w, y0 = min, y1 = min,
                           lwd = lwd, col = col)
        graphics::segments(x0 = loc - w, x1 = loc + w, y0 = max, y1 = max,
                           lwd = lwd, col = col)
    }

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

    data      = data[!is.na(data$CpG_pos), ]
    n_data    = nrow(data)
    p_data    = ncol(data)
    n_CpG     = length(unique(data$CpG_pos))
    n_CpG_plot = min(6, n_CpG)
    n_AMP     = length(unique(data$AMP_group_label))
    AMP       = sort(unique(data$x))
    AMP_level = as.character(data$AMP_group_label)[1 : n_AMP]
    AMP_group = data$AMP_group
    CpG_group = data$CpG_group
    n_predict = n_predict + 1
    x_pred    = seq(from = 0, to = 100, length.out = n_predict)
    title_name = levels(droplevels(data$Target[1]))

    # ---- Figure 1: scatterplot of apparent vs actual methylation ----
    a = lattice::xyplot(
        c(0, 100) ~ c(0, 100), type = "l", lty = 6, lwd = 2, col = "grey",
        strip = lattice::strip.custom(bg = "white"), data = data,
        scales = list(cex = cex_par,
                      x = list(at = AMP, limits = c(-5, 105), rot = 45),
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
        y ~ x, groups = factor(CpG_group), type = "p", pch = 21, lwd = 1.5,
        cex = 1.05, fill = "white", col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    p = a + latticeExtra::as.layer(c)
            name_fig = paste0(dir, "/ExploratoryPlot_",
                          target_safe, "_Fig1.pdf")
        grDevices::pdf(file = name_fig); print(p); grDevices::dev.off()

    # ---- Figure 2: BoxPlot with cubic polynomial fit ----
    graphics::par(mar = c(5, 4, 2, 2) + 0.1)
    boxplot_out = graphics::boxplot(y ~ x, data = data, plot = FALSE)
    m   = boxplot_out$stats[3, ]
    min = boxplot_out$stats[1, ]
    max = boxplot_out$stats[5, ]

    if (opt_BoxPlot == 1)
    {
        m = stats::aggregate(y ~ x, FUN = mean, data = data)$y
    }

    x1 = AMP; x2 = x1 ^ 2; x3 = x1 ^ 3
    fit_lm    = stats::lm(m ~ x1 + x2 + x3)
    c_lm      = stats::coef(summary(fit_lm))[, 1]
    x_pred_plot = x_pred
    predict_c_lm = c_lm[1] + c_lm[2] * x_pred + c_lm[3] * x_pred ^ 2 +
                   c_lm[4] * x_pred ^ 3

            name_fig = paste0(dir, "/ExploratoryPlot_",
                          target_safe, "_Fig2.pdf")
        grDevices::pdf(file = name_fig)

    graphics::plot(
        x_pred_plot, predict_c_lm, type = "l", lty = 2, lwd = 2, col = "black",
        xaxt = "n", yaxt = "n", xlim = c(-1, 101), ylim = c(-1, 101),
        cex.axis = cex_par, cex.lab = cex_par,
        xlab = "% Actual methylation",
        ylab = "% Apparent methylation after PCR", main = title_name
    )
    graphics::lines(x = c(0, 100), y = c(0, 100), type = "l", lty = 6,
                    lwd = 2, col = "grey")
    if (n_AMP <= 5)
    {
        graphics::axis(1, at = AMP, cex.axis = cex_par)
        graphics::axis(2, at = AMP, cex.axis = cex_par)
    } else {
        graphics::axis(1, at = AMP, labels = FALSE, cex.axis = cex_par)
        graphics::axis(2, at = AMP, labels = FALSE, cex.axis = cex_par)
        x_pos = graphics::par()$usr[1] -
                0.05 * (graphics::par()$usr[2] - graphics::par()$usr[1])
        y_pos = graphics::par()$usr[3] -
                0.05 * (graphics::par()$usr[4] - graphics::par()$usr[3])
        graphics::text(cex = cex_par, x = AMP, y = y_pos, labels = AMP,
                       srt = 45, adj = 1, xpd = TRUE)
        graphics::text(cex = cex_par, x = x_pos, y = AMP, labels = AMP,
                       srt = 0, adj = 1, xpd = TRUE)
    }
    plotsegraph(AMP, m, min, max, 5, col = "black")
    graphics::points(AMP, m, type = "p", pch = 21, lwd = 1.5, cex = 1.05,
                     bg = "white", col = "black")
    grDevices::dev.off()

    # ---- Figure 3: apparent methylation vs consecutive CpGs ----
    a = lattice::xyplot(
        x ~ CpG_group, groups = factor(AMP_group), type = "l", lty = 6,
        lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white"),
        data = data,
        scales = list(cex = cex_par,
                      x = list(tick.number = n_CpG_plot),
                      y = list(at = AMP, limits = c(-5, 105))),
        xlab = list(cex = cex_par, label = "CpG"),
        ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"),
        main = title_name
    )
    b = lattice::xyplot(
        y ~ CpG_group, groups = factor(AMP_group), type = "l", lty = 2,
        lwd = 2, col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    c = lattice::xyplot(
        y ~ CpG_group, groups = factor(AMP_group), type = "p", pch = 21,
        lwd = 1.5, cex = 1.05, fill = "white", col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
            name_fig = paste0(dir, "/ExploratoryPlot_",
                          target_safe, "_Fig3.pdf")
        grDevices::pdf(file = name_fig); print(p); grDevices::dev.off()

    # ---- Figure 4: apparent methylation vs genomic position ----
    a = lattice::xyplot(
        x ~ CpG_pos, groups = factor(AMP_group), type = "l", lty = 6,
        lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white"),
        data = data,
        scales = list(cex = cex_par,
                      x = list(tick.number = 3),
                      y = list(at = AMP, limits = c(-5, 105))),
        xlab = list(cex = cex_par, label = "Genomic Position"),
        ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"),
        main = title_name
    )
    b = lattice::xyplot(
        y ~ CpG_pos, groups = factor(AMP_group), type = "l", lty = 2,
        lwd = 2, col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    c = lattice::xyplot(
        y ~ CpG_pos, groups = factor(AMP_group), type = "p", pch = 21,
        lwd = 1.5, cex = 1.05, fill = "white", col = "black",
        strip = lattice::strip.custom(bg = "white"), data = data
    )
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
            name_fig = paste0(dir, "/ExploratoryPlot_",
                          target_safe, "_Fig4.pdf")
        grDevices::pdf(file = name_fig); print(p); grDevices::dev.off()

    # ---- Figure 5: apparent methylation by CpG, faceted by AMP ----
    a = lattice::xyplot(
        x ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
        type = "l", lty = 6, lwd = 2, col = "grey",
        strip = lattice::strip.custom(bg = "white"), data = data,
        scales = list(cex = cex_par,
                      x = list(tick.number = n_CpG_plot),
                      y = list(at = AMP, limits = c(-5, 105))),
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
        type = "p", pch = 21, lwd = 1.5, cex = 1.05, fill = "white",
        col = "black", strip = lattice::strip.custom(bg = "white"), data = data
    )
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
            name_fig = paste0(dir, "/ExploratoryPlot_",
                          target_safe, "_Fig5.pdf")
        grDevices::pdf(file = name_fig); print(p); grDevices::dev.off()
}
