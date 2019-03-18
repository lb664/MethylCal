#' Preliminary data exploration
#'
#' Data visualisation of the standard control experiment.
#'
#' @param data Formatted input data frame obtained from the function
#' \code{\link{Formatting}}.
#' @param Target Name of the target DMR/CpG island/gene to be visualised.
#' @param opt_BoxPlot Boxplot option: If \code{opt_BoxPlot = 0} (default),
#' BoxPlot is centered around the median, whereas if \code{opt_BoxPlot = 1}
#' it is centered around the mean.
#' @param dir In Unix-specific OS, user-specified directory where the
#' plots in \code{\link[grDevices]{pdf}} format are saved. If the directory
#' is not specified, figures are saved in the current working directory.
#' @param cex_par Number indicating the amount by which plotting text
#' and symbols should be scaled relative to the default (\code{cex_par = 1}).
#'
#' @keywords Data exploration, data visualisation, standard control
#' experiment
#'
#' @import INLA
#' @import lattice
#' @import latticeExtra
#'
#' @export
#'
#' @return This function returns five plots. The first scatterplot depicts
#' the values of the recorded apparent methylation levels at each Actual
#' Methylation Percentage (AMP). The second one summirises the recorded
#' apparent methylation levels at each AMP by a BoxPlot, fitting a cubic
#' polynomial regression (Moskalev et al., 2011) on the median (mean)
#' of the data. The third and the four plots show the apparent methylation
#' percentage after PCR at consecutive CpGs and at each CpG physical
#' position on the chromosome. Finally, the fifth plot presents the 
#' apparent methylation levels at consecutive CpGs stratified by AMPs.
#'
#' In Unix-specific OS, figures are saved in the current directory,
#' unless otherwise specified by the user, in \code{\link[grDevices]{pdf}}
#' format. In Windows OS, figures are printed on the screen.
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
#' data = Formatting(BWS_data, AMP = AMP)
#' ExploratoryPlot(data, Target = "KCNQ1OT1")
#'
#' data(Celiac_data)
#' AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#' data = Formatting(Celiac_data, AMP = AMP)
#' ExploratoryPlot(data, Target = "NFKBIA", opt_BoxPlot = 1)


ExploratoryPlot = function(data, Target = NULL, opt_BoxPlot = 0, dir = NULL, cex_par = 1.25)
{   
    plotsegraph = function(loc, mean, min, max, wiskwidth, lwd = 1, col = "back") {
        w = wiskwidth /2
        graphics::segments(x0 = loc, x1 = loc, y0 = min, y1 = max, lwd = lwd, col = col)
        graphics::segments(x0 = loc - w, x1 = loc + w, y0 = min, y1 = min, lwd = lwd, col = col)
        graphics::segments(x0 = loc - w, x1 = loc + w, y0 = max, y1 = max, lwd = lwd, col = col)
    }
    
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
    
    Sys_info = Sys.info()
    
    if (is.null(dir))
    {
        dir = getwd()
    }
    
    cat(rep("\n", 1))
    print(paste("Target", Target, sep = " "))
    cat(rep("\n", 1))
    
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
    AMP_group = data$AMP_group
    CpG_group = data$CpG_group
    n_predict = n_predict + 1
    x_pred = seq(from = 0, to = 100, length.out = n_predict)
    
    title_name = levels(droplevels(data$Target[1]))
    
    a = lattice::xyplot(c(0, 100) ~ c(0, 100), type = "l", lty = 6, lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white"), data = data, 
    scales = list(cex = cex_par, x = list(at = AMP, limits = c(-5, 105), rot = 45), y = list(at = AMP, limits = c(-5, 105))), 
    xlab = list(cex = cex_par, label = "% Actual methylation"), ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"), main = title_name)
    b = lattice::xyplot(y ~ x, groups = factor(CpG_group), type = "l", lty = 2, lwd = 2, col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    c = lattice::xyplot(y ~ x, groups = factor(CpG_group), type = "p", pch = 21, lwd = 1.5, cex = 1.05, fill = "white", col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    p = a + latticeExtra::as.layer(c)
    if (Sys_info[[1]] == "Windows")
    {
        grDevices::dev.new()
        print(p)
    } else {
        name_fig = paste(paste(dir, "/", "ExploratoryPlot_", levels(droplevels(data$Target[1])), "_Fig1", sep = ""), "pdf", sep = ".")
        grDevices::pdf(file = name_fig)
        print(p)
        grDevices::dev.off()
    }
    
    graphics::par(mar = c(5, 4, 2, 2) + 0.1)
    boxplot_out = graphics::boxplot(y ~ x, data = data, plot = FALSE)
    m = boxplot_out$stats[3, ]
    min = boxplot_out$stats[1, ]
    max = boxplot_out$stats[5, ]
    
    if (opt_BoxPlot == 1)
    {
        m = stats::aggregate(y ~ x, FUN = mean, data = data)
        m = m$y
    }
    
    x1 = AMP
    x2 = x1 ^ (2)
    x3 = x1 ^ (3)
    
    fit_lm = stats::lm(m ~ x1 + x2 + x3)
    c_lm = stats::coef(summary(fit_lm))[, 1]
    std_c_lm = stats::coef(summary(fit_lm))[, 2]
    
    fitted_c_lm = stats::fitted(fit_lm)
    x_pred_plot = x_pred
    predict_c_lm = c_lm[1] + c_lm[2] * x_pred + c_lm[3] * x_pred ^2 + c_lm[4] * x_pred ^3
    
    if (Sys_info[[1]] == "Windows")
    {
        grDevices::dev.new()
    } else {
        name_fig = paste(paste(dir, "/", "ExploratoryPlot_", levels(droplevels(data$Target[1])), "_Fig2", sep = ""), "pdf", sep = ".")
        grDevices::pdf(file = name_fig)
    }
    
    graphics::plot(x_pred_plot, predict_c_lm, type = "l", lty = 2, lwd = 2, col = "black", xaxt = "n", yaxt = "n", xlim = c(-1, 101), ylim = c(-1, 101),
    cex.axis = cex_par, cex.lab = cex_par, xlab = "% Actual methylation", ylab = "% Apparent methylation after PCR", main = title_name)
    graphics::lines(x = c(0, 100), y = c(0, 100), type = "l", lty = 6, lwd = 2, col = "grey")
    if (n_AMP <= 5) 
    {
        graphics::axis(1, at = AMP, cex.axis = cex_par)
        graphics::axis(2, at = AMP, cex.axis = cex_par)
    } else {
        graphics::axis(1, at = AMP, labels = FALSE, cex.axis = cex_par)
        graphics::axis(2, at = AMP, labels = FALSE, cex.axis = cex_par)
        x_pos = graphics::par()$usr[1] - 0.05 * (graphics::par()$usr[2] - graphics::par()$usr[1])
        y_pos = graphics::par()$usr[3] - 0.05 * (graphics::par()$usr[4] - graphics::par()$usr[3])
        graphics::text(cex = cex_par, x = AMP, y = y_pos, labels = AMP, srt = 45, adj = 1, xpd = TRUE)
        graphics::text(cex = cex_par, x = x_pos, y = AMP, labels = AMP, srt = 0, adj = 1, xpd = TRUE)
    }
    plot.errbars = plotsegraph(AMP, m, min, max, 5, col = "black")
    graphics::points(AMP, m, type = "p", pch = 21, lwd = 1.5, cex = 1.05, bg = "white", col = "black")
    
    if (Sys_info[[1]] != "Windows")
    {
        grDevices::dev.off()
    }
    
    a = lattice::xyplot(x ~ CpG_group, groups = factor(AMP_group), type = "l", lty = 6, lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white"), data = data, 
    scales = list(cex = cex_par, x = list(tick.number = n_CpG_plot), y = list(at = AMP, limits = c(-5, 105))),
    xlab = list(abbreviate, T, cex = cex_par, label = "CpG"), ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"), main = title_name)
    b = lattice::xyplot(y ~ CpG_group, groups = factor(AMP_group), type = "l", lty = 2, lwd = 2, col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    c = lattice::xyplot(y ~ CpG_group, groups = factor(AMP_group), type = "p", pch = 21, lwd = 1.5, cex = 1.05, fill = "white", col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
    if (Sys_info[[1]] == "Windows")
    {
        grDevices::dev.new()
        print(p)
    } else {
        name_fig = paste(paste(dir, "/", "ExploratoryPlot_", levels(droplevels(data$Target[1])), "_Fig3", sep = ""), "pdf", sep = ".")
        grDevices::pdf(file = name_fig)
        print(p)
        grDevices::dev.off()
    }
    
    a = lattice::xyplot(x ~ CpG_pos, groups = factor(AMP_group), type = "l", lty = 6, lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white"), data = data, 
    scales = list(cex = cex_par, x = list(tick.number = 3), y = list(at = AMP, limits = c(-5, 105))),
    xlab = list(cex = cex_par, label = "Genomic Position"), ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"), main = title_name)
    b = lattice::xyplot(y ~ CpG_pos, groups = factor(AMP_group), type = "l", lty = 2, lwd = 2, col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    c = lattice::xyplot(y ~ CpG_pos, groups = factor(AMP_group), type = "p", pch = 21, lwd = 1.5, cex = 1.05, fill = "white", col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
    if (Sys_info[[1]] == "Windows")
    {
        grDevices::dev.new()
        print(p)
    } else {
        name_fig = paste(paste(dir, "/", "ExploratoryPlot_", levels(droplevels(data$Target[1])), "_Fig4", sep = ""), "pdf", sep = ".")
        grDevices::pdf(file = name_fig)
        print(p)
        grDevices::dev.off()
    }
    
    a = lattice::xyplot(x ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "l", lty = 6, lwd = 2, col = "grey", strip = lattice::strip.custom(bg = "white"), data = data, 
    scales = list(cex = cex_par, x = list(tick.number = n_CpG_plot), y = list(at = AMP, limits = c(-5, 105))), layout = c(n_AMP, 1),
    xlab = list(cex = cex_par, label = "CpG"), ylab = list(cex = cex_par, label = "% Apparent methylation after PCR"), main = title_name)
    b = lattice::xyplot(y ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "l", lty = 2, lwd = 2, col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    c = lattice::xyplot(y ~ CpG_group | factor(AMP_group_label, levels = AMP_level), type = "p", pch = 21, lwd = 1.5, cex = 1.05, fill = "white", col = "black", strip = lattice::strip.custom(bg = "white"), data = data)
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
    if (Sys_info[[1]] == "Windows")
    {
        grDevices::dev.new()
        print(p)
    } else {
        name_fig = paste(paste(dir, "/", "ExploratoryPlot_", levels(droplevels(data$Target[1])), "_Fig5", sep = ""), "pdf", sep = ".")
        grDevices::pdf(file = name_fig)
        print(p)
        grDevices::dev.off()
    }
}
