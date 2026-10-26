#' MethylCal correction of case/control samples
#'
#' Correction of case/control samples using MethylCal calibration.
#'
#' @param data Formatted input data frame obtained from the function
#'   \code{\link{Formatting}}.
#' @param Target Name of the target DMR/CpG island/gene to be visualised.
#' @param prior Prior distribution setup for the random effects and the
#'   Latent Gaussian Field \insertCite{Rue2009}{MethylCal}. Three different
#'   priors are implemented: \itemize{ \item \code{LG}: log-Gamma prior is the
#'   default prior with \code{a = 1} and \code{b = 0.1} parameterization;
#'   \item \code{PC}: Penalized Complexity prior
#'   \insertCite{Simpson2017}{MethylCal}; \item \code{HC}: Half-Cauchy prior
#'   \insertCite{Wang2018}{MethylCal}. }
#' @param n_Control Number of control samples.
#' @param n_Case Number of case samples.
#' @param level_Control Level of significance of the differential methylation
#'   test between case and control corrected samples. Default value is
#'   \code{level_Control = 0.9986501}, which corresponds to the z = 3 quantile
#'   for normally distributed data.
#' @param opt_BoxPlot Boxplot option: If \code{opt_BoxPlot = 0} (default),
#'   BoxPlot is centred around the median, whereas if \code{opt_BoxPlot = 1}
#'   it is centred around the mean.
#' @param dir A user-specified directory where plots are saved in
#'   \code{\link[grDevices]{pdf}} format. If not specified, figures are saved
#'   in the current working directory.
#' @param printing If \code{printing = TRUE} (default), the corrected
#'   methylation levels of the case/control samples using MethylCal
#'   calibration \insertCite{Ochoa2019}{MethylCal} are printed to the screen.
#' @param plotting If \code{plotting = TRUE} (default), the corrected
#'   methylation levels for the control samples as well as their
#'   (1-\code{level_Control})\% confidence interval are depicted. BoxPlots
#'   of the corrected methylation levels are also shown for each control
#'   sample. If cases are included, a second figure presents a BoxPlot for
#'   each case sample as well as the controls' (1-\code{level_Control})\%
#'   confidence interval and the hyper- and hypo-methylated cases (top red
#'   triangles).
#' @param cex_par Number indicating the amount by which plotting text and
#'   symbols should be scaled relative to the default (\code{cex_par = 1}).
#' @keywords Case and control samples, MethylCal calibration
#'
#' @export
#'
#' @return This function returns the corrected methylation level for the
#'   control and (if selected) case samples using MethylCal calibration.
#'   Based on a parametric t-test at (1-\code{level_Control})\%, hyper- and
#'   hypo-methylated cases are also flagged.
#'
#'   Figures are saved in the specified directory (or the current working
#'   directory if none is given) in \code{\link[grDevices]{pdf}} format.
#'
#' @references \insertRef{Ochoa2019}{MethylCal}
#' @references \insertRef{Rue2009}{MethylCal}
#' @references \insertRef{Simpson2017}{MethylCal}
#' @references \insertRef{Wang2018}{MethylCal}
#'
#' @examples
#' \dontrun{
#'    # Apply MethylCal Bayesian calibration to case and control samples and test for
#'    # differential methylation. The Celiac example includes patients from two time
#'    # points and centres boxplots on the sample mean.
#'  
#'    data(BWS_data)
#'    AMP = c(0, 25, 50, 75, 100)
#'    data = Formatting(BWS_data, AMP = AMP, n_Control = 15)
#'    corr_data = MethylCalCorrection(data, Target = "KCNQ1OT1", n_Control = 15)
#'  
#'    data(Celiac_data)
#'    AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#'    data = Formatting(Celiac_data, AMP = AMP, n_Control = 13, n_Case = (2 * 17))
#'    corr_data = MethylCalCorrection(data, Target = "NFKBIA", n_Control = 13,
#'                                    n_Case = (2 * 17), opt_BoxPlot = 1)
#'  
#' }


MethylCalCorrection = function(data, Target = NULL, prior = "LG",
                               n_Control = 0, n_Case = 0,
                               level_Control = 0.9986501, opt_BoxPlot = 0,
                               dir = NULL, printing = TRUE, plotting = TRUE,
                               cex_par = 1.25)
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

    # Invert calibration curve for one sample x CpG cell using a pre-fitted
    # model. INLA is fitted once outside; only the spline inversion is done
    # here, using spline_inv_cpp for the inner argmin search.
    f_MethylCal = function(data, CpG_group, fit_MethylCal,
                           Control_Case_name, n_predict = 10^4)
    {
        CpG_group_idx      = which(data$CpG_group == CpG_group)
        CpG_idx            = CpG_group_idx[1]
        Control_Case_idx   = which(colnames(data) == Control_Case_name)
        fit_spline = stats::spline(
            data$x[CpG_group_idx],
            fit_MethylCal$summary.fitted.values$mean[CpG_group_idx],
            method = "natural", n = n_predict
        )
        m = fit_spline$y - data[CpG_idx, Control_Case_idx]
        if (all(diff(m) >= 0))
        {
            x = spline_inv_cpp(fit_spline$y, data[CpG_idx, Control_Case_idx])
        } else {
            m_sq    = m^2
            x_grid  = seq(0, 100, length.out = n_predict)
            idx_top = sort(m_sq, index.return = TRUE)$ix[1:100]
            idx_min = idx_top[which.min(
                (x_grid[idx_top] - data[CpG_idx, Control_Case_idx])^2)]
            x = x_grid[idx_min]
            message("Note: Calibration spline for CpG ", CpG_group,
                    " is non-monotone for sample \"", Control_Case_name,
                    "\". Using nearest-grid approximation")
        }
        return(x)
    }

    plotsegraph = function(loc, mean, min, max, wiskwidth, lwd = 1,
                           col = "black")
    {
        w = wiskwidth / 2
        graphics::segments(x0 = loc, x1 = loc, y0 = min, y1 = max,
                           lwd = lwd, col = col)
        graphics::segments(x0 = loc - w, x1 = loc + w, y0 = min, y1 = min,
                           lwd = lwd, col = col)
        graphics::segments(x0 = loc - w, x1 = loc + w, y0 = max, y1 = max,
                           lwd = lwd, col = col)
    }

    if (is.null(dir)) dir = getwd()
    dir = sub("/+$", "", dir)   # strip trailing slashes to avoid double // in filenames
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    target_safe = gsub("[/\\:*?\"<>|]", "_", Target)
    if (printing) message("Target: ", Target)

    DIC_remove = 1
    n_predict  = 200

    idx      = which(data$Target == Target)
    data     = data[idx, ]
    data$x1  = data$x
    data$x2  = data$x^2
    data$x3  = data$x^3

    data       = data[!is.na(data$CpG_pos), ]
    n_data     = nrow(data)
    n_CpG      = length(unique(data$CpG_pos))
    n_CpG_plot = min(6, n_CpG)
    n_AMP      = length(unique(data$AMP_group_label))
    AMP        = sort(unique(data$x))
    AMP_level  = as.character(data$AMP_group_label)[1:n_AMP]
    data$CpG_group = rep(1:n_CpG, each = n_AMP)
    n_predict  = n_predict + 1
    x_pred     = seq(from = 0, to = 100, length.out = n_predict)
    CpG_group  = rep(NA, n_CpG * n_predict)
    x_pred_plot = rep(NA, n_CpG * n_predict)
    q           = stats::qnorm(0.95 + 0.025)
    q_Control   = stats::qnorm(level_Control + (1 - level_Control) / 2)

    for (s in 1:n_CpG)
    {
        CpG_group[(1:n_predict) + n_predict * (s - 1)]   = s
        x_pred_plot[(1:n_predict) + n_predict * (s - 1)] = x_pred
    }

    data$xCpG_group = data$CpG_group + max(data$CpG_group)

    if (n_Control == 0 & n_Case == 0)
    {
        stop("No controls or cases specified")
    } else {
        Control    = colnames(data)[(10:(9 + n_Control))]
        Control_ix = sort(Control, index.return = TRUE)$ix
        if (n_Case > 0)
        {
            Case    = colnames(data)[(10 + n_Control):(9 + n_Control + n_Case)]
            Case_ix = sort(Case, index.return = TRUE)$ix
        } else {
            Case = NULL
        }
    }

    Control_Case    = c(Control, Case)
    n_Control       = length(Control)
    n_Case          = length(Case)
    n_Control_Case  = length(Control_Case)

    if (!requireNamespace("INLA", quietly = TRUE))
        stop("Package 'INLA' is required but not installed. ",
             "Install it with:\n",
             "  install.packages(\"INLA\",\n",
             "    repos = c(INLA = \"https://inla.r-inla-download.org/R/stable\",\n",
             "              \"https://cloud.r-project.org\"), dep = TRUE)")

    prior_label = prior   # preserve string before list conversion
    if (prior == "LG")
    {
        prior = list(prec = list(prior = "loggamma", param = c(1, 0.1)))
    } else if (prior == "LGEB") {
        stop("Empirical Bayes LogGamma prior to be implemented")
    } else if (prior == "PC") {
        sdy   = stats::sd(data$y, na.rm = TRUE)
        prior = list(prec = list(prior = "pc.prec", param = c(3 * sdy, 0.01)))
    } else if (prior == "HC") {
        halfcauchy = "expression:
                     lambda = 0.022;
                     precision = exp(log_precision);
                     logdens = -1.5 * log_precision - log(pi * lambda)
                              - log(1 + 1 / (precision * lambda^2));
                     log_jacobian = log_precision;
                     return(logdens + log_jacobian);"
        prior = list(prec = list(prior = halfcauchy))
    }

    formula_MethylCal_1 = y ~ x1 + x2 + x3
    formula_MethylCal_2 = y ~ x1 + x2 + x3 +
        f(AMP_group, model = "iid", hyper = prior, constr = TRUE)
    formula_MethylCal_3 = y ~ x1 + x2 + x3 +
        f(AMP_group, model = "iid", hyper = prior, constr = TRUE) +
        f(CpG_group, model = "iid", hyper = prior, constr = TRUE)
    formula_MethylCal_4 = y ~ x1 + x2 + x3 +
        f(AMP_group, model = "iid", hyper = prior, constr = TRUE) +
        f(CpG_pos,   model = "rw1", hyper = prior, scale.model = TRUE)
    formula_MethylCal_5 = y ~ x1 + x2 + x3 +
        f(AMP_group, model = "iid",   hyper = prior, constr = TRUE) +
        f(CpG_group, model = "iid2d", n = 2 * n_CpG, constr = TRUE) +
        f(xCpG_group, AMP_group, copy = "CpG_group")

    inla_ctrl = list(dic = TRUE, waic = TRUE)
    inla_pred = list(compute = TRUE)

    # Safe wrapper: suppresses all INLA output (C stderr via fd-redirect,
    # R messages via muffleMessage, and cat() output via capture.output),
    # imposes an elapsed timeout to prevent indefinite hangs, and returns
    # NULL with one clean diagnostic message on failure.
    # NULL gets DIC = Inf and is therefore never selected.
    safe_inla = function(formula, label, timeout = getOption("MethylCal.inla_timeout", 300))
    {
        make_tip = function()
            switch(prior_label,
                LG   = "try prior = \"PC\" or prior = \"HC\"",
                LGEB = "try prior = \"PC\" or prior = \"HC\"",
                PC   = "try prior = \"HC\"",
                HC   = "check data quality or model specification",
                       "try prior = \"PC\" or prior = \"HC\""   # default
            )

        # C-level stderr suppression (Unix/macOS only; no-op on Windows).
        # Redirects fd 2 to /dev/null, silencing GSL crash traces from the
        # INLA binary before they reach the terminal.
        on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)

        fit = NULL

        if (.Platform$OS.type == "unix") {
            # ── Unix/macOS: guard process (see inla_guard.R) ─────────────
            # The INLA binary can ignore SIGINT and may leak memory after a
            # numerical failure, so a separate guard process watches it and
            # sends SIGKILL (which cannot be blocked) when needed.
            # The guard stops INLA when it runs too long, uses too much memory, or
            # when this R session disappears (see inla_guard.R). Limits can be set
            # with options(MethylCal.inla_timeout = , MethylCal.inla_max_mem_mb = ).
            max_mem_mb = getOption("MethylCal.inla_max_mem_mb", 4096)
            guard      = .inla_guard_start(timeout, max_mem_mb)
            fit_tmp = NULL
            tryCatch(
                withCallingHandlers(
                    { fit_tmp = INLA::inla(formula, family = "gaussian",
                                     verbose           = FALSE,
                                     control.compute   = inla_ctrl,
                                     control.predictor = inla_pred,
                                     data              = data) },
                    message = function(m) invokeRestart("muffleMessage")
                ),
                interrupt = function(i) {
                    # Manual Ctrl+C
                    if (is.null(fit_tmp))
                        message("\nModel [", label, "] interrupted.  ",
                                "Re-run or try a different prior")
                },
                error = function(e) {
                    reason = .inla_guard_reason(guard)
                    if (reason == "memory") {
                        message("\nModel [", label, "] used more than ", max_mem_mb,
                                " MB of memory and has been terminated\n",
                                "  Tip: ", make_tip(), " for greater stability\n",
                                "  with sparse or weakly-identified data")
                    } else if (reason == "time") {
                        message("\nModel [", label, "] exceeded the ", timeout,
                                "s time limit and has been terminated\n",
                                "  Tip: ", make_tip(), " for greater stability\n",
                                "  with sparse or weakly-identified data")
                    } else {
                        message("\nNumerical instability detected for [", label,
                                "]. The model is excluded from selection\n",
                                "  Tip: ", make_tip(), " for greater stability\n",
                                "  with sparse or weakly-identified data")
                    }
                },
                finally = {
                    .inla_guard_stop(guard)
                    # Stop INLA processes that belong to this R session (re-runs)
                    .inla_kill_children()
                }
            )
            fit = fit_tmp
        } else {
            # ── Windows: direct call with setTimeLimit ────────────────────
            # Timeout detection uses elapsed wall time rather than the error
            # message, because inla.core.safe catches the raw setTimeLimit
            # error internally (via its own tryCatch / try(), which sits
            # above our withCallingHandlers in the handler stack) and
            # re-throws a new error whose message no longer contains
            # "time limit".  Genuine numerical failures take milliseconds;
            # a timeout consumes almost exactly `timeout` seconds.
            if (exists("suppress_stderr_cpp")) suppress_stderr_cpp()
            on.exit(
                if (exists("restore_stderr_cpp")) restore_stderr_cpp(),
                add = TRUE
            )
            t_start = proc.time()[["elapsed"]]
            tryCatch(
                {
                    setTimeLimit(elapsed = timeout, transient = TRUE)
                    run_inla = function()
                        withCallingHandlers(
                            { fit <<- INLA::inla(formula, family = "gaussian",
                                           verbose           = FALSE,
                                           control.compute   = inla_ctrl,
                                           control.predictor = inla_pred,
                                           data              = data) },
                            message = function(m) invokeRestart("muffleMessage")
                        )
                    run_inla()
                },
                error = function(e) {
                    elapsed = proc.time()[["elapsed"]] - t_start
                    is_tl   = elapsed >= 0.9 * timeout ||
                              grepl("time limit|elapsed",
                                    conditionMessage(e), ignore.case = TRUE)
                    if (is_tl) {
                        message("\nModel [", label, "] exceeded the ", timeout,
                                "s time limit and has been terminated\n",
                                "  Tip: ", make_tip(), " for greater stability\n",
                                "  with sparse or weakly-identified data")
                    } else {
                        message("\nNumerical instability detected for [", label,
                                "]. The model is excluded from selection.\n",
                                "  Tip: ", make_tip(), " for greater stability\n",
                                "  with sparse or weakly-identified data")
                    }
                }
            )
        }
        fit
    }

    fit_MethylCal_1 = safe_inla(formula_MethylCal_1, paste(Target, "M1"))
    if (printing) message("Model 1 estimated")
    fit_MethylCal_2 = safe_inla(formula_MethylCal_2, paste(Target, "M2"))
    if (printing) message("Model 2 estimated")
    fit_MethylCal_3 = safe_inla(formula_MethylCal_3, paste(Target, "M3"))
    if (printing) message("Model 3 estimated")
    fit_MethylCal_4 = safe_inla(formula_MethylCal_4, paste(Target, "M4"))
    if (printing) message("Model 4 estimated")
    fit_MethylCal_5 = safe_inla(formula_MethylCal_5, paste(Target, "M5"))
    if (printing) message("Model 5 estimated")

    fits = list(fit_MethylCal_1, fit_MethylCal_2, fit_MethylCal_3,
                fit_MethylCal_4, fit_MethylCal_5)
    mlik_MethylCal = sapply(fits, function(f)
        if (is.null(f)) -Inf else f$mlik[2, 1])
    DIC_MethylCal  = sapply(fits, function(f)
        {
            if (is.null(f)) return(Inf)
            d = f$dic$dic
            if (is.null(d) || !is.finite(d)) Inf else d
        })

    if (all(!is.finite(DIC_MethylCal)))
    {
        warning("All INLA models failed for Target ", Target,
                ". Returning NULL.")
        return(invisible(NULL))
    }

    n_obs = nrow(data)
    fitted_MethylCal = vapply(fits, function(f)
        if (is.null(f)) rep(NA_real_, n_obs)
        else            f$summary.fitted.values$mean,
        numeric(n_obs))
    fitted_MethylCal = clamp_matrix_cpp(fitted_MethylCal, 0, 100)

    RSS_MethylCal = vapply(
        seq_len(5),
        function(k) sum((data$y - fitted_MethylCal[, k])^2, na.rm = TRUE),
        numeric(1)
    )

    if (printing)
    {
        message("Estimated models:")
        message(paste(c("",     "Model 2", "Model 3", "Model 4", "Model 5"), collapse = " "))
        message(paste(c("mlik", round(mlik_MethylCal, 2)), collapse = " "))
        message(paste(c("DIC",  round(DIC_MethylCal,  2)), collapse = " "))
        message(paste(c("RSS",  round(RSS_MethylCal,  2)), collapse = " "))
    }

    DIC_tmp             = DIC_MethylCal
    DIC_tmp[DIC_remove] = NA
    MethylCal_idx       = which(DIC_tmp == min(DIC_tmp, na.rm = TRUE))
    if (length(MethylCal_idx) == 0)
    {
        warning("No valid INLA model for Target ", Target,
                ". Returning NULL.")
        return(invisible(NULL))
    }
    MethylCal_idx = MethylCal_idx[1]

    # ── Model selection summary ────────────────────────────────────────
    {
        m_labels   = paste0("M", seq_len(5))
        failed_idx = which(!is.finite(DIC_MethylCal))
        if (length(failed_idx) > 0) {
            message("Excluded: ",
                    paste(m_labels[failed_idx], collapse = ", "),
                    " \u2014 DIC unavailable.")
        } else {
            message("All models estimated successfully")
        }
        message("Selected model: ", m_labels[MethylCal_idx],
                " (DIC = ", round(DIC_MethylCal[MethylCal_idx], 2), ")\n")
    }

    fit_MethylCal = switch(MethylCal_idx,
        fit_MethylCal_1, fit_MethylCal_2, fit_MethylCal_3,
        fit_MethylCal_4, fit_MethylCal_5
    )

    MethylCal_colnames = paste(Control_Case, "MethylCal_corrected", sep = "_")
    data_Control_Case  = matrix(NA, nrow = n_data, ncol = n_Control_Case)
    colnames(data_Control_Case) = MethylCal_colnames
    data = cbind(data, data_Control_Case)

    if (printing) message("Corrected case/control samples:")

    for (i in 1:n_Control_Case)
    {
        Control_Case_name  = Control_Case[i]
        Control_Case_idx   = which(colnames(data) == Control_Case_name)
        data_colnames_idx  = which(colnames(data) == MethylCal_colnames[i])
        if (printing)
            message(paste("Corrected Target:", Target,
                        "- Case/Control:", Control_Case_name))

        for (s in 1:n_CpG)
        {
            CpG_idx  = which(data$CpG_group == s)
            data_CpG = data[CpG_idx, ]
            y_correct = NA
            if (!all(is.na(data_CpG$y)))
            {
                if (is.na(data[CpG_idx, Control_Case_idx][1]))
                {
                    data[CpG_idx[1], data_colnames_idx] = NA
                } else {
                    y_correct = f_MethylCal(
                        data             = data,
                        CpG_group        = s,
                        fit_MethylCal    = fit_MethylCal,
                        Control_Case_name = Control_Case_name
                    )
                    y_correct = max(0, min(100, y_correct))
                    data[CpG_idx[1], data_colnames_idx] = y_correct
                }
            }
            if (printing)
                message(paste("Original:",
                            round(data[CpG_idx[1], Control_Case_idx], 2),
                            "- MethylCal corrected:", round(y_correct, 2)))
        }
    }


    if (plotting)
    {
        title_name = paste(levels(droplevels(data$Target[1])),
                           "MethylCal: Controls", sep = " - ")
        original_data_idx  = match(Control_Case_name, colnames(data))
        MethylCal_data_idx = match(MethylCal_colnames, colnames(data))
        MethylCal_data     = data[seq(1, n_data, n_AMP), MethylCal_data_idx]
        bxp = graphics::boxplot(MethylCal_data[, Control_ix],
                                data = data, plot = FALSE)
        m   = bxp$stats[3, ]
        mn  = bxp$stats[1, ]
        mx  = bxp$stats[5, ]
        if (opt_BoxPlot == 1)
        {
            m = colMeans(MethylCal_data[, Control_ix], na.rm = TRUE)
            names(m) = NULL
        }
        MethylCal_mean = mean(m, na.rm = TRUE)
        MethylCal_sd   = stats::sd(m,   na.rm = TRUE)

                    name_fig = paste0(dir, "/MethylCalCorrection_",
                              target_safe, "_Fig1.pdf")
            grDevices::pdf(file = name_fig)
        llim = max(0,   MethylCal_mean - q_Control * MethylCal_sd)
        ulim = min(100, MethylCal_mean + q_Control * MethylCal_sd)
        graphics::plot(
            seq_len(n_Control), m, type = "l", lty = 2, lwd = 2,
            col = "white", yaxt = "n", xaxt = "n",
            xlim = c(0.75, n_Control + 0.25), ylim = c(-1, 101),
            cex.axis = cex_par, cex.lab = cex_par, xlab = "",
            ylab = "Corrected methylation degree", main = title_name
        )
        graphics::lines(c(0.75, n_Control + 0.25), c(llim, llim),
                        lty = 6, lwd = 2, col = "grey")
        graphics::lines(c(0.75, n_Control + 0.25), c(ulim, ulim),
                        lty = 6, lwd = 2, col = "grey")
        graphics::axis(1, at = seq_len(n_Control),
                       labels = Control[Control_ix], las = 2,
                       cex.axis = cex_par)
        if (n_AMP <= 5)
        {
            graphics::axis(2, at = unique(data$x), cex.axis = cex_par)
        } else {
            graphics::axis(2, at = unique(data$x), labels = FALSE,
                           cex.axis = cex_par)
            x_pos = graphics::par()$usr[1] -
                    0.05 * (graphics::par()$usr[2] - graphics::par()$usr[1])
            graphics::text(x = x_pos, y = unique(data$x),
                           labels = unique(data$x), srt = 0, adj = 1,
                           xpd = TRUE, cex = cex_par)
        }
        plotsegraph(seq_len(n_Control), m, mn, mx, 0.25, col = "black")
        graphics::points(seq_len(n_Control), m, type = "p", pch = 21,
                         lwd = 1.5, cex = 1.05, bg = "white", col = "black")
        grDevices::dev.off()

        if (n_Case > 0)
        {
            title_name = paste(levels(droplevels(data$Target[1])),
                               "MethylCal: Patients", sep = " - ")
            bxp = graphics::boxplot(MethylCal_data[, Case_ix + n_Control],
                                    data = data, plot = FALSE)
            m   = bxp$stats[3, ]
            mn  = bxp$stats[1, ]
            mx  = bxp$stats[5, ]
            if (opt_BoxPlot == 1)
            {
                m = colMeans(MethylCal_data[, Case_ix + n_Control],
                             na.rm = TRUE)
                names(m) = NULL
            }
                            name_fig = paste0(dir, "/MethylCalCorrection_",
                                  target_safe,
                                  "_Fig2.pdf")
                grDevices::pdf(file = name_fig)
            llim = max(0,   MethylCal_mean - q_Control * MethylCal_sd)
            ulim = min(100, MethylCal_mean + q_Control * MethylCal_sd)
            graphics::plot(
                seq_len(n_Case), m, type = "l", lty = 2, lwd = 2,
                col = "white", yaxt = "n", xaxt = "n",
                xlim = c(0.75, n_Case + 0.25), ylim = c(-1, 101),
                cex.axis = cex_par, cex.lab = cex_par, xlab = "",
                ylab = "Corrected methylation degree", main = title_name
            )
            graphics::lines(c(0.75, n_Case + 0.25), c(llim, llim),
                            lty = 6, lwd = 2, col = "grey")
            graphics::lines(c(0.75, n_Case + 0.25), c(ulim, ulim),
                            lty = 6, lwd = 2, col = "grey")
            graphics::axis(1, at = seq_len(n_Case),
                           labels = Case[Case_ix], las = 2,
                           cex.axis = cex_par * (18 / n_Case))
            if (n_AMP <= 5)
            {
                graphics::axis(2, at = unique(data$x), cex.axis = cex_par)
            } else {
                graphics::axis(2, at = unique(data$x), labels = FALSE,
                               cex.axis = cex_par)
                x_pos = graphics::par()$usr[1] -
                        0.05 * (graphics::par()$usr[2] -
                                graphics::par()$usr[1])
                graphics::text(x = x_pos, y = unique(data$x),
                               labels = unique(data$x), srt = 0, adj = 1,
                               xpd = TRUE, cex = cex_par)
            }
            plotsegraph(seq_len(n_Case), m, mn, mx, 0.25, col = "black")
            graphics::points(seq_len(n_Case), m, type = "p", pch = 21,
                             lwd = 1.5, cex = 1.05, bg = "white",
                             col = "black")
            out_idx   = (llim > m | m > ulim)
            n_out_idx = sum(out_idx, na.rm = TRUE)
            if (n_out_idx > 0)
                graphics::points(which(out_idx), rep(100, n_out_idx),
                                 type = "p", pch = 25, lwd = 1.25, cex = 1.25,
                                 bg = "red", col = "red")
            grDevices::dev.off()

            out_idx = which(out_idx)
            if (printing & length(out_idx) > 0)
            {
                message("MethylCal hyper/hypomethylated samples:")
                message(paste(Case[Case_ix][out_idx], collapse = ", "))
                message(paste(round(m[out_idx], 2), collapse = " "))
            }
        }
    }
    return(data)
}
