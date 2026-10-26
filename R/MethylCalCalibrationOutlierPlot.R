#' MethylCal calibration outlier plots
#'
#' Visualisation of MethylCal calibration of the standard control experiment.
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
#' @param level Level of the posterior predictive region.
#' @param dir A user-specified directory where plots are saved in
#'   \code{\link[grDevices]{pdf}} format. If not specified, figures are saved
#'   in the current working directory.
#' @param printing If \code{printing = TRUE} (default), messages are
#'   printed to the screen regarding the estimated models (mlik =
#'   marginal likelihood, DIC = Deviance Information Criterion, RSS =
#'   Residual Sum of Squares) and the correction of the apparent
#'   methylation levels \insertCite{Ochoa2019}{MethylCal}.
#' @param cex_par Number indicating the amount by which plotting text and
#'   symbols should be scaled relative to the default (\code{cex_par = 1}).
#' @keywords Data visualisation, MethylCal calibration, standard control experiment, outliers
#'
#' @export
#'
#' @return This function detects the presence of outliers in the target
#'   dataset, defined as the values outside the posterior predictive
#'   interval [l, u] with l = Q1 - 1.5 IQR and u = Q3 + 1.5 IQR with IQR =
#'   Q3 - Q1 and Q3 and Q1 being the 75th and 25th percentiles. For a
#'   normal posterior distribution [l, u] = [-2.698 std, 2.698 std] with
#'   std the standard deviation of the posterior predictive distribution.
#'
#'   Three plots are generated with black dots if outliers are detected.
#'   In the first plot, the values of the recorded apparent methylation
#'   levels at each Actual Methylation Percentage (AMP) are depicted with
#'   superimposed MethylCal's calibration curve
#'   \insertCite{Ochoa2019}{MethylCal} for each CpG (red dashed line). The
#'   second plot presents the apparent methylation levels at consecutive
#'   CpGs stratified by AMPs with the predicted values superimposed (red
#'   dashed line) as well as the (1-level)\% posterior predictive region
#'   (dashed-dotted red lines). Finally, the third plot is the
#'   scatterplot of the corrected actual methylation percentage at each
#'   AMP for all CpGs within a DMR/CpG island/gene. If outliers are
#'   detected, the models will be refitted after removing the outliers,
#'   and the plots above will be regenerated.
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
#'    # Fit the MethylCal model with automatic outlier detection. Outliers outside the posterior
#'    # predictive interval are highlighted and the model is refitted after their removal,
#'    # producing up to six diagnostic plots.
#'  
#'    data(BWS_data)
#'    AMP = c(0, 25, 50, 75, 100)
#'    data = Formatting(BWS_data, AMP = AMP)
#'    MethylCalCalibrationOutlierPlot(data, Target = "KCNQ1OT1", prior = "HC")
#'  
#'    data(Celiac_data)
#'    AMP = c(0, 12.5, 25, 37.5, 50, 62.5, 87.5, 100)
#'    data = Formatting(Celiac_data, AMP = AMP)
#'    MethylCalCalibrationOutlierPlot(data, Target = "NFKBIA", level = 0.99, printing = FALSE)
#'    
#' }


MethylCalCalibrationOutlierPlot = function(data, Target = NULL, prior = "LG",
                                           level = 0.95, dir = NULL,
                                           printing = TRUE, cex_par = 1.25)
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
        
    f_MethylCal = function(data, CpG_group, AMP_group, fit_MethylCal,
                           n_predict = 10^4)
    {
        CpG_group_idx = which(data$CpG_group == CpG_group)
        CpG_AMP_idx   = which(data$CpG_group == CpG_group &
                               data$AMP_group == AMP_group)
        fit_spline = stats::spline(
            data$x[CpG_group_idx],
            fit_MethylCal$summary.fitted.values$mean[CpG_group_idx],
            method = "natural", n = n_predict
        )
        m = fit_spline$y - data$y[CpG_AMP_idx]
        if (all(diff(m) >= 0))
        {
            x = spline_inv_cpp(fit_spline$y, data$y[CpG_AMP_idx])
        } else {
            m_sq    = m^2
            x_grid  = seq(0, 100, length.out = n_predict)
            idx_top = sort(m_sq, index.return = TRUE)$ix[1:100]
            idx_min = idx_top[which.min(
                (x_grid[idx_top] - data$x[CpG_AMP_idx])^2)]
            x = x_grid[idx_min]
            message("Note: Calibration spline for CpG ", CpG_group,
                    " is non-monotone at AMP group ", AMP_group,
                    ". Using nearest-grid approximation")
        }
        return(x)
    }
    
    
    if (is.null(dir))
    {
        dir = getwd()
    }
    dir = sub("/+$", "", dir)   # strip trailing slashes to avoid double // in filenames
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    target_safe = gsub("[/\\:*?\"<>|]", "_", Target)
    
    if (printing) message("Target: ", Target)
            
    DIC_remove = 1
    n_predict = 200
    sd_min = 0.75
    
    idx = which(data$Target == Target)
    data = data[idx, ]
    data$x1 = data$x
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
    q_outlier = 2.698
    
    for (s in 1 : n_CpG) 
    {
        CpG_idx = which(data$CpG_group == s)
        CpG_group[(1 : n_predict) + n_predict * (s - 1)] = s
        data_CpG = data[CpG_idx, ]
        x_pred_plot[(1 : n_predict) + n_predict * (s - 1)] = x_pred
    }
    
    data$xCpG_group = data$CpG_group + max(data$CpG_group)
    
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
                                "Re-run or try a different prior.")
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
                                "s time limit and has been terminated.\n",
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
                                "s time limit and has been terminated.\n",
                                "  Tip: ", make_tip(), " for greater stability\n",
                                "  with sparse or weakly-identified data")
                    } else {
                        message("\nNumerical instability detected for [", label,
                                "]. The model is excluded from selection\n",
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
        message(paste(c("mlik", round(mlik_MethylCal[-DIC_remove], 2)), collapse = " "))
        message(paste(c("DIC",  round(DIC_MethylCal[-DIC_remove],  2)), collapse = " "))
        message(paste(c("RSS",  round(RSS_MethylCal[-DIC_remove],  2)), collapse = " "))
    }

    DIC_tmp             = DIC_MethylCal
    DIC_tmp[DIC_remove] = NA
    DIC_tmp[!is.finite(DIC_tmp)] = NA     # exclude failed INLA models
    MethylCal_idx       = which(DIC_tmp == min(DIC_tmp, na.rm = TRUE))
    if (length(MethylCal_idx) == 0)
    {
        # No valid INLA model (M2-M5): fall back to M1 (OLS) if it converged
        if (is.finite(DIC_MethylCal[1]))
        {
            message("No valid INLA model for Target ", Target,
                    ". Falling back to M1 (OLS) for calibration.")
            MethylCal_idx = 1L
        } else {
            warning("No valid INLA model for Target ", Target,
                    ". Returning NULL.")
            return(invisible(NULL))
        }
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
    predict_MethylCal = rep(NA, n_CpG * n_predict)
    
    if (MethylCal_idx == 1)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_1$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
        
        fit_MethylCal_1$summary.fitted.values$sd[fit_MethylCal_1$summary.fitted.values$sd < sd_min] = sd_min
        predict_MethylCal_CI_outlier_low = fit_MethylCal_1$summary.fitted.values$mean
        - q_outlier * fit_MethylCal_1$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_high = fit_MethylCal_1$summary.fitted.values$mean
        + q_outlier * fit_MethylCal_1$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low < 0] = 0
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low > 100] = 100
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high < 0 ] = 0
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high > 100] = 100
    }
        if (MethylCal_idx == 2)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_2$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
        
        fit_MethylCal_2$summary.fitted.values$sd[fit_MethylCal_2$summary.fitted.values$sd < sd_min] = sd_min
        predict_MethylCal_CI_outlier_low = fit_MethylCal_2$summary.fitted.values$mean
        - q_outlier * fit_MethylCal_2$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_high = fit_MethylCal_2$summary.fitted.values$mean
        + q_outlier * fit_MethylCal_2$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low < 0] = 0
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low > 100] = 100
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high < 0 ] = 0
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high > 100] = 100
    }
    if (MethylCal_idx == 3)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_3$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
        
        fit_MethylCal_3$summary.fitted.values$sd[fit_MethylCal_3$summary.fitted.values$sd < sd_min] = sd_min
        predict_MethylCal_CI_outlier_low = fit_MethylCal_3$summary.fitted.values$mean
        - q_outlier * fit_MethylCal_3$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_high = fit_MethylCal_3$summary.fitted.values$mean
        + q_outlier * fit_MethylCal_3$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low < 0] = 0
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low > 100] = 100
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high < 0 ] = 0
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high > 100] = 100
    }
    if (MethylCal_idx == 4)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_4$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
        
        fit_MethylCal_4$summary.fitted.values$sd[fit_MethylCal_4$summary.fitted.values$sd < sd_min] = sd_min
        predict_MethylCal_CI_outlier_low = fit_MethylCal_4$summary.fitted.values$mean
        - q_outlier * fit_MethylCal_4$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_high = fit_MethylCal_4$summary.fitted.values$mean
        + q_outlier * fit_MethylCal_4$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low < 0] = 0
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low > 100] = 100
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high < 0 ] = 0
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high > 100] = 100
    }
    if (MethylCal_idx == 5)
    {
        for (s in unique(data$CpG_group))
        {
            CpG_group_idx = which(data$CpG_group == s)
            fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_5$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
        
        fit_MethylCal_5$summary.fitted.values$sd[fit_MethylCal_5$summary.fitted.values$sd < sd_min] = sd_min
        predict_MethylCal_CI_outlier_low = fit_MethylCal_5$summary.fitted.values$mean
        - q_outlier * fit_MethylCal_5$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_high = fit_MethylCal_5$summary.fitted.values$mean
        + q_outlier * fit_MethylCal_5$summary.fitted.values$sd
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low < 0] = 0
        predict_MethylCal_CI_outlier_low[predict_MethylCal_CI_outlier_low > 100] = 100
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high < 0 ] = 0
        predict_MethylCal_CI_outlier_high[predict_MethylCal_CI_outlier_high > 100] = 100
    }
    
    if (printing == TRUE)
    {
        message("Correcting apparent methylation levels ...")
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
                y_correct = f_MethylCal(data = data, CpG_group = s,
                                        AMP_group = l,
                                        fit_MethylCal = fit_MethylCal)
                y_correct[y_correct < 0] = 0
                y_correct[y_correct > 100] = 100
                data$y_correct_MethylCal[l + (s - 1) * n_AMP] = y_correct
                counter = counter + 1
                
                if (printing == TRUE)
                {
                    message("Correcting apparent methylation levels: ",
                            round(counter / (n_AMP * n_CpG) * 100, 1), "%")
                }
            
            }
         
        }
    
    }
    message("Correcting apparent methylation levels ... DONE")
    
    outlier_idx = !(predict_MethylCal_CI_outlier_low <= data$y & data$y <= predict_MethylCal_CI_outlier_high)
    outlier = cbind(predict_MethylCal_CI_outlier_low, data$y,
                   fitted_MethylCal[, MethylCal_idx],
                   predict_MethylCal_CI_outlier_high, outlier_idx)
    data$y_outlier = NA
    data$y_outlier_correct_MethylCal = NA
    data$y_outlier[outlier_idx] = data$y[outlier_idx]
    data$y_outlier_correct_MethylCal[outlier_idx] = data$y_correct_MethylCal[outlier_idx]
    
    title_name = paste(levels(droplevels(data$Target[1])), "MethylCal", sep = " - ")
    sc_xy = list(cex = cex_par,
                 x = list(at = unique(data$x), limits = c(-5, 105), rot = 45),
                 y = list(at = AMP, limits = c(-5, 105)))
    a = lattice::xyplot(
            c(0, 100) ~ c(0, 100),
            type = "l",
            lty = 6,
            lwd = 2,
            col = "grey",
            strip = lattice::strip.custom(bg = "white"),
            data = data,
            scales = sc_xy,
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
            predict_MethylCal ~ x_pred_plot,
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
    e = lattice::xyplot(
            y_outlier ~ x,
            groups = factor(CpG_group),
            type = "p",
            pch = 21,
            lwd = 1.5,
            cex = 1.05,
            fill = "black",
            col = "black",
            strip = lattice::strip.custom(bg = "white"),
            data = data
    )
    p = a + latticeExtra::as.layer(c) + latticeExtra::as.layer(d) + latticeExtra::as.layer(e)
            name_fig = paste0(dir, "/MethylCalCalibrationOutlierPlot_",
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
            y_outlier ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
            type = "p",
            pch = 21,
            lwd = 1.5,
            cex = 0.85,
            fill = "black",
            col = "black",
            strip = lattice::strip.custom(bg = "white"),
            data = data
    )
    e = lattice::xyplot(
            fitted_MethylCal[, MethylCal_idx] ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
            type = "l",
            lty = 2,
            lwd = 2,
            col = "red",
            strip = lattice::strip.custom(bg = "white"),
            data = data
    )
    f = lattice::xyplot(
            predict_MethylCal_CI_low ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
            type = "l",
            lty = 6,
            lwd = 2,
            col = "red",
            strip = lattice::strip.custom(bg = "white"),
            data = data
    )
    g = lattice::xyplot(
            predict_MethylCal_CI_high ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
            type = "l",
            lty = 6,
            lwd = 2,
            col = "red",
            strip = lattice::strip.custom(bg = "white"),
            data = data
    )
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c) +
        latticeExtra::as.layer(d) + latticeExtra::as.layer(e) +
        latticeExtra::as.layer(f) + latticeExtra::as.layer(g)
            name_fig = paste0(dir, "/MethylCalCalibrationOutlierPlot_",
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
            scales = sc_xy,
            xlab = list(cex = cex_par, label = "% Actual methylation"),
            ylab = list(cex = cex_par,
                        label = "% Apparent methylation after PCR"),
            main = title_name
    )
    b = lattice::xyplot(
            y_correct_MethylCal ~ x,
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
    c = lattice::xyplot(
            y_outlier_correct_MethylCal ~ x,
            groups = factor(CpG_group),
            type = "p",
            pch = 21,
            lwd = 1.5,
            cex = 1.05,
            fill = "black",
            col = "black",
            strip = lattice::strip.custom(bg = "white"),
            data = data
    )
    p = a + latticeExtra::as.layer(b) + latticeExtra::as.layer(c)
            name_fig = paste0(dir, "/MethylCalCalibrationOutlierPlot_",
            target_safe,
            "_Fig3.pdf")
        grDevices::pdf(file = name_fig)
        print(p)
        grDevices::dev.off()

    data$y[outlier_idx] = NA

    n_obs_valid = sum(!is.na(data$y))
    if (n_obs_valid == 0)
    {
        warning("All observations flagged as outliers for Target ",
                Target, ". Skipping post-outlier refitting.")
    } else if (any(is.na(data$y)))
    {
        fit_MethylCal_1 = safe_inla(formula_MethylCal_1,
                                    paste(Target, "M1 (post-outlier)"))
        if (printing) message("Model 1 after removing outliers estimated")
        fit_MethylCal_2 = safe_inla(formula_MethylCal_2,
                                    paste(Target, "M2 (post-outlier)"))
        if (printing) message("Model 2 after removing outliers estimated")
        fit_MethylCal_3 = safe_inla(formula_MethylCal_3,
                                    paste(Target, "M3 (post-outlier)"))
        if (printing) message("Model 3 after removing outliers estimated")
        fit_MethylCal_4 = safe_inla(formula_MethylCal_4,
                                    paste(Target, "M4 (post-outlier)"))
        if (printing) message("Model 4 after removing outliers estimated")
        fit_MethylCal_5 = safe_inla(formula_MethylCal_5,
                                    paste(Target, "M5 (post-outlier)"))
        if (printing) message("Model 5 after removing outliers estimated")

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
            warning("All post-outlier INLA models failed for Target ",
                    Target, ". Skipping refitted plots.")
        } else {
        n_obs = nrow(data)
        fitted_MethylCal = vapply(fits, function(f)
            if (is.null(f)) rep(NA_real_, n_obs)
            else            f$summary.fitted.values$mean,
            numeric(n_obs))
        fitted_MethylCal = clamp_matrix_cpp(fitted_MethylCal, 0, 100)

        RSS_MethylCal = vapply(
            seq_len(5),
            function(k) sum((data$y - fitted_MethylCal[, k])^2,
                            na.rm = TRUE),
            numeric(1)
        )

        if (printing == TRUE)
        {
            message("Estimated models after removing outliers:")
            message(paste(c("",     "Model 2", "Model 3", "Model 4", "Model 5"), collapse = " "))
            message(paste(c("mlik", round(mlik_MethylCal[-DIC_remove], 2)), collapse = " "))
            message(paste(c("DIC",  round(DIC_MethylCal[-DIC_remove],  2)), collapse = " "))
            message(paste(c("RSS",  round(RSS_MethylCal[-DIC_remove],  2)), collapse = " "))
        }

        DIC_MethylCal_tmp = DIC_MethylCal
        DIC_MethylCal_tmp[DIC_remove] = NA
        DIC_MethylCal_tmp[!is.finite(DIC_MethylCal_tmp)] = NA  # exclude failed INLA models
        MethylCal_idx = which(DIC_MethylCal_tmp == min(DIC_MethylCal_tmp,
                                                        na.rm = TRUE))
        # If no valid INLA model (M2-M5), fall back to M1 (OLS)
        if (length(MethylCal_idx) == 0 && is.finite(DIC_MethylCal[1]))
        {
            message("No valid post-outlier INLA model for Target ", Target,
                    ". Falling back to M1 (OLS) for post-outlier correction.")
            MethylCal_idx = 1L
        }
        if (length(MethylCal_idx) == 0)
        {
            warning("No valid post-outlier model for Target ",
                    Target, ". Skipping refitted plots.")
        } else {
        MethylCal_idx = MethylCal_idx[1]
        fit_MethylCal = switch(MethylCal_idx,
            fit_MethylCal_1, fit_MethylCal_2, fit_MethylCal_3,
            fit_MethylCal_4, fit_MethylCal_5
        )
        predict_MethylCal = rep(NA, n_CpG * n_predict)
        
        if (MethylCal_idx == 1)
        {
            for (s in unique(data$CpG_group))
            {
                CpG_group_idx = which(data$CpG_group == s)
                fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_1$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
                fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_2$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
                fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_3$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
                fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_4$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
                fit_spline = stats::spline(
                data$x[CpG_group_idx],
                fit_MethylCal_5$summary.fitted.values$mean[CpG_group_idx],
                method = "natural", n = n_predict
            )
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
        
        # ── Model selection summary ─────────────────────────────────────────
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
        if (printing == TRUE)
        {
            message("Correcting apparent methylation levels after removing outliers...")
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
                    CpG_AMP_idx = which(data$CpG_group == s & data$AMP_group == l)
                    
                    if (!is.na(data$y[CpG_AMP_idx]))
                    {
                        y_correct = f_MethylCal(
                            data = data, CpG_group = s,
                            AMP_group = l,
                            fit_MethylCal = fit_MethylCal
                        )
                        y_correct[y_correct < 0] = 0
                        y_correct[y_correct > 100] = 100
                        data$y_correct_MethylCal[l + (s - 1) * n_AMP] = y_correct
                        counter = counter + 1
                        
                        if (printing == TRUE)
                        {
                            message("Correcting: after outlier removal ",
                                    round(counter / (n_AMP * n_CpG) * 100, 1), "%")
                        }
                    }
                
                }
            
            }
        
        }
        message("Correcting apparent methylation levels after removing outliers DONE")
        
        title_name = paste(levels(droplevels(data$Target[1])), "MethylCal", sep = " - ")
        
        a = lattice::xyplot(
                c(0, 100) ~ c(0, 100),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "grey",
                strip = lattice::strip.custom(bg = "white"),
                data = data,
                scales = sc_xy,
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
                predict_MethylCal ~ x_pred_plot,
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
                    name_fig = paste0(dir, "/MethylCalCalibrationOutlierPlot_",
            target_safe,
            "_Fig4.pdf")
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
                fitted_MethylCal[, MethylCal_idx] ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 2,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        e = lattice::xyplot(
                predict_MethylCal_CI_low ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
                type = "l",
                lty = 6,
                lwd = 2,
                col = "red",
                strip = lattice::strip.custom(bg = "white"),
                data = data
        )
        f = lattice::xyplot(
                predict_MethylCal_CI_high ~ CpG_group | factor(AMP_group_label, levels = AMP_level),
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
                    name_fig = paste0(dir, "/MethylCalCalibrationOutlierPlot_",
            target_safe,
            "_Fig5.pdf")
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
                scales = sc_xy,
                xlab = list(cex = cex_par, label = "% Actual methylation"),
                ylab = list(cex = cex_par,
                            label = "% Apparent methylation after PCR"),
                main = title_name
        )
        b = lattice::xyplot(
                y_correct_MethylCal ~ x,
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
        p = a + latticeExtra::as.layer(b)
                    name_fig = paste0(dir, "/MethylCalCalibrationOutlierPlot_",
            target_safe,
            "_Fig6.pdf")
            grDevices::pdf(file = name_fig)
            print(p)
            grDevices::dev.off()
    } # end else (post-outlier models converged)
    } # end else (valid MethylCal_idx found)

    }

}
