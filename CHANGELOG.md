# News/Changelog

Changes in **MethylCal** package version **0.2.0** (14/10/2026):
<ul>
    <li> First release after GitHub version 0.1.0. It consolidates all the work
         previously tracked under the internal numbers 1.0.0, 1.1.0 and 2.0.0,
         now numbered 0.1.1, 0.1.2 and 0.1.3 (listed below)
    <li> Tested with R-INLA <b>26.08.07</b>. A model that fails to fit no longer stops
         the analysis. This was the <em>"The inla-program exited with an error"</em>
         failure of version 0.1.0
    <li> Faster: INLA is fitted once per target instead of once per CpG x AMP cell
    <li> INLA safety guard (Unix/macOS): a separate guard process stops INLA when it exceeds a
         memory limit (default 4096 MB, real memory: physical footprint on macOS), when a model takes
         longer than the time limit (default 300 s), or when the R session that started it is gone.
         It only touches INLA processes of its own R session (the previous clean-up killed every INLA
         process of the user). Reason: with R-INLA 26.08.07 on macOS arm64 the INLA executable leaked
         about 30 MB/s after a numerical failure (about 130 GB in 70 minutes) and could outlive R. Limits:
         <code>options(MethylCal.inla_max_mem_mb = , MethylCal.inla_timeout = )</code>
    <li> Backward compatible with 0.1.0. New arguments (<code>printing</code> in
         <code>ExploratoryPlot</code>, <code>MoskalevCalibrationPlot</code> and
         <code>WarneckeCalibrationPlot</code>, <code>ModA</code>/<code>ModB</code> in
         <code>MoskalevCorrection</code>) are placed at the end of the argument lists, so
         positional calls written for 0.1.0 still work. The <code>TableFormatting</code>
         argument <code>writing</code> is renamed <code>write_out</code>, and <code>writing</code>
         still works with a deprecation warning. Return values and file names match 0.1.0.
         Corrected values differ by at most 0.01 percentage points
    <li> Fixed: <code>MethylCalCalibrationPlot</code> wrote a single file named
         <code>MethylCal[CalibrationPlot_...</code> in 0.1.0. It now writes the three figures
         with the correct names, and targets containing <code>/</code> or <code>:</code> no
         longer fail when saving plots
</ul>

Development history before 0.2.0 (versions 0.1.1 to 0.1.3 were not released on GitHub)

Changes in **MethylCal** package version **0.1.3** (09/06/2026):
<ul>
    <li> Platform-aware INLA timeout. macOS uses a <code>processx</code> background
         timer with <code>pkill -9</code> to hard-kill the INLA binary after
         <code>timeout</code> seconds. Windows uses <code>setTimeLimit()</code>
         with elapsed-time-based classification of timeouts vs. numerical failures
    <li> All user-facing <code>cat()</code> and <code>print()</code> calls replaced
         with <code>message()</code>. A two-line model-selection summary (excluded
         models and selected model with DIC) is now printed before each correction
         step
</ul>

Changes in **MethylCal** package version **0.1.2** (01/06/2026):
<ul>
    <li> All three INLA model-fitting functions now wrap each <code>inla()</code>
         call in a <code>safe_inla()</code> helper. A GMRFLib crash no longer
         aborts the pipeline. The failing model receives <code>DIC = Inf</code>
         and is excluded from selection, and if all five models fail the function
         returns <code>NULL</code> with a warning
    <li> Context-aware prior tip messages: the instability warning now suggests
         the next-stronger prior (e.g. PC to HC) rather than a fixed list
    <li> PDF output enabled on Windows. Plots are now saved to disk on all
         platforms instead of being displayed on screen
</ul>

Changes in **MethylCal** package version **0.1.1** (01/03/2026):
<ul>
    <li> Rcpp back-end added: <code>src/MethylCal_utils.cpp</code> provides
         compiled routines for calibration curve inversion
         (<code>spline_inv_cpp</code>, <code>poly_invert_cpp</code>) and value
         clamping (<code>clamp_vector_cpp</code>, <code>clamp_matrix_cpp</code>)
    <li> INLA models now fitted once per target outside the correction loop.
         The pre-fitted object is passed into <code>f_MethylCal</code> rather
         than being refitted per CpG x AMP cell
    <li> Critical bug fixes: <code>$summary.fitted$mean</code> corrected to
         <code>$summary.fitted.values$mean</code>, <code>$mlik[2]</code> to
         <code>$mlik[2, 1]</code>, inverted printing guard, and figure filenames
         with literal <code>[</code> and duplicate <code>_Fig1</code> suffix
</ul>
