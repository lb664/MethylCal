#include <Rcpp.h>
using namespace Rcpp;

//' Clamp a numeric vector to [lo, hi]
//'
//' Values below \code{lo} are set to \code{lo}; values above \code{hi} are
//' set to \code{hi}; \code{NA} values are preserved.
//'
//' @param x   Numeric vector.
//' @param lo  Lower bound (inclusive).
//' @param hi  Upper bound (inclusive).
//' @return    Clamped numeric vector of the same length as \code{x}.
//' @keywords internal
// [[Rcpp::export]]
NumericVector clamp_vector_cpp(NumericVector x, double lo, double hi)
{
    int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++) {
        double v = x[i];
        if (ISNA(v))      out[i] = NA_REAL;
        else if (v < lo)  out[i] = lo;
        else if (v > hi)  out[i] = hi;
        else              out[i] = v;
    }
    return out;
}

//' Clamp a numeric matrix to [lo, hi] (column-major, in-place copy)
//'
//' @param x   Numeric matrix.
//' @param lo  Lower bound (inclusive).
//' @param hi  Upper bound (inclusive).
//' @return    Clamped numeric matrix with the same dimensions as \code{x}.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix clamp_matrix_cpp(NumericMatrix x, double lo, double hi)
{
    int n = x.nrow();
    int p = x.ncol();
    NumericMatrix out(n, p);
    for (int j = 0; j < p; j++) {
        for (int i = 0; i < n; i++) {
            double v = x(i, j);
            if (ISNA(v))      out(i, j) = NA_REAL;
            else if (v < lo)  out(i, j) = lo;
            else if (v > hi)  out(i, j) = hi;
            else              out(i, j) = v;
        }
    }
    return out;
}

//' Find the x-grid value where a spline best approximates a target y
//'
//' Performs a fast C++ linear scan over the pre-computed spline values
//' \code{y_spline} (generated in R via \code{stats::spline()}) and returns
//' the x-grid coordinate at which the squared residual
//' \eqn{(y_{\text{spline}} - y_{\text{target}})^2} is minimised.
//'
//' @param y_spline  Numeric vector of spline values on a fine equi-spaced
//'   grid of length \eqn{n}.
//' @param y_target  Scalar target value to match.
//' @param x_start   Left endpoint of the x-grid (default \code{0}).
//' @param x_end     Right endpoint of the x-grid (default \code{100}).
//' @return          Scalar x value at which the spline best matches
//'   \code{y_target}.
//' @keywords internal
// [[Rcpp::export]]
double spline_inv_cpp(NumericVector y_spline, double y_target,
                      double x_start = 0.0, double x_end = 100.0)
{
    int    n       = y_spline.size();
    double best_sq = R_PosInf;
    int    best_i  = 0;

    for (int i = 0; i < n; i++) {
        double d  = y_spline[i] - y_target;
        double sq = d * d;
        if (sq < best_sq) {
            best_sq = sq;
            best_i  = i;
        }
    }
    return x_start + (x_end - x_start) * (double)best_i / (n - 1.0);
}

//' Invert a cubic polynomial evaluated on a grid (Moskalev inversion helper)
//'
//' Computes \eqn{(a_0 + a_1 x + a_2 x^2 + a_3 x^3 - y)^2} for each element
//' of \code{x_grid} and returns the x value that minimises this expression.
//' Used as a faster alternative to \code{stats::optimize()} for the
//' Moskalev cubic-polynomial inversion step.
//'
//' @param par     Numeric vector of length 4: \eqn{(a_0, a_1, a_2, a_3)}.
//' @param x_grid  Numeric vector of candidate x values.
//' @param y_obs   Observed apparent methylation value to invert.
//' @return        Scalar x value that minimises the squared residual.
//' @keywords internal
// [[Rcpp::export]]
double poly3_inv_cpp(NumericVector par, NumericVector x_grid, double y_obs)
{
    double a0 = par[0], a1 = par[1], a2 = par[2], a3 = par[3];
    int    n  = x_grid.size();
    double best_sq = R_PosInf;
    int    best_i  = 0;

    for (int i = 0; i < n; i++) {
        double xi  = x_grid[i];
        double val = a0 + a1 * xi + a2 * xi * xi + a3 * xi * xi * xi - y_obs;
        double sq  = val * val;
        if (sq < best_sq) {
            best_sq = sq;
            best_i  = i;
        }
    }
    return x_grid[best_i];
}

//' Invert a cubic polynomial on an interval via dense grid search
//'
//' Finds x in [\code{lower}, \code{upper}] minimising
//' \eqn{(a_0 + a_1 x + a_2 x^2 + a_3 x^3 - y_{\text{obs}})^2}.
//' Used by \code{\link{MoskalevCorrection}} to invert the fitted cubic
//' calibration curve for each sample.
//'
//' @param y_obs  Observed apparent methylation value to invert.
//' @param par    Numeric vector of length 4: \eqn{(a_0, a_1, a_2, a_3)}.
//' @param lower  Lower bound of the search interval (default \code{-10}).
//' @param upper  Upper bound of the search interval (default \code{110}).
//' @return       Scalar x minimising the squared residual.
//' @keywords internal
// [[Rcpp::export]]
double poly_invert_cpp(double y_obs, NumericVector par,
                       double lower = -10.0, double upper = 110.0)
{
    double a0 = par[0], a1 = par[1], a2 = par[2], a3 = par[3];
    int    n  = 10000;
    double best_sq = R_PosInf;
    double best_x  = lower;

    for (int i = 0; i < n; i++) {
        double x   = lower + (upper - lower) * (double)i / (n - 1.0);
        double val = a0 + a1 * x + a2 * x * x + a3 * x * x * x - y_obs;
        double sq  = val * val;
        if (sq < best_sq) {
            best_sq = sq;
            best_x  = x;
        }
    }
    return best_x;
}

// ── stderr suppression (silences INLA native-binary GSL output) ──────────────
// Uses POSIX dup/dup2 to redirect fd 2 to /dev/null before launching INLA
// and restore it afterwards.  On Windows the functions are no-ops since
// Windows INLA uses a different subprocess model.

#ifndef _WIN32
#include <unistd.h>
#include <fcntl.h>
static int saved_stderr_fd_g = -1;
#endif

//' Redirect C-level stderr to /dev/null
//'
//' Called immediately before each \code{inla()} invocation to suppress
//' GSL diagnostic messages and segfault traces from INLA's native binary.
//' Must be paired with a call to \code{\link{restore_stderr_cpp}}.
//' On Windows this function is a no-op.
//' @keywords internal
// [[Rcpp::export]]
void suppress_stderr_cpp()
{
#ifndef _WIN32
    saved_stderr_fd_g = dup(2);
    int null_fd = open("/dev/null", O_WRONLY);
    if (null_fd != -1) { dup2(null_fd, 2); close(null_fd); }
#endif
}

//' Restore C-level stderr after suppression
//'
//' Must be called after every \code{\link{suppress_stderr_cpp}} call,
//' including on early exit (use \code{on.exit}).
//' On Windows this function is a no-op.
//' @keywords internal
// [[Rcpp::export]]
void restore_stderr_cpp()
{
#ifndef _WIN32
    if (saved_stderr_fd_g != -1) {
        dup2(saved_stderr_fd_g, 2);
        close(saved_stderr_fd_g);
        saved_stderr_fd_g = -1;
    }
#endif
}
