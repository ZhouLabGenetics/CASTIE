#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double binLogL(double cur_tao, const arma::sp_mat& fam, const arma::sp_mat& W, const arma::vec& Y, const arma::mat& X) {
    // Create V matrix
    arma::sp_mat V = W + cur_tao * fam;
    
    // Cholesky decomposition of V
    arma::sp_mat L;
    bool success = arma::chol(L, V);
    if (!success) {
        Rcpp::Rcerr << "the V matrix is not invertible." << std::endl;
        return NA_REAL;  // Return an invalid value if not invertible
    }

    // Solve for ViX
    arma::mat ViX = solve(L, X);  // L^-1 * X
    if (ViX.has_nan()) {
        Rcpp::Rcerr << "the ViX matrix is not invertible." << std::endl;
        return NA_REAL;
    }

    // Log determinant of V (using Cholesky)
    double logdet_V = sum(log(arma::diagvec(L)));

    // Solve for XtViX and its determinant
    arma::mat XtViX = X.t() * ViX;
    double logdet_XtVX = log(det(XtViX));  // log of the determinant

    // Solve for inv(XtVX)ViX^T
    arma::mat inv_XtVX_ViX = solve(XtViX, ViX.t());

    // Solve for PY
    arma::vec PY = solve(L, Y) - ViX * (inv_XtVX_ViX * Y);

    // Return the log-likelihood
    return -0.5 * (logdet_V + logdet_XtVX + arma::dot(Y, PY));
}

