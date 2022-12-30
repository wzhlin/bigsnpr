// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bedXPtr
SEXP bedXPtr(std::string path, int n, int p);
RcppExport SEXP _bigsnpr_bedXPtr(SEXP pathSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(bedXPtr(path, n, p));
    return rcpp_result_gen;
END_RCPP
}
// bed_colstats
List bed_colstats(Environment obj_bed, const IntegerVector& ind_row, const IntegerVector& ind_col, int ncores);
RcppExport SEXP _bigsnpr_bed_colstats(SEXP obj_bedSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(bed_colstats(obj_bed, ind_row, ind_col, ncores));
    return rcpp_result_gen;
END_RCPP
}
// bed_col_counts_cpp
IntegerMatrix bed_col_counts_cpp(Environment obj_bed, const IntegerVector& ind_row, const IntegerVector& ind_col, int ncores);
RcppExport SEXP _bigsnpr_bed_col_counts_cpp(SEXP obj_bedSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(bed_col_counts_cpp(obj_bed, ind_row, ind_col, ncores));
    return rcpp_result_gen;
END_RCPP
}
// bed_row_counts_cpp
arma::Mat<int> bed_row_counts_cpp(Environment obj_bed, const IntegerVector& ind_row, const IntegerVector& ind_col, int ncores);
RcppExport SEXP _bigsnpr_bed_row_counts_cpp(SEXP obj_bedSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(bed_row_counts_cpp(obj_bed, ind_row, ind_col, ncores));
    return rcpp_result_gen;
END_RCPP
}
// read_bed_scaled
NumericMatrix read_bed_scaled(Environment obj_bed, const IntegerVector& ind_row, const IntegerVector& ind_col, const NumericVector& center, const NumericVector& scale);
RcppExport SEXP _bigsnpr_read_bed_scaled(SEXP obj_bedSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP centerSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(read_bed_scaled(obj_bed, ind_row, ind_col, center, scale));
    return rcpp_result_gen;
END_RCPP
}
// prod_and_rowSumsSq
List prod_and_rowSumsSq(Environment obj_bed, const IntegerVector& ind_row, const IntegerVector& ind_col, const NumericVector& center, const NumericVector& scale, const NumericMatrix& V);
RcppExport SEXP _bigsnpr_prod_and_rowSumsSq(SEXP obj_bedSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(prod_and_rowSumsSq(obj_bed, ind_row, ind_col, center, scale, V));
    return rcpp_result_gen;
END_RCPP
}
// bed_pMatVec4
NumericVector bed_pMatVec4(Environment obj_bed, const IntegerVector& ind_row, const IntegerVector& ind_col, const NumericVector& center, const NumericVector& scale, const NumericVector& x, int ncores);
RcppExport SEXP _bigsnpr_bed_pMatVec4(SEXP obj_bedSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP xSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(bed_pMatVec4(obj_bed, ind_row, ind_col, center, scale, x, ncores));
    return rcpp_result_gen;
END_RCPP
}
// bed_cpMatVec4
NumericVector bed_cpMatVec4(Environment obj_bed, const IntegerVector& ind_row, const IntegerVector& ind_col, const NumericVector& center, const NumericVector& scale, const NumericVector& x, int ncores);
RcppExport SEXP _bigsnpr_bed_cpMatVec4(SEXP obj_bedSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP xSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(bed_cpMatVec4(obj_bed, ind_row, ind_col, center, scale, x, ncores));
    return rcpp_result_gen;
END_RCPP
}
// bed_clumping_chr
void bed_clumping_chr(Environment obj_bed, Environment BM2, const IntegerVector& ind_row, const IntegerVector& ind_col, const NumericVector& center, const NumericVector& scale, const IntegerVector& ordInd, const IntegerVector& rankInd, const NumericVector& pos, double size, double thr, int ncores);
RcppExport SEXP _bigsnpr_bed_clumping_chr(SEXP obj_bedSEXP, SEXP BM2SEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP ordIndSEXP, SEXP rankIndSEXP, SEXP posSEXP, SEXP sizeSEXP, SEXP thrSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM2(BM2SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ordInd(ordIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rankInd(rankIndSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pos(posSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    bed_clumping_chr(obj_bed, BM2, ind_row, ind_col, center, scale, ordInd, rankInd, pos, size, thr, ncores);
    return R_NilValue;
END_RCPP
}
// clumping_chr_cached
arma::sp_mat clumping_chr_cached(Environment BM, Environment BM2, arma::sp_mat sqcor, const IntegerVector& spInd, const IntegerVector& rowInd, const IntegerVector& colInd, const IntegerVector& ordInd, const IntegerVector& rankInd, const NumericVector& pos, const NumericVector& sumX, const NumericVector& denoX, double size, double thr, int ncores);
RcppExport SEXP _bigsnpr_clumping_chr_cached(SEXP BMSEXP, SEXP BM2SEXP, SEXP sqcorSEXP, SEXP spIndSEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP ordIndSEXP, SEXP rankIndSEXP, SEXP posSEXP, SEXP sumXSEXP, SEXP denoXSEXP, SEXP sizeSEXP, SEXP thrSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM2(BM2SEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type sqcor(sqcorSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type spInd(spIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ordInd(ordIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rankInd(rankIndSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pos(posSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sumX(sumXSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type denoX(denoXSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(clumping_chr_cached(BM, BM2, sqcor, spInd, rowInd, colInd, ordInd, rankInd, pos, sumX, denoX, size, thr, ncores));
    return rcpp_result_gen;
END_RCPP
}
// clumping_chr
void clumping_chr(Environment BM, Environment BM2, const IntegerVector& rowInd, const IntegerVector& colInd, const IntegerVector& ordInd, const IntegerVector& rankInd, const NumericVector& pos, const NumericVector& sumX, const NumericVector& denoX, double size, double thr, int ncores);
RcppExport SEXP _bigsnpr_clumping_chr(SEXP BMSEXP, SEXP BM2SEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP ordIndSEXP, SEXP rankIndSEXP, SEXP posSEXP, SEXP sumXSEXP, SEXP denoXSEXP, SEXP sizeSEXP, SEXP thrSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM2(BM2SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ordInd(ordIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rankInd(rankIndSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pos(posSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sumX(sumXSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type denoX(denoXSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    clumping_chr(BM, BM2, rowInd, colInd, ordInd, rankInd, pos, sumX, denoX, size, thr, ncores);
    return R_NilValue;
END_RCPP
}
// snp_colstats
ListOf<NumericVector> snp_colstats(Environment BM, const IntegerVector& rowInd, const IntegerVector& colInd, int ncores);
RcppExport SEXP _bigsnpr_snp_colstats(SEXP BMSEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(snp_colstats(BM, rowInd, colInd, ncores));
    return rcpp_result_gen;
END_RCPP
}
// replaceSNP
void replaceSNP(Environment BM, Environment BM2, const IntegerVector& rowInd, const IntegerVector& colInd);
RcppExport SEXP _bigsnpr_replaceSNP(SEXP BMSEXP, SEXP BM2SEXP, SEXP rowIndSEXP, SEXP colIndSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM2(BM2SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    replaceSNP(BM, BM2, rowInd, colInd);
    return R_NilValue;
END_RCPP
}
// corMat
List corMat(Environment obj, const IntegerVector& rowInd, const IntegerVector& colInd, double size, const NumericVector& thr, const NumericVector& pos, bool fill_diag, int ncores);
RcppExport SEXP _bigsnpr_corMat(SEXP objSEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP sizeSEXP, SEXP thrSEXP, SEXP posSEXP, SEXP fill_diagSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pos(posSEXP);
    Rcpp::traits::input_parameter< bool >::type fill_diag(fill_diagSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(corMat(obj, rowInd, colInd, size, thr, pos, fill_diag, ncores));
    return rcpp_result_gen;
END_RCPP
}
// impute
void impute(Environment BM, int method, int ncores);
RcppExport SEXP _bigsnpr_impute(SEXP BMSEXP, SEXP methodSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    impute(BM, method, ncores);
    return R_NilValue;
END_RCPP
}
// lassosum2
List lassosum2(Environment corr, const NumericVector& beta_hat, const NumericVector& lambda, const NumericVector& delta_plus_one, const IntegerVector& ind_sub, double dfmax, int maxiter, double tol);
RcppExport SEXP _bigsnpr_lassosum2(SEXP corrSEXP, SEXP beta_hatSEXP, SEXP lambdaSEXP, SEXP delta_plus_oneSEXP, SEXP ind_subSEXP, SEXP dfmaxSEXP, SEXP maxiterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type corr(corrSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type delta_plus_one(delta_plus_oneSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_sub(ind_subSEXP);
    Rcpp::traits::input_parameter< double >::type dfmax(dfmaxSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(lassosum2(corr, beta_hat, lambda, delta_plus_one, ind_sub, dfmax, maxiter, tol));
    return rcpp_result_gen;
END_RCPP
}
// ld_scores_sfbm
NumericVector ld_scores_sfbm(Rcpp::Environment X, const IntegerVector& ind_sub, int ncores);
RcppExport SEXP _bigsnpr_ld_scores_sfbm(SEXP XSEXP, SEXP ind_subSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Environment >::type X(XSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_sub(ind_subSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(ld_scores_sfbm(X, ind_sub, ncores));
    return rcpp_result_gen;
END_RCPP
}
// ld_scores
NumericVector ld_scores(Environment obj, const IntegerVector& rowInd, const IntegerVector& colInd, double size, const NumericVector& pos, int ncores);
RcppExport SEXP _bigsnpr_ld_scores(SEXP objSEXP, SEXP rowIndSEXP, SEXP colIndSEXP, SEXP sizeSEXP, SEXP posSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pos(posSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(ld_scores(obj, rowInd, colInd, size, pos, ncores));
    return rcpp_result_gen;
END_RCPP
}
// MLE_alpha
arma::vec& MLE_alpha(arma::vec& par, const std::vector<int>& ind_causal, const NumericVector& log_var, const NumericVector& curr_beta, const NumericVector& alpha_bounds, bool boot, bool verbose);
RcppExport SEXP _bigsnpr_MLE_alpha(SEXP parSEXP, SEXP ind_causalSEXP, SEXP log_varSEXP, SEXP curr_betaSEXP, SEXP alpha_boundsSEXP, SEXP bootSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type ind_causal(ind_causalSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type log_var(log_varSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type curr_beta(curr_betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha_bounds(alpha_boundsSEXP);
    Rcpp::traits::input_parameter< bool >::type boot(bootSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(MLE_alpha(par, ind_causal, log_var, curr_beta, alpha_bounds, boot, verbose));
    return rcpp_result_gen;
END_RCPP
}
// ldpred2_gibbs_auto
List ldpred2_gibbs_auto(Environment corr, const NumericVector& beta_hat, const NumericVector& n_vec, const NumericVector& log_var, const IntegerVector& ind_sub, double p_init, double h2_init, int burn_in, int num_iter, int report_step, bool no_jump_sign, double shrink_corr, bool use_mle, const NumericVector& alpha_bounds, double mean_ld, bool verbose);
RcppExport SEXP _bigsnpr_ldpred2_gibbs_auto(SEXP corrSEXP, SEXP beta_hatSEXP, SEXP n_vecSEXP, SEXP log_varSEXP, SEXP ind_subSEXP, SEXP p_initSEXP, SEXP h2_initSEXP, SEXP burn_inSEXP, SEXP num_iterSEXP, SEXP report_stepSEXP, SEXP no_jump_signSEXP, SEXP shrink_corrSEXP, SEXP use_mleSEXP, SEXP alpha_boundsSEXP, SEXP mean_ldSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type corr(corrSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type n_vec(n_vecSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type log_var(log_varSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_sub(ind_subSEXP);
    Rcpp::traits::input_parameter< double >::type p_init(p_initSEXP);
    Rcpp::traits::input_parameter< double >::type h2_init(h2_initSEXP);
    Rcpp::traits::input_parameter< int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< int >::type report_step(report_stepSEXP);
    Rcpp::traits::input_parameter< bool >::type no_jump_sign(no_jump_signSEXP);
    Rcpp::traits::input_parameter< double >::type shrink_corr(shrink_corrSEXP);
    Rcpp::traits::input_parameter< bool >::type use_mle(use_mleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha_bounds(alpha_boundsSEXP);
    Rcpp::traits::input_parameter< double >::type mean_ld(mean_ldSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(ldpred2_gibbs_auto(corr, beta_hat, n_vec, log_var, ind_sub, p_init, h2_init, burn_in, num_iter, report_step, no_jump_sign, shrink_corr, use_mle, alpha_bounds, mean_ld, verbose));
    return rcpp_result_gen;
END_RCPP
}
// ldpred2_gibbs_one_sampling
NumericMatrix ldpred2_gibbs_one_sampling(Environment corr, const NumericVector& beta_hat, const NumericVector& n_vec, const IntegerVector& ind_sub, double h2, double p, bool sparse, int burn_in, int num_iter);
RcppExport SEXP _bigsnpr_ldpred2_gibbs_one_sampling(SEXP corrSEXP, SEXP beta_hatSEXP, SEXP n_vecSEXP, SEXP ind_subSEXP, SEXP h2SEXP, SEXP pSEXP, SEXP sparseSEXP, SEXP burn_inSEXP, SEXP num_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type corr(corrSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type n_vec(n_vecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_sub(ind_subSEXP);
    Rcpp::traits::input_parameter< double >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< int >::type num_iter(num_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(ldpred2_gibbs_one_sampling(corr, beta_hat, n_vec, ind_sub, h2, p, sparse, burn_in, num_iter));
    return rcpp_result_gen;
END_RCPP
}
// ldpred2_gibbs_one
NumericVector ldpred2_gibbs_one(Environment corr, const NumericVector& beta_hat, const NumericVector& n_vec, const IntegerVector& ind_sub, double h2, double p, bool sparse, int burn_in, int num_iter);
RcppExport SEXP _bigsnpr_ldpred2_gibbs_one(SEXP corrSEXP, SEXP beta_hatSEXP, SEXP n_vecSEXP, SEXP ind_subSEXP, SEXP h2SEXP, SEXP pSEXP, SEXP sparseSEXP, SEXP burn_inSEXP, SEXP num_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type corr(corrSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type n_vec(n_vecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_sub(ind_subSEXP);
    Rcpp::traits::input_parameter< double >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< int >::type num_iter(num_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(ldpred2_gibbs_one(corr, beta_hat, n_vec, ind_sub, h2, p, sparse, burn_in, num_iter));
    return rcpp_result_gen;
END_RCPP
}
// multLinReg
NumericMatrix multLinReg(SEXP obj, const IntegerVector& ind_row, const IntegerVector& ind_col, const NumericMatrix& U, int ncores);
RcppExport SEXP _bigsnpr_multLinReg(SEXP objSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP USEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(multLinReg(obj, ind_row, ind_col, U, ncores));
    return rcpp_result_gen;
END_RCPP
}
// extract_submat_bgen
arma::mat& extract_submat_bgen(std::string filename, const std::vector<size_t>& offsets, arma::mat& X, const IntegerVector& ind_row, const NumericVector& decode, bool dosage, int N, int ncores);
RcppExport SEXP _bigsnpr_extract_submat_bgen(SEXP filenameSEXP, SEXP offsetsSEXP, SEXP XSEXP, SEXP ind_rowSEXP, SEXP decodeSEXP, SEXP dosageSEXP, SEXP NSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type decode(decodeSEXP);
    Rcpp::traits::input_parameter< bool >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_submat_bgen(filename, offsets, X, ind_row, decode, dosage, N, ncores));
    return rcpp_result_gen;
END_RCPP
}
// prod_bgen2
arma::mat& prod_bgen2(std::string filename, const NumericVector& offsets, arma::mat& XY, const arma::mat& Y, const IntegerVector& ind_row, const NumericVector& decode, bool dosage, int N, int max_size, int ncores);
RcppExport SEXP _bigsnpr_prod_bgen2(SEXP filenameSEXP, SEXP offsetsSEXP, SEXP XYSEXP, SEXP YSEXP, SEXP ind_rowSEXP, SEXP decodeSEXP, SEXP dosageSEXP, SEXP NSEXP, SEXP max_sizeSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type XY(XYSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type decode(decodeSEXP);
    Rcpp::traits::input_parameter< bool >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type max_size(max_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(prod_bgen2(filename, offsets, XY, Y, ind_row, decode, dosage, N, max_size, ncores));
    return rcpp_result_gen;
END_RCPP
}
// read_bgen
List read_bgen(std::string filename, const NumericVector& offsets, const Environment& BM, const IntegerVector& ind_row, const IntegerVector& ind_col, const RawVector& decode, bool dosage, int N, int ncores);
RcppExport SEXP _bigsnpr_read_bgen(SEXP filenameSEXP, SEXP offsetsSEXP, SEXP BMSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP decodeSEXP, SEXP dosageSEXP, SEXP NSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< const Environment& >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< const RawVector& >::type decode(decodeSEXP);
    Rcpp::traits::input_parameter< bool >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(read_bgen(filename, offsets, BM, ind_row, ind_col, decode, dosage, N, ncores));
    return rcpp_result_gen;
END_RCPP
}
// readbina
bool readbina(const char * filename, Environment BM, const RawMatrix& tab);
RcppExport SEXP _bigsnpr_readbina(SEXP filenameSEXP, SEXP BMSEXP, SEXP tabSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const RawMatrix& >::type tab(tabSEXP);
    rcpp_result_gen = Rcpp::wrap(readbina(filename, BM, tab));
    return rcpp_result_gen;
END_RCPP
}
// readbina2
void readbina2(Environment BM, Environment obj_bed, const IntegerVector& ind_row, const IntegerVector& ind_col, int ncores);
RcppExport SEXP _bigsnpr_readbina2(SEXP BMSEXP, SEXP obj_bedSEXP, SEXP ind_rowSEXP, SEXP ind_colSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< Environment >::type obj_bed(obj_bedSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_row(ind_rowSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ind_col(ind_colSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    readbina2(BM, obj_bed, ind_row, ind_col, ncores);
    return R_NilValue;
END_RCPP
}
// sp_colSumsSq_sym
NumericVector sp_colSumsSq_sym(std::vector<size_t> p, const IntegerVector& i, const NumericVector& x);
RcppExport SEXP _bigsnpr_sp_colSumsSq_sym(SEXP pSEXP, SEXP iSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<size_t> >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_colSumsSq_sym(p, i, x));
    return rcpp_result_gen;
END_RCPP
}
// get_L
List get_L(std::vector<size_t> p, const IntegerVector& i, const NumericVector& x, double thr_r2, double max_r2);
RcppExport SEXP _bigsnpr_get_L(SEXP pSEXP, SEXP iSEXP, SEXP xSEXP, SEXP thr_r2SEXP, SEXP max_r2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<size_t> >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type thr_r2(thr_r2SEXP);
    Rcpp::traits::input_parameter< double >::type max_r2(max_r2SEXP);
    rcpp_result_gen = Rcpp::wrap(get_L(p, i, x, thr_r2, max_r2));
    return rcpp_result_gen;
END_RCPP
}
// get_C
List get_C(const arma::sp_mat& L, int min_size, int max_size, int max_K, double max_cost, const NumericVector& pos_scaled);
RcppExport SEXP _bigsnpr_get_C(SEXP LSEXP, SEXP min_sizeSEXP, SEXP max_sizeSEXP, SEXP max_KSEXP, SEXP max_costSEXP, SEXP pos_scaledSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type min_size(min_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type max_size(max_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type max_K(max_KSEXP);
    Rcpp::traits::input_parameter< double >::type max_cost(max_costSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pos_scaled(pos_scaledSEXP);
    rcpp_result_gen = Rcpp::wrap(get_C(L, min_size, max_size, max_K, max_cost, pos_scaled));
    return rcpp_result_gen;
END_RCPP
}
// get_perc
double get_perc(const NumericVector& p, const IntegerVector& i, const IntegerVector& all_last);
RcppExport SEXP _bigsnpr_get_perc(SEXP pSEXP, SEXP iSEXP, SEXP all_lastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type all_last(all_lastSEXP);
    rcpp_result_gen = Rcpp::wrap(get_perc(p, i, all_last));
    return rcpp_result_gen;
END_RCPP
}
// writebina
void writebina(const char * filename, Environment BM, const RawVector& tab, const IntegerVector& rowInd, const IntegerVector& colInd);
RcppExport SEXP _bigsnpr_writebina(SEXP filenameSEXP, SEXP BMSEXP, SEXP tabSEXP, SEXP rowIndSEXP, SEXP colIndSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    Rcpp::traits::input_parameter< const RawVector& >::type tab(tabSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type rowInd(rowIndSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type colInd(colIndSEXP);
    writebina(filename, BM, tab, rowInd, colInd);
    return R_NilValue;
END_RCPP
}
// testWrite
void testWrite(const RawVector& v, const char * filename);
RcppExport SEXP _bigsnpr_testWrite(SEXP vSEXP, SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RawVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP);
    testWrite(v, filename);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bigsnpr_bedXPtr", (DL_FUNC) &_bigsnpr_bedXPtr, 3},
    {"_bigsnpr_bed_colstats", (DL_FUNC) &_bigsnpr_bed_colstats, 4},
    {"_bigsnpr_bed_col_counts_cpp", (DL_FUNC) &_bigsnpr_bed_col_counts_cpp, 4},
    {"_bigsnpr_bed_row_counts_cpp", (DL_FUNC) &_bigsnpr_bed_row_counts_cpp, 4},
    {"_bigsnpr_read_bed_scaled", (DL_FUNC) &_bigsnpr_read_bed_scaled, 5},
    {"_bigsnpr_prod_and_rowSumsSq", (DL_FUNC) &_bigsnpr_prod_and_rowSumsSq, 6},
    {"_bigsnpr_bed_pMatVec4", (DL_FUNC) &_bigsnpr_bed_pMatVec4, 7},
    {"_bigsnpr_bed_cpMatVec4", (DL_FUNC) &_bigsnpr_bed_cpMatVec4, 7},
    {"_bigsnpr_bed_clumping_chr", (DL_FUNC) &_bigsnpr_bed_clumping_chr, 12},
    {"_bigsnpr_clumping_chr_cached", (DL_FUNC) &_bigsnpr_clumping_chr_cached, 14},
    {"_bigsnpr_clumping_chr", (DL_FUNC) &_bigsnpr_clumping_chr, 12},
    {"_bigsnpr_snp_colstats", (DL_FUNC) &_bigsnpr_snp_colstats, 4},
    {"_bigsnpr_replaceSNP", (DL_FUNC) &_bigsnpr_replaceSNP, 4},
    {"_bigsnpr_corMat", (DL_FUNC) &_bigsnpr_corMat, 8},
    {"_bigsnpr_impute", (DL_FUNC) &_bigsnpr_impute, 3},
    {"_bigsnpr_lassosum2", (DL_FUNC) &_bigsnpr_lassosum2, 8},
    {"_bigsnpr_ld_scores_sfbm", (DL_FUNC) &_bigsnpr_ld_scores_sfbm, 3},
    {"_bigsnpr_ld_scores", (DL_FUNC) &_bigsnpr_ld_scores, 6},
    {"_bigsnpr_MLE_alpha", (DL_FUNC) &_bigsnpr_MLE_alpha, 7},
    {"_bigsnpr_ldpred2_gibbs_auto", (DL_FUNC) &_bigsnpr_ldpred2_gibbs_auto, 16},
    {"_bigsnpr_ldpred2_gibbs_one_sampling", (DL_FUNC) &_bigsnpr_ldpred2_gibbs_one_sampling, 9},
    {"_bigsnpr_ldpred2_gibbs_one", (DL_FUNC) &_bigsnpr_ldpred2_gibbs_one, 9},
    {"_bigsnpr_multLinReg", (DL_FUNC) &_bigsnpr_multLinReg, 5},
    {"_bigsnpr_extract_submat_bgen", (DL_FUNC) &_bigsnpr_extract_submat_bgen, 8},
    {"_bigsnpr_prod_bgen2", (DL_FUNC) &_bigsnpr_prod_bgen2, 10},
    {"_bigsnpr_read_bgen", (DL_FUNC) &_bigsnpr_read_bgen, 9},
    {"_bigsnpr_readbina", (DL_FUNC) &_bigsnpr_readbina, 3},
    {"_bigsnpr_readbina2", (DL_FUNC) &_bigsnpr_readbina2, 5},
    {"_bigsnpr_sp_colSumsSq_sym", (DL_FUNC) &_bigsnpr_sp_colSumsSq_sym, 3},
    {"_bigsnpr_get_L", (DL_FUNC) &_bigsnpr_get_L, 5},
    {"_bigsnpr_get_C", (DL_FUNC) &_bigsnpr_get_C, 6},
    {"_bigsnpr_get_perc", (DL_FUNC) &_bigsnpr_get_perc, 3},
    {"_bigsnpr_writebina", (DL_FUNC) &_bigsnpr_writebina, 5},
    {"_bigsnpr_testWrite", (DL_FUNC) &_bigsnpr_testWrite, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_bigsnpr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
