/******************************************************************************/

#include <bigsparser/SFBM.h>
#include <bigstatsr/utils.h>

/******************************************************************************/

// [[Rcpp::export]]
NumericVector ldpred2_gibbs_one(Environment corr,
                                const NumericVector& beta_hat,
                                const NumericVector& n_vec,
                                const IntegerVector& ind_sub,
                                double h2,
                                double p,
                                bool sparse,
                                int burn_in,
                                int num_iter) {

  XPtr<SFBM> sfbm = corr["address"];

  int m = beta_hat.size();
  myassert_size(n_vec.size(), m);
  NumericVector curr_beta(m), avg_beta(m);  // only for the subset
  int m2 = sfbm->ncol();
  NumericVector dotprods(m2);  // for the full corr

  double h2_per_var = h2 / (m * p);
  double inv_odd_p = (1 - p) / p;
  double gap0 = 2 *
    std::inner_product(beta_hat.begin(), beta_hat.end(), beta_hat.begin(), 0.0);

  for (int k = -burn_in; k < num_iter; k++) {

    double gap = 0;

    for (int j = 0; j < m; j++) {

      int j2 = ind_sub[j];
      double res_beta_hat_j = beta_hat[j] - (dotprods[j2] - curr_beta[j]);

      double C1 = h2_per_var * n_vec[j];
      double C2 = 1 / (1 + 1 / C1);
      double C3 = C2 * res_beta_hat_j;
      double C4 = C2 / n_vec[j];

      double post_p_j = 1 /
        (1 + inv_odd_p * ::sqrt(1 + C1) * ::exp(-C3 * C3 / C4 / 2));

      double diff = -curr_beta[j];
      if (sparse && (post_p_j < p)) {
        curr_beta[j] = 0;
      } else {
        if (post_p_j > ::unif_rand()) {
          curr_beta[j] = ::Rf_rnorm(C3, ::sqrt(C4));
          diff += curr_beta[j];
          gap += curr_beta[j] * curr_beta[j];
        } else {
          curr_beta[j] = 0;
        }
        if (k >= 0) avg_beta[j] += C3 * post_p_j;
      }
      if (diff != 0) dotprods = sfbm->incr_mult_col(j2, dotprods, diff);
    }

    if (gap > gap0) { avg_beta.fill(NA_REAL); return avg_beta; }
  }

  return avg_beta / num_iter;
}

/******************************************************************************/
