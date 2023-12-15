/******************************************************************************/

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

inline double square(double x) { return x * x; }

/******************************************************************************/

// [[Rcpp::export]]
List find_highld(Environment corr,
                 int j,
                 LogicalVector& keep,
                 const NumericVector& thr) {

  bool keep_save = keep[j];
  keep[j] = false;  // discard itself

  Rcpp::XPtr<SFBM> sfbm = corr["address"];
  if (sfbm->is_compact()) Rcpp::stop("This does not work with a compact SFBM");
  const NumericVector p = corr["p"];
  const double * data = sfbm->i_x();

  std::vector<int>    highld_inds;
  std::vector<double> highld_vals;

  size_t lo = 2 * p[j];
  size_t up = 2 * p[j + 1];

  for (size_t k = lo; k < up; k += 2) {
    int i_k = data[k];
    if (keep[i_k]) {
      double x_k = data[k + 1];
      if (square(x_k) > thr[i_k]) {
        highld_inds.push_back(i_k);
        highld_vals.push_back(x_k);
      }
    }
  }

  keep[j] = keep_save;  // reset

  IntegerVector inds = wrap(highld_inds);
  NumericVector vals = wrap(highld_vals);

  return List::create(inds, vals);
}

/******************************************************************************/

// [[Rcpp::export]]
void test_ld_score(Environment corr,
                   const IntegerVector& ord,
                   const IntegerVector& ind,
                   LogicalVector& keep,
                   const NumericVector& thr) {

  Rcpp::XPtr<SFBM> sfbm = corr["address"];
  if (sfbm->is_compact()) Rcpp::stop("This does not work with a compact SFBM");
  const NumericVector p = corr["p"];
  const double * data = sfbm->i_x();

  for (auto& j : ind) keep[j] = true;

  for (auto& o : ord) {

    int j = ind[o];
    size_t lo = 2 * p[j];
    size_t up = 2 * p[j + 1];

    double ld_score = 0;

    for (size_t k = lo; k < up; k += 2) {
      int i_k = data[k];
      if (keep[i_k]) {
        double x_k = data[k + 1];
        ld_score += square(x_k);
        if (ld_score > thr[o]) {
          keep[j] = false;
          break;
        }
      }
    }
  }
}

/******************************************************************************/
