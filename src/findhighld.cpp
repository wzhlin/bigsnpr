#include <Rcpp.h>
using namespace Rcpp;

inline double square(double x) { return x * x; }

// [[Rcpp::export]]
IntegerVector find_highld(const IntegerVector& I,
                          const NumericVector& X,
                          double thr_highld,
                          int min_nb_highld) {

  std::vector<int> highld_inds;

  int K = X.size();
  for (int k = 0; k < K; k++) {
    double X2_k = square(X[k]);
    if (X2_k > thr_highld)
      highld_inds.push_back(I[k] + 1);
  }

  if (int(highld_inds.size()) >= min_nb_highld) {

    IntegerVector top_inds = wrap(highld_inds);
    top_inds.attr("thresholded") = true;
    return top_inds;

  } else {

    IntegerVector top_inds(min_nb_highld, -1);
    NumericVector top_vals(min_nb_highld, -1.0);

    for (int k = 0; k < K; k++) {
      double X2_k = square(X[k]);
      int I_k = I[k] + 1;
      if (X2_k > top_vals[min_nb_highld-1]) {  // larger than at least one
        // place it at the end of the list
        int j = min_nb_highld-1;
        top_inds[j] = I_k;
        top_vals[j] = X2_k;
        j--;
        // check whether can move it further up
        for(; j >= 0 && X2_k > top_vals[j]; j--) {
          // move this one down the list
          top_inds[j+1] = top_inds[j];
          top_vals[j+1] = top_vals[j];
          // and move k up
          top_inds[j] = I_k;
          top_vals[j] = X2_k;
        }
      }
    }

    top_inds.attr("thresholded") = false;
    return top_inds;
  }

}
