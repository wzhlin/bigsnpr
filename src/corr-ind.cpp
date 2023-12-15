/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include "bed-acc.h"

using namespace Rcpp;

/******************************************************************************/

template <class C>
List corMatInd0(C macc,
                const std::vector<size_t>& P,
                const std::vector<int>& I,
                int ncores) {

  int n = macc.nrow();
  int m = macc.ncol();

  List res(m);

  int chunk_size = ceil(m / (10.0 * ncores));

  #pragma omp parallel for schedule(dynamic, chunk_size) num_threads(ncores)
  for (int j0 = 0; j0 < m; j0++) {

    size_t lo = P[j0];
    size_t up = P[j0 + 1];
    if (lo == up) continue;

    size_t K = up - lo;
    NumericVector X(K);

    // pre-computation
    double xSum0 = 0, xxSum0 = 0;
    for (int i = 0; i < n; i++) {
      double x = macc(i, j0);
      if (x != 3) {
        xSum0  += x;
        xxSum0 += x * x;
      }
    }

    // main computation
    for (size_t k = 0; k < K; k++) {

      int j = I[lo + k];

      int nona = 0;
      double xSum = xSum0, xxSum = xxSum0;
      double ySum = 0, yySum = 0, xySum = 0;
      for (int i = 0; i < n; i++) {

        double x = macc(i, j0);
        if (x == 3) continue;

        double y = macc(i, j);
        if (y == 3) {
          // y is missing, but not x
          xSum  -= x;
          xxSum -= x * x;
        } else {
          // none missing
          nona++;
          ySum  += y;
          yySum += y * y;
          xySum += x * y;
        }
      }

      double num = xySum - xSum * ySum / nona;
      double deno_x = xxSum - xSum * xSum / nona;
      double deno_y = yySum - ySum * ySum / nona;
      double r = num / ::sqrt(deno_x * deno_y);
      if (r > 1) { r = 1; } else if (r < -1) { r = -1; }  // e.g. 1.000...04

      X[k] = r;
    }

    #pragma omp critical
    res[j0] = X;
  }

  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
List corMatInd(Environment obj,
               const IntegerVector& rowInd,
               const IntegerVector& colInd,
               const std::vector<size_t>& P,
               const std::vector<int>& I,
               int ncores) {

  myassert_size(colInd.size(), P.size() - 1);

  if (obj.exists("code256")) {
    XPtr<FBM> xpBM = obj["address"];
    NumericVector code = clone(as<NumericVector>(obj["code256"]));
    code[is_na(code)] = 3;
    SubBMCode256Acc macc(xpBM, rowInd, colInd, code, 1);
    return corMatInd0(macc,P, I, ncores);
  } else if (obj.exists("bedfile")) {
    XPtr<bed> xp_bed = obj["address"];
    bedAcc macc(xp_bed, rowInd, colInd);
    return corMatInd0(macc, P, I, ncores);
  } else {
    throw Rcpp::exception("Unknown object type.");
  }
}

/******************************************************************************/
