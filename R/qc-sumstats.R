################################################################################

#' Eigen decomposition to get x^-0.5
#'
#' @param x A square symmetric matrix.
#' @param prop_eigs The proportion of eigen components to keep.
#' @param tol The minimum eigenvalue allowed; others are discarded.
#'
#' @return An approximation of x^-0.5.
#'
#' @examples
#' x <- tcrossprod(matrix(1:4, 2))
#' eig_scaled <- bigsnpr:::eigen_halfinv(x, prop_eigs = 1)
#' all.equal(solve(x), tcrossprod(eig_scaled))
#'
eigen_halfinv <- function(x, prop_eigs, tol = 1e-4) {

  e <- eigen(x, symmetric = TRUE)  # eigen decomposition
  e_value  <- e$values
  e_vector <- e$vectors

  nrank <- sum(e_value > tol)
  nkeep <- min(ceiling(prop_eigs * nrow(x)), nrank)

  ind_keep <- seq_len(nkeep)
  e_vector_keep <- e_vector[, ind_keep, drop = FALSE]
  e_value_keep <- e_value[ind_keep]

  sweep(e_vector_keep, 2, sqrt(e_value_keep), '/')
}

################################################################################

select_useful <- function(ld_high, ld_current_high, thr = 0.01) {

  w <- ld_current_high^2
  keep <- rep(TRUE, length(w))
  ld_score <- Matrix::colSums(ld_high^2)

  for (j in order(w)) {
    if ((w[j]^sqrt(ld_score[j])) < thr) {
      keep[j] <- FALSE
      ld_score <- update_ldscore(ld_score, ld_high@i, ld_high@x, ld_high@p, j - 1)
      # in the end, ld_score should be 0 for those not kept
    }
  }

  keep
}

################################################################################

#' Impute z-score using nearby variants
#'
#' @param ld_current_high LD between variant to impute and its high LD variants.
#' @param ld_high_high LD matrix between the high LD variants.
#' @param z_high z-scores of `id_high` variants.
#' @inheritParams eigen_halfinv
#'
#' @return The imputed z-score, as well as the denominator of the chi-squared
#'   statistic (z - z_imp)^2 / deno.
#'
#' @references Chen, Wenhan, et al. "Improved analyses of GWAS summary statistics
#'   by reducing data heterogeneity and errors." Nature Communications 12.1
#'   (2021): 7117. \doi{10.1038/s41467-021-27438-7}.
#'
#' @examples
#' bigsnp <- snp_attachExtdata()
#' G <- bigsnp$genotypes
#' Z <- big_univLogReg(G, bigsnp$fam$affection - 1, ind.col = 1:500)$score
#' ld <- snp_cor(G, ind.col = 1:500)
#' id_current <- which.max(Matrix::colSums(ld^2))
#' ld_current <- ld[, id_current]
#' id_high <- setdiff(which(ld_current^2 > 0.02), id_current)
#' z_imp <- bigsnpr:::impute_z(ld_current_high = ld_current[id_high],
#'                             ld_high_high = ld[id_high, id_high],
#'                             z_high = Z[id_high],
#'                             prop_eigs = 0.4)
#' c(Z[id_current], z_imp[1, ])
#' chi2 <- (Z[id_current] - z_imp[1, ])^2 / z_imp[2, ]
#' pchisq(chi2, df = 1, lower.tail = FALSE)
#'
impute_z <- function(ld_current_high, ld_high_high, z_high, prop_eigs) {

  # eigendecomposition to get R_tt^-0.5
  eig_scaled <- eigen_halfinv(ld_high_high, prop_eigs = prop_eigs)

  # R_it R_tt^-0.5
  r_half <- ld_current_high %*% eig_scaled

  # R_it R_tt^-1 z_t
  z_imputed <- cumsum(r_half * as.vector(crossprod(eig_scaled, z_high)))

  # 1 - R_it R_tt^-1 R_ti
  deno <- 1 - cumsum(r_half^2)

  rbind(z_imputed, deno)
}

################################################################################

#' Quality control of summary statistics
#'
#' Quality control of summary statistics, by comparing z-scores with z-scores
#' imputed using nearby correlated variants.
#'
#' @inheritParams snp_ldpred2_auto
#' @param z_sumstats Initial z-scores from GWAS summary statistics.
#' @param thr_highld Threshold on squared correlations to decide which variants
#'   to use for imputation. Default is `0.05`.
#' @param min_nb_highld Minimum number of variants to use for imputation.
#'   In case not enough variants pass `thr_highld`, the most correlated are used.
#'   Default is `20`.
#' @inheritParams eigen_halfinv
#' @param max_run Maximum number of iterations (and therefore the maximum number
#'   of variants that can be removed). Default is 10% of all variants.
#' @param pval_thr P-value threshold for QCing. Default is `5e-8`.
#' @param print_iter Whether to print the iteration numbers. Default is `FALSE`.
#'
#' @return A tibble (data frame) with 6 variables for each variants:
#'   - `$z_ss`: the initial vector of `z_sumstats`,
#'   - `$z_ss_imp`: the imputed z-scores,
#'   - `$deno_qc`: the denominators of the chi-squared statistics for QCing,
#'   - `$chi2_qc`: the chi-squared statistics,
#'   - `$pval_qc`: the corresponding p-values,
#'   - `$rm_qc`: whether the p-values are significant or not.
#' @export
#'
#' @examples
#' bigsnp <- snp_attachExtdata()
#' G <- bigsnp$genotypes
#' Z <- big_univLogReg(G, bigsnp$fam$affection - 1, ind.col = 1:500)$score
#' ind <- which.max(abs(Z)); Z[ind] <- -Z[ind]
#' ld <- snp_cor(G, ind.col = 1:500)
#' corr <- as_SFBM(ld, compact = TRUE)
#'
#' res <- snp_qc_sumstats(corr, Z)
#'
#' res$pval_qc[ind]
#' summary(res$pval_qc[-ind])
#'
snp_qc_sumstats <- function(corr, z_sumstats,
                            thr_highld = 0.05,
                            min_nb_highld = 20,
                            prop_eigs = 0.4,
                            max_run = 0.1 * length(z_sumstats),
                            pval_thr = 5e-8,
                            print_iter = FALSE,
                            ncores = 1) {

  assert_lengths(rows_along(corr), cols_along(corr), z_sumstats)
  assert_nona(z_sumstats)

  chi2_thr <- stats::qchisq(pval_thr, df = 1, lower.tail = FALSE)

  # preallocate results
  all_chi2    <- rep(NA_real_, length(z_sumstats))
  all_id_high <- list()

  # for the first iteration, impute all, using all other variants
  removed <- rep(FALSE, length(z_sumstats))
  id_this <- seq_along(z_sumstats)

  bigparallelr::register_parallel(ncores)

  for (i_run in seq_len(max_run)) {

    if (print_iter) cat(i_run, "")

    FUNs <- c("find_highld", "select_useful", "impute_z")
    all_res_this <- foreach(id_current = id_this, .export = FUNs) %dopar% {

      ld_current <- corr[, id_current]
      ld_current[id_current] <- 0
      ld_current[removed]    <- 0

      id_high <- find_highld(ld_current@i, ld_current@x, thr_highld, min_nb_highld)
      if (min_nb_highld > 0 && id_high[min_nb_highld] < 0)
        id_high <- id_high[id_high > 0]
      if (length(id_high) < 2) return(list(NA_real_, id_high))

      keep <- select_useful(corr[id_high, id_high], ld_current[id_high])

      id_high <- id_high[keep]
      if (length(id_high) < 2) return(list(NA_real_, id_high))

      res_impute <- impute_z(ld_current_high = ld_current[id_high],
                             z_high = z_sumstats[id_high],
                             ld_high_high = corr[id_high, id_high],
                             prop_eigs = 1)

      all_chi2 <- (z_sumstats[id_current] - res_impute[1, ])^2 /
        pmax(res_impute[2, ], 0.01)

      list(min(all_chi2[-1]), id_high)
    }

    # save current results
    all_chi2[id_this]    <- sapply(all_res_this, function(x) x[[1]])
    all_id_high[id_this] <- lapply(all_res_this, function(x) x[[2]])

    # get the worst variant and remove it
    id_worst <- which.max(all_chi2 * !removed)
    if (all_chi2[id_worst] > chi2_thr) {
      removed[id_worst] <- TRUE
      # find which variants used id_worst for imputation, and update them
      run_next <- sapply(all_id_high, function(id_high) id_worst %in% id_high)
      id_this <- which(run_next)
    } else break

  }

  pval <- stats::pchisq(all_chi2, df = 1, lower.tail = FALSE)
  # note that a variant that has been removed in one iteration can still have
  # a non-significant p-value in the end; this prevents many false positives

  tibble::tibble(
    chi2_qc  = all_chi2,
    pval_qc  = pval,
    rm_qc    = pval < pval_thr
  )
}

################################################################################
