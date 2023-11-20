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
                            removed = rep(FALSE, length(z_sumstats)),
                            max_run = length(z_sumstats),
                            pval_thr = 5e-8,
                            print_iter = FALSE,
                            thr_select = 0.01,
                            ncores = 1) {

  assert_lengths(rows_along(corr), cols_along(corr), z_sumstats, removed)
  assert_nona(z_sumstats)

  chi2_thr <- stats::qchisq(pval_thr, df = 1, lower.tail = FALSE)

  # preallocate results
  all_chi2    <- rep(-1, length(z_sumstats))
  all_id_high <- list()

  keep_high <- rep(FALSE, length(z_sumstats))  # placeholder

  # for the first iteration, impute all, using all other variants
  id_this <- seq_along(z_sumstats)
  keep <- !removed

  bigparallelr::register_parallel(ncores)

  for (i_run in seq_len(max_run)) {

    if (print_iter) cat(i_run, "")

    if (length(id_this) > 0) {

      FUNs <- c("find_highld", "test_ld_score")
      all_res_this <- foreach(id_current = id_this, .export = FUNs) %dopar% {

        high_ld <- find_highld(corr, id_current - 1L, keep, thr_highld)
        ind_high_r <- high_ld[[1]] + 1L
        if (length(ind_high_r) < 2) return(list(NA_real_, ind_high_r))

        # select useful -> test for w^sqrt(ld_score) < thr
        r <- pchisq(all_chi2[ind_high_r], df = 20, lower.tail = FALSE)
        w <- high_ld[[2]]^2 * r  # use r to prioritize better variants
        test_ld_score(corr, ord = order(w) - 1L, ind = high_ld[[1]],
                      keep = keep_high, thr = (log(thr_select) / log(w))^2)
        id_high <- which(keep_high)
        keep_high[id_high] <- FALSE  # reset

        ld_current_high <- high_ld[[2]][match(id_high, ind_high_r)]
        ld_high_high <- as.matrix(corr[id_high, id_high])

        e <- eigen(ld_high_high, symmetric = TRUE)  # eigen decomposition
        e_val  <- e$values
        ind_keep <- which(e_val > (e_val[1] * 1e-10))
        e_vec_keep <- e$vectors[, ind_keep, drop = FALSE]
        e_val_keep <- e_val[ind_keep]

        # pick some regularization (add to the diagonal)
        r_half_noscale_sq <- drop(ld_current_high %*% e_vec_keep)^2
        if (sum(r_half_noscale_sq / e_val_keep) < 0.99) {
          lam <- 0
        } else {
          lam <- uniroot(interval = c(0, 10), extendInt = "no", function(lam)
            0.99 - sum(r_half_noscale_sq / (e_val_keep + lam)))$root
        }
        # R_tt^-0.5
        eig_scaled <- sweep(e_vec_keep, 2, sqrt(e_val_keep + lam), '/')
        # R_it R_tt^-0.5
        r_half <- drop(ld_current_high %*% eig_scaled)
        # R_it R_tt^-1 z_t
        z_imputed <- cumsum(r_half * crossprod(eig_scaled, z_sumstats[id_high]))
        # 1 - R_it R_tt^-1 R_ti
        deno <- 1 - cumsum(r_half^2)

        multi_chi2 <- (z_sumstats[id_current] - z_imputed)^2 / deno

        list(median(head(multi_chi2, length(multi_chi2) / 2)), id_high)
      }

      # save current results
      all_chi2[id_this]    <- sapply(all_res_this, function(x) x[[1]])
      all_id_high[id_this] <- lapply(all_res_this, function(x) x[[2]])
    }

    # get the worst variant and remove it
    id_worst <- which.max(all_chi2 * keep)
    if (all_chi2[id_worst] > chi2_thr) {
      keep[id_worst] <- FALSE
      # find which variants used id_worst for imputation, and update them
      if (i_run > 1) {  # redo them all in the 2nd iter, for updating 'r'
        run_next <- sapply(all_id_high, function(id_high) id_worst %in% id_high)
        id_this <- which(run_next)
      }
    } else break

  }

  pval <- stats::pchisq(all_chi2, df = 1, lower.tail = FALSE)
  # note that a variant that has been removed in one iteration can still have
  # a non-significant p-value in the end; this prevents many false positives

  tibble::tibble(
    chi2_qc = all_chi2,
    pval_qc = pval,
    rm_qc   = (pval < pval_thr) | removed
  )
}

################################################################################
