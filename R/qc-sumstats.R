################################################################################

#' Eigen value decomposition to get sqrt root of inverse of x
#'
#' @param x A symmetric matrix
#' @param prop_eig Proportion of eigen values to be attained
#' @param eig_tol Min eigenvalue attained
#'
#' @return A pseudo sqrt root of inverse x
#'
#' @examples
#'
#' x <- crossprod(matrix(1:4, 2))
#' x_inv_half <- bigsnpr:::eigen_halfinv(x, prop_eig = 1)
#' all.equal(solve(x), tcrossprod(x_inv_half))
#'
eigen_halfinv <- function(x, prop_eig, eig_tol = 1e-4) {

  e_object <- eigen(x, symmetric = TRUE) # do eigen decomp
  e_value <- e_object$values
  e_vector <- e_object$vectors

  nrank <- sum(e_value > eig_tol)
  nkeep <- min(ceiling(prop_eig * nrow(x)), nrank) # get top prop_eig of eigen val

  ind_keep <- seq_len(nkeep)
  e_vector_keep <- e_vector[, ind_keep, drop = FALSE]
  e_value_keep <- e_value[ind_keep]

  sweep(e_vector_keep, 2, sqrt(e_value_keep), "/")
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

#' Impute z-score using highly correlated variants
#'
#' @param id_current current id of variant to be impute
#' @param id_high id of highly correlated variants
#' @param ld_current_high LD matrix for id_current and id_high
#' @param z_high z-scores for highly correlated variants
#' @param prop_eig proportion of eigenvalues to be used in eigenvalue decomp
#'
#' @return Imputed z-score, and denominator of chisq statistics
#'
#' @examples
#' obj_bigsnp <- snp_attachExtdata()
#' G <- obj_bigsnp$genotypes
#' Z <- big_univLogReg(G, obj_bigsnp$fam$affection - 1, ind.col = 1:500)$score
#' ld <- snp_cor(G, ind.col = 1:500)
#' id_current <- which.max(Matrix::colSums(ld^2))
#' ld_current <- ld[, id_current]
#' id_high <- setdiff(which(ld_current^2 > 0.05), id_current)
#' ld_high_high <- ld[id_high, id_high]
#' z_imputed <- bigsnpr:::impute_z(ld_high_high = ld_high_high,
#'                                 ld_current_high = ld_current[id_high],
#'                                 z_high = Z[id_high],
#'                                 prop_eig = 0.4)
#' chi2 <- (Z[id_current] - z_imputed[1, ])^2 / z_imputed[2, ]
#' pchisq(chi2, df = 1, lower.tail = FALSE)
#'
#'
impute_z <- function(ld_high_high,
                     ld_current_high,
                     z_high,
                     prop_eig = 0.4) {

  eig_scaled <- eigen_halfinv(x = ld_high_high, prop_eig = prop_eig) # get eig

  r_half <- ld_current_high %*% eig_scaled

  z_imputed <- cumsum(r_half * as.vector(crossprod(eig_scaled, z_high)))

  deno <- 1 - cumsum(r_half^2)

  rbind(z_imputed, deno)
}


#' Quality Control of summary statistics
#'
#' @param z_sumstats z_sumstats from GWAS summary statistics
#' @param corr LD matrix
#' @param thr_highld Threshold for highly correlated variants, `0.05` by default
#' @param num_highld Minimum number of highly correlated variants, `20` by default
#' @param prop_eig Proportion of eigenvalues to attain, `0.4` by default
#' @param max_run Maximum round of iterations to run.
#' @param p_val_thr P-value Threshold for QC, `5e-8` by default
#' @param print_info Print_info or not, `FALSE` by default
#' @param ncores Number of cores to run in parallel
#'
#' @export
#'
#' @return A dataframe wth three variables of each variants:
#' - `$chi2_final`: chi-squared statistics
#' - `$p_val`: corresponding p-value
#' - `$remove`: corresponding variant being removed or not
#'
#'
#' @examples
#'
#'
#' obj_bigsnp <- snp_attachExtdata()
#' G <- obj_bigsnp$genotypes
#' Z <- big_univLogReg(G, obj_bigsnp$fam$affection - 1, ind.col = 1:500)$score
#' ind <- which.max(abs(Z))
#' Z[ind] <- -Z[ind]
#' ld <- as_SFBM(snp_cor(G, ind.col = 1:500), compact = TRUE)
#' res <- snp_qc_sumstats(z_sumstats = Z, corr = ld)
#' cat(res$p_val[ind], res$remove[ind])
#'

snp_qc_sumstats <- function(z_sumstats,
                            corr,
                            thr_highld = 0.05,
                            num_highld = 20,
                            prop_eig = 0.4,
                            max_run = 0.1 * length(z_sumstats),
                            p_val_thr = 5e-8,
                            print_info = FALSE,
                            ncores = 1) {

  bigassertr::assert_lengths(z_sumstats, rows_along(corr), cols_along(corr))

  chi2_gwide_thr <- stats::qchisq(p_val_thr, df = 1, lower.tail = FALSE)

  all_chi2      <- rep(NA_real_, length(z_sumstats))
  all_id_high   <- list()

  removed <- rep(FALSE, length(z_sumstats))
  id_this <- seq_along(z_sumstats)

  bigparallelr::register_parallel(ncores)

  st_time_all <- Sys.time()

  for (i_run in seq_len(max_run)) {

    st_time_this <- Sys.time()

    id_rm <- which(removed)

    if (print_info) {
      cat("=== Iteration #", i_run, "length for this iteration:", length(id_this), "===\n")
    }

    if (length(id_this) > 0) {


      fun_use <- c("find_highld", "select_useful", "impute_z")
      all_res_this <- foreach(id_current = id_this, .export = fun_use) %dopar% {

        ld_current <- corr[, id_current]
        ld_current[id_current] <- 0
        ld_current[removed]    <- 0

        id_high <- find_highld(ld_current@i, ld_current@x, thr_highld, num_highld)


        if (num_highld > 0 && id_high[num_highld] < 0) {
          id_high <- id_high[id_high > 0]
        }
        if (length(id_high) < 2)
          return(list(NA_real_, id_high))

        ld_high_high <- corr[id_high, id_high]
        keep <-  select_useful(ld_high_high, ld_current[id_high])

        id_high <- id_high[keep]
        if (length(id_high) < 2)
          return(list(NA_real_, id_high))


        res_impute <-  impute_z(ld_high_high = corr[id_high, id_high],
                                ld_current_high = ld_current[id_high],
                                z_high = z_sumstats[id_high],
                                prop_eig = prop_eig)

        all_chi2 <- (z_sumstats[id_current] - res_impute[1, ])^2 /
          pmax(res_impute[2, ], 0.01)

        list(min(all_chi2[-1]), id_high)
      }

      all_chi2[id_this]    <- sapply(all_res_this, function(x) x[[1]])
      all_id_high[id_this] <- lapply(all_res_this, function(x) x[[2]])

      if (print_info) {
        timing <- difftime(Sys.time(), st_time_this, units = "secs")
        cat("Time for this iteration:", round(timing, 1), "seconds.\n")
      }

    }

    id_rm_this <- which.max(all_chi2 * !removed)

    if (all_chi2[id_rm_this] > chi2_gwide_thr) {
      removed[id_rm_this] <- TRUE
      cat("   ", id_rm_this, " will be removed.\n")

      run_next <- sapply(all_id_high, function(vec) {
        id_rm_this %in% vec
      })
      id_this <- which(run_next)
    } else break

  } # end iteration

  if (print_info) {
    timing_all <- difftime(Sys.time(), st_time_all, units = "mins")
    cat("Time for all iterations:", round(timing_all, 1), "minutes.\n")
  }

  pval_final <- stats::pchisq(all_chi2, df = 1, lower.tail = FALSE)

  cat("End for Quality Control.\n")

  data.frame(
    chi2_final = all_chi2,
    p_val = pval_final,
    remove = pval_final < p_val_thr
  )
}

################################################################################
