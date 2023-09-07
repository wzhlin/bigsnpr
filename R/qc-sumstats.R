################################################################################

# still have some problems


#' Eigen value decomposition to get sqrt root of inverse of x
#'
#' @param x A symmetric matrix
#' @param prop_eig Proportion of eigen values to be attained
#' @param eig_tol Min eigenvalue attained
#'
#' @return A pseudo sqrt root of inverse of x
#'
#' @example
#'
eigen_halfinv <- function(x, prop_eig, eig_tol = 1e-4) {

  e_object <- eigen(x, symmetric = TRUE) # do eigen decomp
  e_value <- e_object$values
  e_vector <- e_object$vector

  nrank <- sum(e_value > eig_tol)
  nkeep <- min(ceiling(prop_eig * nrow(x)), nrank) # get top prop_eig of eigen val

  ind_keep <- seq_len(nkeep)
  e_vector_keep <- e_vector[, ind_keep, drop = FALSE]
  e_value_keep <- e_value[ind_keep]

  sweep(e_vector_keep, 2, sqrt(e_value_keep), "/")
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
#' @example
#'
impute_z <- function(id_current,
                     id_high,
                     ld_current_high,
                     z_high,
                     prop_eig) {

  eig_scaled <- eigen_halfinv(ld[id_high, id_high], prop_eig = prop_eig) # get eig

  z_imputed <-
    (ld_current_high %*% eig_scaled) %*% crossprod(eig_scaled, z_high)

  deno <- 1 - tcrossprod(ld_current_high %*% eig_scaled)

  c(z_imputed, deno)
}


#' Quality Control of summary statistics
#'
#' @param z_scores z_scores from GWAS summary statistics
#' @param ld LD matrix
#' @param thr_highld Threshold for highly correlated variants, `0.05` by default
#' @param num_highld Minimum number of highly correlated variants, `20` by default
#' @param prop_eig Proportion of eigenvalues to attain, `0.4` by default
#' @param max_run Maximum round of iterations to run.
#' @param ncores Number of cores to run in parallel
#'
#' @return A dataframe of original z-score, imputed z-score, chisq statistics, p-value, and whether to be removed or not
#'
#' @export
#'
#' @example
#'

snp_qc_sumstats <- function(z_scores,
                            ld,
                            thr_highld = 0.05,
                            num_highld = 20,
                            prop_eig = 0.4,
                            max_run = 0.1 * length(z_scores),
                            ncores = 1) {

  chi2_gwide_thr <- qchisq(5e-8, df = 1, lower.tail = FALSE)

  all_imputed_z <- rep(NA_real_, length(z_scores))
  all_deno      <- rep(NA_real_, length(z_scores))
  all_id_high   <- list()

  removed <- rep(FALSE, length(z_scores))

  id_this <- seq_len(length(z_scores))

  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  st_time_all <- Sys.time()

  for (i_run in seq_len(max_run)) {

    st_time_this <- Sys.time()

    id_rm <- which(removed)

    cat("=== i_run =", i_run, "-- length for this:", length(id_this), "===\n")
    all_res_this <- foreach(id_current = id_this) %dopar% {
      ld_current <- ld[, id_current]
      ld_current[c(id_current, id_rm)] <- 0

      id_high <- find_highld(ld_current@i, ld_current@x, thr_highld, num_highld)

      res_impute <- impute_z(id_current, id_high,
                             ld_current_high = ld_current[id_high],
                             z_high = z_scores[id_high],
                             prop_eig = prop_eig)
      list(res_impute, id_high)
    }

    all_imputed_z[id_this] <- sapply(all_res_this, function(x) x[[1]][[1]])
    all_deno[id_this]      <- pmax(sapply(all_res_this, function(x) x[[1]][[2]]),
                                   .Machine$double.eps)
    all_id_high[id_this]   <- lapply(all_res_this, function(x) x[[2]])

    timing <- difftime(Sys.time(), st_time_this, units = "secs")

    cat("Time for this:", round(timing, 1), "seconds.\n")

    chi2 <- (z_scores - all_imputed_z)^2 / all_deno
    id_rm_this <- which.max(chi2 * !removed)

    if (chi2[id_rm_this] > chi2_gwide_thr) {
      removed[id_rm_this] <- TRUE
      run_next <- sapply(all_id_high, function(vec) {
        any(id_rm_this %in% vec)
      })
      id_this <- which(run_next)
    } else break

  } # end of iteration

  timing_all <- difftime(Sys.time(), st_time_all, units = "mins")
  cat("Time for all:", round(timing_all, 1), "minutes.\n")

  stopCluster(cl)

  chi2_final <- (z_scores - all_imputed_z)^2 / all_deno
  pval_final <- pchisq(chi2_final, df = 1, lower.tail = FALSE)


  data.frame(
    z_scores = z_scores,
    z_imputed = all_imputed_z,
    chi2_final = chi2_final,
    p_val = pval_final,
    remove = pval_final < 5e-8
  )

}

################################################################################
