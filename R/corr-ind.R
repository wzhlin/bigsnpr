################################################################################

corInd0 <- function(Gna, list_ind,
                    ind.row = rows_along(Gna),
                    ind.col = cols_along(Gna),
                    ncores = 1) {

  assert_lengths(list_ind, ind.col)

  list_ind <- lapply(seq_along(list_ind), function(k) {
    ind <- list_ind[[k]]
    unique(sort(ind[k >= ind])) - 1L
  })

  corr <- new("dsCMatrix", uplo = "U")
  m <- length(list_ind)
  corr@Dim <- c(m, m)
  corr@i <- unlist(list_ind)
  corr@p <- c(0L, cumsum(lengths(list_ind)))
  rm(list_ind); gc()

  corr@x <- unlist(corMatInd(
    obj    = Gna,
    rowInd = ind.row,
    colInd = ind.col,
    P      = corr@p,
    I      = corr@i,
    ncores = ncores
  ))

  if (anyNA(corr@x))
    warning2("NA or NaN values in the resulting correlation matrix.")

  corr
}

################################################################################

#' Correlation matrix
#'
#' TODO: CHANGE DESC AND DOCU.
#' Get significant (Pearson) correlations between nearby SNPs of the same chromosome
#' (p-values are computed using a two-sided t-test).
#'
#' @inheritParams bigsnpr-package
#' @param size For one SNP, window size around this SNP to compute correlations.
#' Default is `500`. If not providing `infos.pos` (`NULL`, the default), this is
#' a window in number of SNPs, otherwise it is a window in kb (genetic distance).
#' @param alpha Type-I error for testing correlations.
#'   Default is `1` (no threshold is applied).
#' @param thr_r2 Threshold to apply on squared correlations. Default is `0`.
#' @param fill.diag Whether to fill the diagonal with 1s (the default)
#' or to keep it as 0s.
#'
#' @return The (Pearson) correlation matrix. This is a sparse symmetric matrix.
#'
#' @import Matrix
#'
#' @examples
#' test <- snp_attachExtdata()
#' G <- test$genotypes
#'
#' keep <- Matrix::rsparsematrix(100, 100, 0.05, symmetric = TRUE) != 0
#' ind <- Matrix::which(keep, arr.ind = TRUE)
#' list_ind <- split(ind[, 1], factor(1:100)[ind[, 2]])
#'
#' corr <- snp_corInd(G, ind.col = 1:100, list_ind = list_ind)
#' ind2 <- Matrix::which(corr != 0, arr.ind = TRUE)
#' all.equal(ind, ind2)
#'
#' @export
#'
snp_corInd <- function(Gna, list_ind,
                       ind.row = rows_along(Gna),
                       ind.col = cols_along(Gna),
                       ncores = 1) {

  args <- as.list(environment())

  check_args()

  do.call(corInd0, args)
}

################################################################################

#' @rdname snp_corInd
#' @export
bed_corInd <- function(obj.bed, list_ind,
                       ind.row = rows_along(obj.bed),
                       ind.col = cols_along(obj.bed),
                       ncores = 1) {

  args <- as.list(environment())
  names(args)[names(args) == "obj.bed"] <- "Gna"

  check_args()

  do.call(corInd0, args)
}

################################################################################
