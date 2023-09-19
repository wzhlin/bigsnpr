
################################################################################

context("qc")

################################################################################

test_that("this is a test", {

  # local host at 10.0.0.158
  library(bigsnpr)
  ncores <- nb_cores()*2    # for parallel computing
  id_chr <- 22
  prop_eig <- 0.4
  thr_highld <- 0.05      # change thr for ld         # change
  prop_error <- 0.01      # error proportion          # change
  num_highld <- 20        # min number of high ld
  max_run <- 500          # max run for while loop

  setwd("C:/Users/linwe/Desktop/dentist-main")
  obj_bigsnp <- readRDS("C:/Users/linwe/Desktop/dentist-main/1000g/chr_22.rds")
  geno <- obj_bigsnp$genotypes
  pos <- obj_bigsnp$map$physical.pos
  chr <- obj_bigsnp$map$chromosome

  ld <- runonce::save_run(
    as_SFBM(as(readRDS(paste0("C:/Users/linwe/Desktop/dentist-main/1000g/ld.chr_", id_chr, ".rds")), "generalMatrix"),
            backingfile = paste0("C:/Users/linwe/Desktop/dentist-main/1000g/ld_sfbm.chr_", id_chr), compact = TRUE),
    file = paste0("C:/Users/linwe/Desktop/dentist-main/1000g/ld_sfbm.chr_", id_chr, ".rds")
  )


    # generate phenotype
  causal_rate  <- 0.1
  m_causal <- round(ncol(geno) * causal_rate)
  set.seed(123); causal_pos <- sort(sample(cols_along(geno), m_causal))
  



  # set.seed(100); ind_wrong <- sort(sample(causal_pos, 0.1 * length(causal_pos)))
  set.seed(100); ind_wrong <- sort(sample(sample(cols_along(geno), m_causal), 0.1 * length(causal_pos)))
  
  geno_wrong <- big_copy(geno)
  temp <- sapply(ind_wrong, function(i) {
    geno_wrong[, i] <- sample(3, nrow(geno_wrong), replace = TRUE) - 1
    0
  })

  # bigassertr::assert_dir("tmp-data")
  # pos2 <- snp_asGeneticPos(chr, pos, dir = "tmp-data")
  # ld_wrong <- snp_cor(geno_wrong, size = 3/1000, infos.pos = pos2, ncores = ncores)
  # ld_wrong <- as_SFBM(ld_wrong, compact = TRUE)





  # sumstat <- bigreadr::fread2("C:/Users/linwe/Desktop/dentist-main/1000g/sumstat.chr_22.txt")


  h2  <- 0.5
  prob <- numeric(ncol(geno))
  prob[causal_pos] <- 1
  obj_pheno_wrong <- snp_simuPheno(G = geno_wrong, h2 = h2, M = m_causal, prob = prob)
  y_wrong <- obj_pheno_wrong$pheno
  # gwas summary statistics
  gwas <- big_univLinReg(X = geno_wrong, y.train = y_wrong, ncores = ncores)
  p_value <- pchisq((gwas$estim / gwas$std.err)^2, df = 1, lower.tail = FALSE)

  # allele freq
  af <- big_colstats(geno_wrong, ncores = ncores)$sum / (2 * nrow(geno_wrong))

  sumstat_wrong <- data.frame(SNP = obj_bigsnp$map$marker.ID,
                          A1 = obj_bigsnp$map$allele1,
                          A2 = obj_bigsnp$map$allele2,
                          freq = af,
                          beta = gwas$estim,
                          se = gwas$std.err,
                          p = p_value,
                          N = nrow(geno_wrong))



  sumstat_now <- sumstat_wrong
  ld_run <- ld

  z_sumstats <- with(sumstat_now, beta / se)


  res <- snp_qc_sumstats(z_sumstats = z_sumstats, ld = ld_run,
                         thr_highld = thr_highld, prop_eig = prop_eig,
                         max_run = 500, ncores = ncores, print_info = TRUE)

  chi2_final <- res$chi2_final

  removed <- res$remove
  ind_wrong_logic <- numeric(ncol(geno_wrong))
  ind_wrong_logic[ind_wrong] <- 1
  df <- data.frame(ind_wrong = as.logical(ind_wrong_logic),
                    removed = as.logical(removed))
  print(with(df, table(removed, ind_wrong, exclude = NULL)))

# in selected causal snps
#          ind_wrong
# removed FALSE  TRUE
#   FALSE 17285   114
#   TRUE      0    60

# in all snps
#        ind_wrong
# removed FALSE  TRUE
#   FALSE 17223    95
#   TRUE     62    79


  expect_equal(1,1, check.attributes = FALSE, tolerance = 4e-5)

})
