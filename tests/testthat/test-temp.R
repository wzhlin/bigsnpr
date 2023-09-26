
################################################################################

context("qc")

################################################################################

test_that("this is a test", {

  # local host at 10.0.0.158
  library(bigsnpr)
  library(dplyr)
  ncores <- 12    # for parallel computing
  id_chr <- 22
  prop_eig <- 0.4
  thr_highld <- 0.05      # change thr for ld         # change
  prop_error <- 0.01      # error proportion          # change
  num_highld <- 20        # min number of high ld
  max_run <- 500          # max run for while loop

  ld <- runonce::save_run(
    as_SFBM(as(readRDS(paste0("C:/Users/linwe/Desktop/dentist-main/1000g/ld.chr_", id_chr, ".rds")), "generalMatrix"),
            backingfile = paste0("C:/Users/linwe/Desktop/dentist-main/1000g/ld_sfbm.chr_", id_chr), compact = TRUE),
    file = paste0("C:/Users/linwe/Desktop/dentist-main/1000g/ld_sfbm.chr_", id_chr, ".rds")
  )

  # read or run error
  sumstat_error <- bigreadr::fread2(paste0(
    "C:/Users/linwe/Desktop/dentist-main/sumstat/sumstat.error.", prop_error, "_chr_", id_chr, ".txt"))

  z_sumstats <- with(sumstat_error, beta / se)

  res <- snp_qc_sumstats(z_sumstats = z_sumstats, ld = ld,
                         thr_highld = thr_highld, prop_eig = prop_eig,
                         print_info = TRUE,
                         max_run = 500, ncores = ncores)

  chi2_final <- res$chi2_final
  sumstat <- bigreadr::fread2(paste0("C://Users/linwe/Desktop/dentist-main/1000g/sumstat.chr_22.txt"))
  
  df_res <- data.frame(
    is_error = !near(sumstat$beta, sumstat_error$beta),
    qc = res$remove
  )
  
  print(with(df_res, table(is_error, qc)))
#           qc
# is_error FALSE  TRUE
#    FALSE 17268    11
#    TRUE    140    34
  
  ind_gwas_low <- which(sumstat_error$p < 0.01)
  
  print(with(df_res[ind_gwas_low,], table(is_error, qc)))
#           qc
# is_error FALSE TRUE
#    FALSE  1046    0
#    TRUE      2    7
  

  expect_equal(chi2_final[1], 0.7524024, check.attributes = FALSE, tolerance = 4e-5)

})
