library(bigsnpr)
bigsnp <- snp_attachExtdata()
G <- bigsnp$genotypes
Z <- big_univLogReg(G, bigsnp$fam$affection - 1, ind.col = 1:500)$score
ld <- snp_cor(G, ind.col = 1:500)
id_current <- which.max(Matrix::colSums(ld^2))
ld_current <- ld[, id_current]
id_high <- setdiff(which(ld_current^2 > 0.02), id_current)
ld_current_high = ld_current[id_high]
ld_high_high = ld[id_high, id_high]

id_high2 <- c(id_current, id_high)
ld_high2 <- ld[id_high2, id_high2]
eig_scaled <- bigsnpr:::eigen_halfinv(ld_high2, prop_eigs = 1)

# these three are the same
1 - ld_current_high %*% solve(ld_high_high) %*% ld_current_high
1 / solve(ld_high2)[1, 1]
1 / crossprod(eig_scaled[1, ])
