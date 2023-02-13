library(dplyr)
sim_blk <- function(n_total, n_blk) {
  if (n_total / n_blk < 3) {
    blks <- rep(floor(n_total / n_blk), n_blk)
    n_extra <- n_total - sum(blks)
    blks[1 : n_extra] <- blks[1 : n_extra] + 1
    return(blks)
  }
  blks <- 0
  while(any(blks <= 0)) {
    fracs <- runif(n_blk)
    fracs <- fracs / sum(fracs)
    blks <- floor(n_total * fracs)
    # print(blks)
  }
  blks <- c(blks[-length(blks)], n_total - sum(blks[-length(blks)]))
  return(blks)
}
simulate_ld <- function(J, n_blk, rank_frac = 0.1) {
  # rank_frac controls the amount of LD. 
  # The rank of an LD block is equal to size of the block x rank_frac
  # 
  r <- round(rank_frac * J)
  J_blks <- sim_blk(J, n_blk)  # J's (number of snps) per block
  r_blks <- sim_blk(r, n_blk)  # r's (rank) per block
  # Z <- matrix(rnorm(N * r), nrow = N)
  Rs <- list()
  for(i in 1 : n_blk) {
    Rs[[i]] <- matrix(runif(J_blks[i] * r_blks[i], -0.5, 0.5), nrow = J_blks[i])
  }
  R <- do.call(Matrix::bdiag, Rs) %>% as.matrix
  R <- R + matrix(runif(J * r, -0.1, 0.1), nrow = J)  # add additional long range LD
  LD <- R %*% t(R)
  sd_diag_ld <- sqrt(diag(LD))
  LD <- sweep(LD, 1, sd_diag_ld, '/')
  LD <- sweep(LD, 2, sd_diag_ld, '/')
  return(list(LD = LD, rank = r, blk.size = J_blks))
}
simulate_genotype <- function(N, LD) {
  r <- LD$rank
  kk <- eigen(LD$LD)
  L <- kk$vectors[, 1 : r]
  v <- kk$values[1 : r]
  L <- sweep(L, 2, sqrt(v), '*')
  Z <- matrix(rnorm(r * N), nrow = N)
  return(Z %*% t(L))
}
# simulate_snp_effects <- function(M) {
#   return(rnorm(M))
# }
get_error_var <- function(y, h2) {
  return(var(y) / h2 * (1 - h2))
}
# simulate_phenotype <- function(genotypes, beta, h2) {
#   N <- nrow(genotypes)
#   y <- genotypes %*% beta
#   evar <- get_error_var(y, h2)
#   y <- y + rnorm(N, sd = sqrt(evar))
#   return(y)
# }
# simulate_mediation_phenotype <- function(mediator_ys, b, pve) {
#   y <- mediator_ys %*% b 
#   evar <- get_error_var(y, pve)
#   y <- y + rnorm(nrow(mediator_ys), sd = sqrt(evar))
#   return(y)
# } 
decomp_ld <- function(geno) {
  tmp <- svd(geno)
  m <- MASS::Null(tmp$v)
  ld_eig_vec <- cbind(tmp$v, m)
  ld_eig_val <- c(tmp$d, rep(0, ncol(m)))
  return(list(ld_eig_vec = ld_eig_vec, ld_eig_val = ld_eig_val))
}