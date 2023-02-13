run_gwas <- function(phenotype, genotype) {
  res <- fast_linear_regression(as.vector(phenotype), genotype, rep(1, nrow(genotype)))
  return(do.call(cbind, res))
}

run_s_brainxcan <- function(betahat, se_betahat, ld_eig_vec, ld_eig_val, weights) {
  tmp <- ld_eig_vec %*% weights
  sigma2_g <- sum(tmp ^ 2 * ld_eig_val)
  sigma2_l <- get_ld_diag(ld_eig_vec, ld_eig_val)
  beta <- sum(weights * betahat * sigma2_l) / sigma2_g
  z <- sum(weights * sqrt(sigma2_l / sigma2_g) * betahat / se_betahat)
  return(list(beta = beta, z = z))
}

fast_linear_regression = function(y, x, covariate = NULL) {
  x = as.matrix(x)
  if(is.null(covariate)) {
    y_ = y
    x_ = x
    dof = length(y) - 1
  } else {
    covariate = as.matrix(covariate)
    res = qr(covariate)
    Q_ = qr.Q(res)
    x_ = x - Q_ %*% (t(Q_) %*% x)
    y_ = as.numeric(y - Q_ %*% (t(Q_) %*% y))
    dof = length(y) - ncol(covariate) - 1
  }
  bhat = colMeans(y_ * x_) / colMeans(x_ ^ 2)
  sigma2 = colSums((y_ - x_ * bhat) ^ 2) / dof
  se = sqrt(1 / colSums(x_ ^ 2) * sigma2 )
  pval = pt(abs(bhat / se), dof, lower.tail = F, log.p = T)
  list(bhat = bhat, pval = exp(pval) * 2, se = se)
}