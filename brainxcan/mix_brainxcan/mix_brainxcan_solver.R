# main function
solve_mixed_brainxcan <- function(betahat, se_betahat, gwas_n, ld_eig_vec, ld_eig_val, weights,
                       tol = 1e-5, max_iter = 5000) {
  d <- get_ld_diag(ld_eig_vec, ld_eig_val)
  tbeta <- t(ld_eig_vec) %*% (betahat * d)
  sigma2_Y <- gwas_n * mean(d * se_betahat ^ 2)
  a <- sigma2_Y
  b <- 1e-6
  beta <- 0
  tb <- t(ld_eig_vec) %*% weights
  diff <- 1
  niter <- 0
  tD <- calc_tD(b, gwas_n, ld_eig_val)
  # diffs <- c()
  # as <- c()
  # bs <- c()
  # betas <- c()
  # loglik <- c(eval_log_lik(sigma2_Y, gwas_n, a, b, beta, ld_eig_val, tD, tbeta, tb))
  while(diff > tol & niter < max_iter) {
    a0 <- a
    b0 <- b
    beta0 <- beta
    # a <- sigma2_Y - b * gwas_n * sum(tbeta / tD * tbeta) - 2 * sum(tbeta / tD * tb) * beta + sum(tb * ld_eig_val / tD * tb) * beta ^ 2
    a <- eval_rt_Ainv_r_div_n(sigma2_Y, gwas_n, b, beta, ld_eig_val, tD, tbeta, tb)
    b <- solve_b(gwas_n, ld_eig_val, tbeta, tb, beta, a, b0)
    if(is.na(b)) {
      return(list(a = a, b = b, beta = beta, z = NA, conv = FALSE))
    }
    tD <- calc_tD(b, gwas_n, ld_eig_val)
    beta <- sum(tbeta / tD * tb) / sum(tb / tD * ld_eig_val * tb)
    diff <- sqrt(sum((a - a0) ^ 2) + sum((b - b0) ^ 2) + sum((beta - beta0) ^ 2))
    # diffs <- c(diffs, diff)
    # as <- c(as, a)
    # bs <- c(bs, b)
    # betas <- c(betas, beta)
    niter <- niter + 1
    # loglik <- c(loglik, eval_log_lik(sigma2_Y, gwas_n, a, b, beta, ld_eig_val, tD, tbeta, tb))
    
  }
  var_beta <- a / sum(tb * ld_eig_val / tD * tb) / gwas_n 
  return(list(a = a, b = b, beta = beta, z = beta / sqrt(var_beta), conv = diff < tol))
}


# helper functions
eval_rt_Ainv_r_div_n <- function(sigma2_Y, n, b, beta, D, tD, tbeta, tb) {
  sigma2_Y - b * n * sum(tbeta / tD * tbeta) - 2 * sum(tbeta / tD * tb) * beta + sum(tb * D / tD * tb) * beta ^ 2
}

get_ld_diag <- function(ld_eig_vec, ld_eig_val) {
  return(diag(ld_eig_vec %*% diag(ld_eig_val) %*% t(ld_eig_vec)))
}

calc_tD <- function(b, n, D) {
  return(1 + b * n * D)
}

func_dev_b <- function(b, n, D, tbeta, tb, beta, a) {
  tD <- calc_tD(b, n, D)
  return(
    - sum(D / tD) + n / a * (sum(tbeta / (tD ^ 2) * tbeta)) - 
      2 * sum(tb / (tD ^ 2) * D * tbeta) * beta + 
      sum(tb / (tD ^ 2) * (D ^ 2) * tb) * beta ^ 2)
}

func_obj_b <- function(b, n, D, tbeta, tb, beta, a, sigma2_Y) {
  tD <- calc_tD(b, n, D)
  return(
    -1 / 2 * sum(log(tD)) - 
      n / 2 / a * (- n * b * sum(tbeta / tD * tbeta) - 2 * sum(tbeta / tD * tb) * beta + sum(tb * D / tD * tb) * beta ^ 2))
}

solve_b <- function(n, D, tbeta, tb, beta, a, b0, func = func_obj_b, func_dev = func_dev_b) {
  # tmp0 <- func_dev_b(0, n, D, tbeta, tb, beta, a)
  # up_bound <- b0
  # tmp1 <- func_dev_b(up_bound, n, D, tbeta, tb, beta, a)
  # niter <- 0
  # while(tmp0 * tmp1 > 0 & niter < 100) {
  #   up_bound <- up_bound * 2
  #   tmp1 <- func_dev_b(up_bound, n, D, tbeta, tb, beta, a)
  #   niter <- niter + 1
  # }
  # h <- up_bound
  # l <- 0
  # if(tmp1 * tmp0 < 0) {
  #   res <- uniroot(func_dev_b, c(l, h), n = n, D = D, tbeta = tbeta, tb = tb, beta = beta, a = a)
  #   return(res$root)
  # } else {
  #   return(NA) 
  # }
  return(
    optim(
      b0, func, func_dev, 
      n = n, D = D, tbeta = tbeta, tb = tb, 
      beta = beta, a = a, sigma2_Y = sigma2_Y, 
      method = 'Brent', lower = 0, upper = 1, 
      control = list(fnscale = -1))$par)
}