# example run of mix-BrainXcan

# setwd('brainxcan/mix_brainxcan/')
source('mix_brainxcan_solver.R')
source('data_simulator.R')
source('gwas_and_s_brainxcan.R')

# simulation procedure
# null model:
# Y = X beta_Y + e_Y
# 
# alternative model:
# Y = (X beta_E) beta_mediate + X beta_Y + e_Y
# 
# calculating S-BrainXcan and mix-BrainXcan 
# weights = beta_E (i.e. assuming a perfect genetic predictor of E)
# GWAS on Y is obtained from Y and X


set.seed(2023)
h2_Y <- 0.1
n_snp <- 500
n_sample <- 1200
ld <- simulate_ld(n_snp, n_blk = 5)
X <- simulate_genotype(n_sample, ld)
cov_X <- cov(X)
ld_eig <- decomp_ld(cov_X)
beta_E <- rnorm(n_snp)
beta_Y <- rnorm(n_snp)
beta_mediate <- rnorm(1)
E_genetic <- X %*% beta_E
Y_genetic_null <- X %*% beta_Y
Y_genetic_alt <- E_genetic * beta_mediate + X %*% beta_Y
Y_null <- Y_genetic_null + rnorm(n_sample, sd = sqrt(get_error_var(Y_genetic_null, h2_Y)))
Y_alt <- Y_genetic_alt + rnorm(n_sample, sd = sqrt(get_error_var(Y_genetic_alt, h2_Y)))

gwas_Y_null <- run_gwas(Y_null, X)
gwas_Y_alt <- run_gwas(Y_alt, X)

results <- list()
results[['sbrainxcan_null']] <- run_s_brainxcan(
  betahat = gwas_Y_null[, 1], se_betahat = gwas_Y_null[, 3], 
  ld_eig_vec = ld_eig$ld_eig_vec, ld_eig_val = ld_eig$ld_eig_val,
  weights = beta_E)
results[['sbrainxcan_alt']] <- run_s_brainxcan(
  betahat = gwas_Y_alt[, 1], se_betahat = gwas_Y_alt[, 3], 
  ld_eig_vec = ld_eig$ld_eig_vec, ld_eig_val = ld_eig$ld_eig_val,
  weights = beta_E)
results[['mix_brainxcan_null']] <- solve_mixed_brainxcan(
  gwas_Y_null[, 1], gwas_Y_null[, 3], n_sample, 
  ld_eig$ld_eig_vec, ld_eig$ld_eig_val, beta_E)
results[['mix_brainxcan_alt']] <- solve_mixed_brainxcan(
  gwas_Y_alt[, 1], gwas_Y_alt[, 3], n_sample, 
  ld_eig$ld_eig_vec, ld_eig$ld_eig_val, beta_E)
cols <- c('beta', 'z')
do.call(rbind, lapply(results, function(x) x[cols] %>% as.data.frame))
message('true mediation effect size = ', beta_mediate)
