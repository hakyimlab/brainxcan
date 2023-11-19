library(optparse)

option_list <- list(
    make_option(c("-i", "--mr_rds"), type="character", default=NULL,
                help="Input MR RDS (output of run_mr_local.R)",
                metavar="character"),
    make_option(c("-b", "--bxcan"), type="character", default=NULL,
                help="BrainXcan result (single line).",
                metavar="character"),
    make_option(c("-p", "--output_plot"), type="character", default=NULL,
                help="Output MR plot.",
                metavar="character"),
    make_option(c("-t", "--output_table"), type="character", default=NULL,
                help="Output MR table.",
                metavar="character"),
    make_option(c("-k", "--bxcan_pval_col"), type="character", default=NULL,
                help="The name of the p-value column to be used in BrainXcan results.",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

srcpath = Sys.getenv("R_VIS_SRC")
source(paste0(srcpath, '/', 'vis_helper.R'))

library(dplyr)
library(ggplot2)
library(patchwork)
options(stringsAsFactors = F)
theme_set(theme_classic(base_size = 12))

# logging config
logging::basicConfig(level = opt$log_level)

logging::loginfo('Loading MR RDS and BrainXcan 1-line result.')
mr_res = readRDS(opt$mr_rds)
bxcan = read.csv(opt$bxcan)

logging::loginfo('Plotting.')
p1 = mr_res$idp2pheno$data %>% ggplot() + geom_hline(yintercept = 0, color = 'gray') + 
  geom_vline(xintercept = 0, color = 'gray') + 
  geom_point(aes(x = beta.exposure, y = beta.outcome), alpha = 0.5) + 
  geom_errorbar(aes(x = beta.exposure, ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.5) +
  geom_errorbarh(aes(y = beta.outcome, xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.5) +
  xlab('SNP estimated effect in IDP') + ylab('SNP estimated effect in phenotype') + ggtitle('MR: IDP -> Phenotype')
p2 = mr_res$pheno2idp$data %>% ggplot() + geom_hline(yintercept = 0, color = 'gray') + 
  geom_vline(xintercept = 0, color = 'gray') + 
  geom_point(aes(x = beta.exposure, y = beta.outcome), alpha = 0.5) + 
  geom_errorbar(aes(x = beta.exposure, ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.5) +
  geom_errorbarh(aes(y = beta.outcome, xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.5) +
  xlab('SNP estimated effect in phenotype') + ylab('SNP estimated effect in IDP') + ggtitle('MR: Phenotype -> IDP')
ggsave(opt$output_plot, p1 + p2, height = 4, width = 8)

logging::loginfo('Collecting summary of results.')
mr_methods = c('Inverse variance weighted', 'Weighted median', 'MR Egger')
df_mr = rbind(
  mr_res$idp2pheno$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'IDP -> Phenotype'),
  mr_res$pheno2idp$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'Phenotype -> IDP') 
) %>% select(direction, method, nsnp, b, pval)
bxcan = bxcan %>% mutate(direction = NA, method = 'BrainXcan') %>% rename(nsnp = nsnp_used, b = bhat)
bxcan = bxcan[, c("direction", "method", "nsnp", "b", opt$bxcan_pval_col)]
colnames(bxcan)[5] = 'pval'
df_res = rbind(
  df_mr, 
  bxcan
)

# add ACAT based meta-analysis to combine MR tests
bxcan_direction = sign(df_res %>% filter(method == 'BrainXcan') %>% pull(b))
acat_res = list()
for(dd in c('IDP -> Phenotype', 'Phenotype -> IDP')) {
  sub = df_mr %>% filter(!is.na(direction) & direction == dd)
  p_acat = acat_signed(sub$pval, sub$b, bxcan_direction)
  acat_res[[length(acat_res) + 1]] = data.frame(direction = dd, method = 'ACAT meta-analysis', nsnp = NA, b = bxcan_direction, pval = p_acat)
}
acat_res = do.call(rbind, acat_res)
df_res = rbind(df_res, acat_res)

# add heterogeneity test
df_het = list()
if(nrow(mr_res$idp2pheno$heterogeneity) > 0) {
  df_het[[length(df_het) + 1]] <- mr_res$idp2pheno$heterogeneity %>% filter(method %in% mr_methods) %>% mutate(direction = 'IDP -> Phenotype')
}
if(nrow(mr_res$pheno2idp$heterogeneity) > 0) {
  df_het[[length(df_het) + 1]] <- mr_res$pheno2idp$heterogeneity %>% filter(method %in% mr_methods) %>% mutate(direction = 'Phenotype -> IDP')
}
if(length(df_het) > 0) {
  df_het = do.call(rbind, df_het) %>% select(direction, method, Q, Q_df, Q_pval)
  df_res = left_join(df_res, df_het, by = c('method', 'direction'))
}

write.table(df_res, opt$output_table, sep = '\t', col = T, row = F, quo = F)

logging::loginfo('Done.')

