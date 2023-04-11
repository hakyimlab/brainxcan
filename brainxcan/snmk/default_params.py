from collections import OrderedDict 

CHRS = [ i for i in range(1, 23) ]
PYTHON_EXE = 'python'
RSCRIPT_EXE = 'Rscript'
PLINK_EXE = 'plink'
LD_CLUMP_YAML = '{datadir}/mr/ld_clump.yaml'
BXCAN_VIS_DATADIR = '{datadir}/bxcan_vis'
BXCAN_IDP_META = '{datadir}/bxcan_vis/idp_meta_data.csv'
BXCAN_COLOR_CODE = '{datadir}/bxcan_vis/report_color_code.yaml'
SPEARMAN_CUTOFF = 0.1
BXCAN_REGION_VIS = False
BXCAN_EMPIRICAL_Z = False
BXCAN_EMPIRICAL_Z_SEED = 1
BXCAN_EMPIRICAL_Z_NREPEAT = 1000
BXCAN_LDBLOCK_PERM = None
BXCAN_LDBLOCK_PERM_SEED = 1
BXCAN_LDBLOCK_PERM_NREPEAT = 10
BXCAN_LDBLOCK_PERM_ELASTIC_NET = False


IDP_TYPE = [
    'residual',  # default
    ['original', 'residual']  # options 
]
MODEL_TYPE = [ 
    'ridge', # default
    ['ridge', 'elastic_net'] # options
]
GWAS_POP = 'EUR'
GWAS_POPS = ['EUR', 'SAS', 'AMR', 'EAS', 'AFR']

GENO_COV_PATTERN = '{datadir}/geno_covar/chr{chr_num}.banded.npz'
IDP_WEIGHTS_PATTERN = '{datadir}/idp_weights/{model_type}/{idp_type}.{idp_modality}.parquet'
IDP_GWAS_PATTERN = '{datadir}/idp_gwas/{idp_type}.{idp_modality}.chr{chr_num}/{idp_code}.parquet'
IDP_GWAS_SNP_PATTERN = '{datadir}/idp_gwas/snp_bim/chr{chr_num}.bim'
MR_LD_PANEL_PATTERN = '{datadir}/mr/ieugwasr/{gwas_pop}'
CORRECTION_FACTOR_EMP = 1
CORRECTION_FACTOR_PERM = 1.1

BXCAN_SIGNIF = OrderedDict([
    ('signif_pval', 5e-2),
    ('signif_max_idps', 10)
])

IDP_WEIGHTS_COLS = OrderedDict([
    ('snpid', 'variant_id'),
    ('effect_allele', 'a1'),
    ('non_effect_allele', 'a0'),
    ('chr', 'chr')
])

