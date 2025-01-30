import re
import pandas as pd
import numpy as np
import pathlib
from collections import OrderedDict 
from brainxcan.sbxcan.util.misc import read_table

BASE_PAIR = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}
GC_NUMBER = 0.456
N_IDP_META_COLS = 4
N_GWAS_META_COLS = 5

def check_flip(a1, a2, b1, b2):
    res = []
    for _a1, _a2, _b1, _b2 in zip(a1, a2, b1, b2):
        res.append(_check_flip(_a1, _a2, _b1, _b2))
    return np.array(res)

def _check_flip(a0, a1, b0, b1):
    '''
    check if (a0, a1) and (b0, b1) are of the same direction.
    If there is nan or they don't match at all or ambiguious return nan
    Else if they are in the same direction, return 1
    Else return -1
    '''
    if a0 is np.nan or a1 is np.nan or b0 is np.nan or b1 is np.nan:
        return np.nan
    # remove ambiguious first.
    if a0 == BASE_PAIR[a1] or b0 == BASE_PAIR[b1]:
        return np.nan
    # exact match
    if a0 == b0 and a1 == b1:
        return 1
    # flip
    if a0 == b1 and a1 == b0:
        return -1    
    # compliment match
    if a0 == BASE_PAIR[b0] and a1 == BASE_PAIR[b1]:
        return 1
    # compliment flip
    if a0 == BASE_PAIR[b1] and a1 == BASE_PAIR[b0]:
        return -1  
    # if all above does not return, it has to be invalid.
    return np.nan

def rearrage_df_by_target(df, target, df_value_cols):
    df_res = target[['snpid', 'chr', 'effect_allele', 'non_effect_allele']]
    df_res = pd.merge(
        df_res, df, 
        on=['snpid', 'chr'], 
        suffixes=['_res', '_df'],
        how='left'
    )
    flip_factor = check_flip(
        a1=df_res.effect_allele_res, 
        a2=df_res.non_effect_allele_res,
        b1=df_res.effect_allele_df, 
        b2=df_res.non_effect_allele_df
    )
    # we need to carry the missingness when we move on
    with np.errstate(invalid='ignore'):
        df_res[df_value_cols] = df_res[df_value_cols] * flip_factor[:, np.newaxis]
    df_res.drop(
        columns=['effect_allele_df', 'non_effect_allele_df'], inplace=True
    )
    df_res.rename(
        columns={
            'effect_allele_res': 'effect_allele',
            'non_effect_allele_res': 'non_effect_allele'
        }, 
        inplace=True
    )
    return df_res

def harmonize_gwas_and_weight(gwas, weight):
    '''
    Harmonize GWAS to weight SNP set.
    But only keep the ones that present in both.
    '''
    df_common = pd.merge(
        gwas[['snpid', 'chr', 'position', 'effect_allele', 'non_effect_allele']],
        weight[['snpid', 'chr', 'effect_allele', 'non_effect_allele']],
        on=['snpid', 'chr'],
        suffixes=['_gwas', '_weight']
    )
    flip_factor = check_flip(
        a1=df_common.effect_allele_gwas, 
        a2=df_common.non_effect_allele_gwas,
        b1=df_common.effect_allele_weight, 
        b2=df_common.non_effect_allele_weight
    )
    
    # need to remove the invalid variant before moving on
    to_keep_ind = np.logical_not(np.isnan(flip_factor))
    df_common = df_common[ to_keep_ind ].reset_index(drop=True)
    flip_factor = flip_factor[ to_keep_ind ]
    
    df_common.drop(columns=['effect_allele_gwas', 'non_effect_allele_gwas'], inplace=True)
    df_common.rename(columns={'effect_allele_weight': 'effect_allele', 'non_effect_allele_weight': 'non_effect_allele'}, inplace=True)

    df_gwas = pd.merge(
        df_common[['snpid', 'chr', 'position', 'effect_allele', 'non_effect_allele']], 
        gwas.drop(columns=['effect_allele', 'non_effect_allele']), 
        on=['snpid', 'chr']
    )
    df_gwas.effect_size = df_gwas.effect_size * flip_factor
    
    df_weight = pd.merge(
        df_common[['snpid', 'chr', 'position', 'effect_allele', 'non_effect_allele']], 
        weight.drop(columns=['effect_allele', 'non_effect_allele']), 
        on=['snpid', 'chr']
    )
    return df_gwas, df_weight

def _parse_args(args_list, desired_cols=None, no_raise=False):
    fn = args_list[0]
    if not pathlib.Path(fn).is_file():
        raise ValueError('Filename is wrong. Cannot find the file.')
    dict = {}
    snpid_name = None
    desired_cols_tmp = []
    for i in args_list[1:]:
        tmp = i.split(':')
        if len(tmp) != 2:
            raise ValueError('Wrong args list. Need [col]:[name] pairs.')
        col, name = tmp
        if desired_cols is None:
            desired_cols_tmp.append(col) 
        elif col not in desired_cols:
            if no_raise is True:
                continue
            else:
                raise ValueError(f'Wrong col = {col}.')
        dict[col] = name
    rename_dict = OrderedDict()
    if desired_cols is None:
        desired_cols = desired_cols_tmp
    for dd in desired_cols:
        if dd not in dict:
            if no_raise is True:
                continue
            else:
                raise ValueError(f'Need to have col = {dd}.')
        rename_dict[dict[dd]] = dd
    return fn, rename_dict
    
def _parse_gwas_args(args_list, mode='effect_size'):
    if mode == 'effect_size':
        have_effect_size = True
    elif mode == 'zscore':
        have_effect_size = False
    else:
        raise ValueError(f'Wrong loading mode for GWAS file: mode = {mode}')
    # for kk in args_list:
    #     if 'effect_size:' in kk:
    #         have_effect_size = True
    if have_effect_size is True:
        desired_cols = [
            'snpid', 'non_effect_allele', 'effect_allele', 
            'effect_size', 'effect_size_se', 'chr', 'position'
        ]
    else:
        desired_cols = [
            'snpid', 'non_effect_allele', 'effect_allele', 
            'zscore', 'allele_frequency', 'sample_size', 'chr', 'position'
        ]
    fn, rename_dict = _parse_args(args_list, desired_cols, no_raise=True)
    for k, v in rename_dict.items():
        if v == 'snpid':
            snpid_name = k
            break
    return fn, rename_dict, snpid_name

def get_snpid_col(gwas_args_list):
    for i in gwas_args_list:
        if 'snpid:' in i:
            _, tmp = i.split(':')
            return tmp

def impute_b_from_z(zscore, af, n):
    if (n < 0).sum() > 0:
        print('Warning: There are sample size < 0. These variants will be discarded')
    if (af >= 1).sum() or (af <= 0).sum():
        print('Warning: There are af outside (0, 1). These variants will be discarded')
    se = 1 / np.sqrt(2 * n * af * (1 - af))
    bhat = zscore * se
    return bhat, se

def clean_up_chr(ll):
    for i in range(len(ll)):
        ll[i] = re.sub('chr', '', ll[i])
    return ll

def get_key_by_val(val, dict_):
    for i in dict_.keys():
        if dict_[i] == val:
            return i
    return None

def load_gwas(gwas_args_list, fill_pos_col=False):
    snpid_col = get_snpid_col(gwas_args_list[1:])
    # fn = gwas_args_list[0]
    fn, rename_dict = _parse_args(gwas_args_list, desired_cols=None)
    df = read_table(fn, indiv_col=snpid_col)
    if fill_pos_col:
        df['position'] = np.nan
        gwas_args_list.append('position:position')
    k_effect_size = get_key_by_val('effect_size', rename_dict)
    k_zscore = get_key_by_val('zscore', rename_dict)
    if k_effect_size is not None and k_effect_size in df.columns:
        _, rename_dict, snpid_col = _parse_gwas_args(gwas_args_list, mode='effect_size')
    elif k_zscore is not None and k_zscore in df.columns:
        _, rename_dict, snpid_col = _parse_gwas_args(gwas_args_list, mode='zscore')
    else:
        raise ValueError('We need either effect_size or zscore in GWAS file.')
    df.rename(columns={'indiv': snpid_col}, inplace=True)
    df.rename(columns=rename_dict, inplace=True)
    df.drop_duplicates('snpid', inplace=True)
    df.chr = clean_up_chr(list(df.chr.astype(str)))
    if 'effect_size' not in rename_dict.values():
        df['effect_size'], df['effect_size_se'] = impute_b_from_z(df.zscore, df.allele_frequency, df.sample_size)
    # some qc on gwas
    # remove se with 0 or inf
    # remove effect size with na
    df.effect_size_se.replace([0, np.inf, -np.inf], np.nan, inplace=True)
    df = df[ (~ df.effect_size.isna()) & (~ df.effect_size_se.isna()) ].reset_index(drop=True)
    
    desired_cols = [
        'snpid', 'non_effect_allele', 'effect_allele', 
        'effect_size', 'effect_size_se', 'chr', 'position'
    ]
    return df[desired_cols]

def _parse_idp_args(args_list):
    desired_cols = [
        'snpid', 'non_effect_allele', 'effect_allele', 'chr'
    ]
    fn, rename_dict = _parse_args(args_list, desired_cols)
    return fn, rename_dict
        
def load_idp(args_list, spearman_cutoff=0.1):
    fn, rename_dict = _parse_idp_args(args_list)
    df = pd.read_parquet(fn)
    df.rename(columns=rename_dict, inplace=True)
    df.chr = df.chr.astype(str)
    # load perf.tsv.gz
    fn2 = re.sub('.parquet', '.perf.tsv.gz', fn)
    df_perf = pd.read_csv(fn2, compression='gzip', sep='\t')
    df_perf.rename(columns={'phenotype': 'IDP'}, inplace=True)
    perf_cols = ['R2', 'Pearson', 'Spearman']
    df_perf.rename(columns={ k : f'CV_{k}' for k in perf_cols }, inplace=True)
    df_perf = df_perf[ df_perf.CV_Spearman >= spearman_cutoff ].reset_index(drop=True)
    cols_to_keep = list(df.columns[:N_IDP_META_COLS]) + list(df_perf.IDP)
    df = df[ cols_to_keep ].copy()
    df = df[ df.iloc[:, N_IDP_META_COLS:].values.sum(axis=1) != 0 ].reset_index(drop=True)
    return df, df_perf

def load_cov_meta(fn):
    fn = '.'.join(fn.split('.')[:-1])
    fn = fn + '.snp_meta.parquet'
    return pd.read_parquet(fn)

def genomic_control(zscore):
    chisq = np.power(zscore, 2)
    lambda_gc = np.median(chisq) / GC_NUMBER
    chisq_adj = chisq / lambda_gc
    z_adj = np.sqrt(chisq_adj) * np.sign(zscore)
    return z_adj, lambda_gc    
    
def variance_control(zscore, phi, gwas_n, h2):
    varz = 1 + phi * gwas_n * h2 
    z_adj = zscore / np.sqrt(varz)
    return z_adj, varz

def _cleanup_ldblock(df):
    df = df.sort_values(by=['chr', 'start', 'end'])
    df = df[ df.chr.isin([ str(i) for i in range(1, 23) ])]
    df = df.reset_index(drop=True)
    for i in range(1, 23):
        kk = df[df.chr == str(i)].reset_index(drop=True)
        for j in range(kk.shape[0] - 1):
            if kk.end[j] != kk.start[j + 1]:
                raise ValueError(f'''
                    LD block file has non-concatenating blocks at 
                    chr = {i}, end({j}) = {kk.end[j]}, start({j + 1}) = {kk.start[j + 1]}
                ''')
    return df

def load_ldblock(fn):
    df = pd.read_csv(fn, sep='\s+')
    df.chr = [ re.sub('^chr', '', i) for i in df.chr ]
    df.rename(columns={'stop': 'end'}, inplace=True)
    # BED file use base0
    df.start = df.start.astype(int) + 1
    df.end = df.end.astype(int) + 1
    df = _cleanup_ldblock(df)
    return df

def get_idxs_by_block(meta, block):
    res = []
    if len(list(meta.chr.unique())) > 1:
        raise ValueError('There are more than one chromosome in meta. Exit!')
    chrm = meta.chr[0]
    # for chrm in range(1, 23):
    block_sub = block[block.chr == str(chrm)].reset_index(drop=True)
    if block_sub.shape[0] == 0:
        raise ValueError('No desired chromosome in LD block')
    meta_w_idx = pd.DataFrame({'chr': meta.chr, 'pos': meta.position, 'idx': [ i for i in range(meta.shape[0])]})
    chrm = meta_w_idx.chr[0]
    meta_i = meta_w_idx[ meta_w_idx.pos < block_sub.start.values[0] ]
    if meta_i.shape[0] > 0:
        res.append(list(meta_i.idx))
    for i in range(block_sub.shape[0]):
        meta_i = meta_w_idx[ (meta_w_idx.pos < block_sub.end.values[i]) & (meta_w_idx.pos >= block_sub.start.values[i]) ]
        if meta_i.shape[0] > 0:
            res.append(list(meta_i.idx))
    meta_i = meta_w_idx[ meta_w_idx.pos >= block_sub.end.values[-1] ]
    if meta_i.shape[0] > 0:
        res.append(list(meta_i.idx))
    return res

def simulate_weights(weight, nrepeat):
    null_weight = np.random.normal(size=(weight.shape[0], nrepeat))
    null_weight[np.isnan(weight).sum(axis=1) != 0, :] = np.nan
    return null_weight

def permute_weights(weight, weight_meta, ld_blocks, nrepeat):
    n = weight.shape[1]
    not_nan = np.isnan(weight).sum(axis=1) == 0
    perm_weight = np.zeros((weight.shape[0], nrepeat * n))
    idxs_by_block = get_idxs_by_block(weight_meta.iloc[not_nan, :].reset_index(drop=True), ld_blocks)
    for i in range(nrepeat):
        idxs = []
        block_idxs = np.random.permutation(len(idxs_by_block))
        for j in block_idxs:
            idxs += idxs_by_block[j]
        perm_weight[not_nan, i * n : (i + 1) * n] = weight[idxs, :]
    perm_weight[~not_nan, :] = np.nan
    return perm_weight
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_sbrainxcan.py', description='''
        Run S-BrainXcan with pre-computed genotype covariance.
    ''')
    parser.add_argument('--genotype_covariance', help='''
        The genotype covariance computed in build_genotype_covariance.py
        Accept wildcard {chr_num}.
        Will automatically search for the corresponding meta SNP file.
    ''')
    parser.add_argument('--gwas', nargs='+', help='''
        Need to have column names for: 
            snpid, non_effect_allele, effect_allele, 
            effect_size, effect_size_se, chr, position.
        If there is no effect_size avaliable, it could 
        impute effect_size from zscore, allele_frequency, 
        sample_size.
        The format is: snpid:rsid_col, ..., chr:chr, position:position.
        Note that position column is NOT required if --ldblock_perm is not 
        specified.
    ''')
    parser.add_argument('--idp_weight', nargs='+', help='''
        The IDP weight table is in parquet format.
        It contains columns:
            snpid, effect_allele, non_effect_allele, chr.
        Along with all other columns for the IDPs.
        Specify the column names, e.g.: snpid:rsID, ..., chr:chr
    ''')
    parser.add_argument('--spearman_cutoff', type=float, default=0.1, help='''
        The CV Spearman cutoff applied to models. 
    ''')
    parser.add_argument('--ldblock_perm', default=None, help='''
        If want to obtain adjusted BrainXcan by LD block-based permutation, 
        use this argument and specific the LD block BED file here.
        IMPORTANT NOTE: Currently, the LD block-based permutation only works 
        for dense models (e.g. ridge models). 
    ''')
    parser.add_argument('--ldblock_perm_seed', default=1, 
        type=int, help='''
        The random seed when obtaining the LD block-based permutation.
    ''')
    parser.add_argument('--ldblock_perm_nrepeat', default=10, 
        type=int, help='''
        The number of repeats when obtaining the LD block-based permutation
    ''')
    parser.add_argument('--empirical_null', action='store_true', help='''
        If specified, will report z-score adjusted by empirical null (based on 
        random IDP prediction models).  
        The default number of repeats is 1000 and the random seed is 1.
        Set these values in --empirical_null_nrepeat and --empirical_null_seed.
        IMPORTANT NOTE: Currently, the empirical null only works for dense 
        models (e.g. ridge models).
    ''')
    parser.add_argument('--empirical_null_nrepeat', default=1000, 
        type=int, help='''
        The number of repeats when obtaining the empirical null
    ''')
    parser.add_argument('--empirical_null_seed', default=1, 
        type=int, help='''
        The random seed when obtaining the empirical null.
        IMPORTANT NOTE: Random seed will be over-written by --ldblock_perm_seed
    ''')
    parser.add_argument('--output_prefix', help='''
        The output CSV file prefix.
        Will return marginal test result.
    ''')
    # parser.add_argument('--z_ld_weight', type=float, default=1e-4, help='''
    #     LD = (1 - z_ld_weight) * LD + z_ld_weight * (Z @ Z.T)
    #     to avoid mis-specified LD.
    # ''')
    args = parser.parse_args()
    
    from tqdm import tqdm
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    from brainxcan.sbxcan.util.CovConstructor import CovMatrix
    # from util.susie_wrapper import run_susie_wrapper
    from brainxcan.sbxcan.util.misc import z2p
    
    need_pos_col = False
    if args.ldblock_perm is not None:
        need_pos_col = True
    
    logging.info('Loading GWAS.')
    df_gwas = load_gwas(args.gwas, fill_pos_col=not need_pos_col)
    # df_gwas columns: 
    # snpid, non_effect_allele, effect_allele, 
    # effect_size, effect_size_se, chr
    logging.info('GWAS SNP = {}'.format(df_gwas.shape[0]))
    
    logging.info('Loading IDP weights.')
    df_weight, df_perf = load_idp(
        args.idp_weight, 
        spearman_cutoff=args.spearman_cutoff
    )
    idp_names = list(df_weight.columns[N_IDP_META_COLS:])
    nidp = len(idp_names)
    nsnp_total = (df_weight.iloc[:, N_IDP_META_COLS:].values != 0).sum(axis=0)
    logging.info('IDP SNP = {} and number of IDPs = {}'.format(df_weight.shape[0], nidp))
   
    nrepeat_null = 0
    if args.empirical_null is True:
        nrepeat_null = args.empirical_null_nrepeat
        logging.info(f'''
            Generating IDP weights for the empirical null:
            nrepeat = {nrepeat_null} and seed = {args.empirical_null_seed}.''')
        np.random.seed(args.empirical_null_seed)
        idp_names += [ f'null_{i}' for i in range(nrepeat_null) ]
    
    nrepeat_perm = 0
    if args.ldblock_perm is not None:
        nrepeat_perm = args.ldblock_perm_nrepeat
        logging.info(f'''
            Generating permuted weights by LD block:
            nrepeat = {nrepeat_perm} and seed = {args.ldblock_perm_seed}''')
        np.random.seed(args.ldblock_perm_seed)
        ld_blocks = load_ldblock(args.ldblock_perm)
        idp_names += [ f'perm{j}_idp{k}' for j in range(nrepeat_perm) for k in range(nidp) ]
    
    nperm = nrepeat_perm * nidp
    nnull = nrepeat_null
     
    logging.info('Harmonizing GWAS and IDP weights.')
    # harmonize GWAS and IDP weight table so that they have the same set of 
    # SNPs (including direction).
    df_gwas, df_weight = harmonize_gwas_and_weight(df_gwas, df_weight)
    logging.info('{} SNPs left after harmonizing GWAS and IDP weights.'.format(df_gwas.shape[0]))
    
    # please refer to https://github.com/hakyimlab/yanyu-notebook/blob/master/notes/date_112420.Rmd
    # for the details of the S-BrainXcan formula
    # to take the following procedure.
    # 0. subset IDP and GWAS SNPs.
    # 1. Per chromosome
    #   1.1 obtain D(chr), S_R(chr), and var_R(chr).
    #   1.2 compute numer_b(chr) = Gamma(chr).T @ (var_R(chr) * b_gwas(chr)) 
    #   1.3 compute numer_z(chr) = Gamma(chr).T @ (S_R(chr) * z_gwas(chr)) 
    # 2. compute marginal test.
    #   2.1 D = sum_chr D(chr), var_D = diag(D), S_D = sqrt(var_D)
    #   2.2 beta_brainxcan = ( sum_chr numer_b(chr) ) / var_D
    #   2.3 z_brainxcan =  ( sum_chr numer_z(chr) ) / S_D
    # 3. run susieR.
    #   3.1 Sigma = D / S_D[:, np.newaxis] / S_D[np.newaxis, :]
    
    ntotal = nidp + nperm + nnull
    D = np.zeros((ntotal, ntotal))
    numer_b = np.zeros((ntotal))
    numer_z = np.zeros((ntotal))
    nsnp_used = np.zeros((ntotal)).astype(int)
    for i in range(1, 23):
        
        df_gwas_sub = df_gwas[ df_gwas.chr == str(i) ].reset_index(drop=True)
        df_weight_sub = df_weight[ df_weight.chr == str(i) ].reset_index(drop=True)
        if df_gwas_sub.shape[0] == 0:
            continue
        
        logging.info(f'Chromosome {i}: Loading genotype covariance meta information.')
        df_cov_meta = load_cov_meta(args.genotype_covariance.format(chr_num=i))
        
        # step0
        n0 = df_weight_sub.shape[0]  # for book keeping
        # we enforce the GWAS table and the IDP weights to have 
        # the same SNPs as genotype covariance
        # the weights of the missing ones are set to NaN.
        df_gwas_sub = rearrage_df_by_target(
            df=df_gwas_sub, 
            target=df_cov_meta,
            df_value_cols=['effect_size']
        )
        df_weight_sub = rearrage_df_by_target(
            df=df_weight_sub, 
            target=df_cov_meta,
            df_value_cols=list(df_weight.columns[N_GWAS_META_COLS:])
        )
        n1 = df_gwas_sub.effect_size.notna().sum()
        logging.info('Step0 Chromosome {}: {} out of {} SNPs in IDP/GWAS are used.'.format(i, n1, n0))
        
        logging.info(f'Step1 Chromosome {i}: Working with genotype covariance.')
        
        weight_all = []
        weight = df_weight_sub.iloc[:, N_GWAS_META_COLS: ].to_numpy(copy=True)
        weight_all.append(weight)
        if args.empirical_null is True:
            weight_null = simulate_weights(weight, nrepeat_null)
            weight_all.append(weight_null)
        if args.ldblock_perm is not None:
            weight_perm = permute_weights(
                weight, df_weight_sub.iloc[:, :N_GWAS_META_COLS], 
                ld_blocks[ld_blocks.chr == str(i)].reset_index(drop=True), 
                nrepeat_perm)
            weight_all.append(weight_perm)
        weight = np.concatenate(weight_all, axis=1)
        weight[np.isnan(weight)] = 0
        nsnp_used += (weight != 0).sum(axis=0)
        
        b_gwas = df_gwas_sub.effect_size.to_numpy(copy=True)
        b_gwas[np.isnan(b_gwas)] = 0
        se_gwas = df_gwas_sub.effect_size_se.to_numpy(copy=True)
        se_gwas[np.isnan(se_gwas)] = 1
        z_gwas = b_gwas / se_gwas
        
        cov_mat = CovMatrix(args.genotype_covariance.format(chr_num=i))
        cov_x_weight, diag_cov = cov_mat.eval_matmul_on_left(weight, param=100)  # param is used in case of naive geno cov in h5 format
        D_chr = weight.T @ cov_x_weight
        del cov_x_weight
        var_R_chr = diag_cov
        numer_b_chr = weight.T @ ( var_R_chr * b_gwas )
        numer_z_chr = weight.T @ ( np.sqrt(var_R_chr) * z_gwas )
        D += D_chr
        numer_b += numer_b_chr
        numer_z += numer_z_chr
    
    # add a checker to handle traits with D == 0
    n_bad_pheno = np.sum(D.diagonal() == 0)
    if n_bad_pheno > 0:
        good_pheno = D.diagonal() != 0
        D = D[:, good_pheno][good_pheno, :]
        numer_b = numer_b[good_pheno]
        numer_z = numer_z[good_pheno]
        idp_names = list(np.array(idp_names)[good_pheno])
        nsnp_used = nsnp_used[good_pheno]
        nsnp_total = nsnp_total[good_pheno]
        
    logging.info('Step2: Computing marginal test.')
    S_D = np.sqrt(D.diagonal())
    beta_brainxcan = numer_b / np.power(S_D, 2)
    z_brainxcan = numer_z / S_D
    
    # logging.info('Step3: Running susieR.')
    # Sigma =  D / S_D[:, np.newaxis] / S_D[np.newaxis, ]
    # susie_pip, susie_cs = run_susie_wrapper(z_brainxcan, Sigma, params={'z_ld_weight': args.z_ld_weight})
         
    logging.info('Saving outputs.')
    df_res = pd.DataFrame({
        'IDP': idp_names[: nidp],
        'bhat': beta_brainxcan[: nidp],
        'pval': z2p(z_brainxcan[: nidp]),
        'z_brainxcan': z_brainxcan[: nidp],
        'nsnp_used': nsnp_used[: nidp],
        'nsnp_total': nsnp_total[: nidp]
        # 'pip': susie_pip,
        # 'cs95': susie_cs
    })
    
    if args.empirical_null is True:
        df_adj = pd.DataFrame({
            'name': idp_names[nidp : (nidp + nnull)], 
            'value': z_brainxcan[nidp : (nidp + nnull)] })
        df_adj['rand'] = args.empirical_null_seed
        df_adj = df_adj.to_csv(args.output_prefix + '.emp_null.csv', index=False)
    if args.ldblock_perm is not None:
        df_adj = pd.DataFrame({
            'name': idp_names[(nidp + nnull) : (nidp + nnull + nperm)], 
            'value': z_brainxcan[(nidp + nnull) : (nidp + nnull + nperm)] })
        df_adj['rand'] = args.ldblock_perm_seed
        df_adj = df_adj.to_csv(args.output_prefix + '.perm_null.csv', index=False)
            
    df_res = pd.merge(df_res, df_perf, on='IDP', how='left')
    df_res.fillna('NA', inplace=True)
    df_res.sort_values(by='pval').to_csv(args.output_prefix + '.csv', index=False)
    
    logging.info('Done.')
    
