load_color_code_yaml = function(fn) {
  kk = yaml::read_yaml(fn)
  res = c()
  for(n in names(kk)) {
    if(length(kk[n]) > 1) {
      message('More than one color code, will load the first one.')
    }
    res = c(res, kk[[n]][1])
  }
  names(res) = names(kk)
  res
}

gen_shape_map = function() {
  data.frame(
    lr = c('left', 'right', 'NA', 'vermis'),
    # code = c('\u25C0', '\u25B6', '\u25A0', '\u25A0'),
    code = c(60, 62, 20, 20),
    stringsAsFactors = F
  )
}

plot_bxcan_ordered <- function(xcandf, color_map, z_thres)
{
  shape_map = gen_shape_map()
  tmp = xcandf %>% mutate(lr = as.character(side))
  tmp$lr[is.na(tmp$lr)] = 'NA'
  tmp = tmp %>%
    left_join(shape_map, by = 'lr') %>% select(-lr) %>%
    mutate(kk = reorder(region, zscore ^ 2, FUN = max))
  tmp2 = tmp %>%
    group_by(kk) %>% summarize(zmax = max(zscore), zmin = min(zscore)) %>% ungroup()
  pp = tmp %>%  
    ggplot() + 
    geom_hline(yintercept = c(-z_thres, z_thres), col = 'gray') + 
    geom_segment(data = tmp2, aes(x = kk, xend = kk, y = zmin, yend = zmax), color = 'black', size = 0.1) + 
    # geom_text(aes(x = kk, y = zscore, color = subtype, label = code), family = "Arial Unicode MS") + 
    geom_point(aes(x = kk, y = zscore, color = subtype, shape = code), size = 4, alpha = 0.7) +
    scale_shape_identity() + 
    coord_flip() +
    scale_color_manual(values = color_map) + 
    ylab('BrainXcan z-score') + 
    theme(axis.title.y = element_blank()) + 
    theme(legend.title = element_blank())
  pp
}

plot_bxcan_ordered_heatmap <- function(xcandf0, color_map, z_thres) {
  shape_map = gen_shape_map()
  tmp = xcandf0 %>% mutate(lr = as.character(side))
  tmp$lr[is.na(tmp$lr)] = 'NA'
  lr2 <- paste0('-', tmp$lr)
  lr2[lr2 == '-NA'] <- ''
  region2 <- tmp$region
  region2[substr(region2, 1, 2) == 'PC'] <- 'PC'
  region2 <- paste0(region2, lr2)
  subtype2 <- tmp$subtype
  pcs <- subtype2[substr(subtype2, 1, 2) == 'PC'] 
  pcs <- unlist(lapply(strsplit(pcs, '-'), function(x) x[2]))
  subtype2[substr(subtype2, 1, 2) == 'PC'] <- pcs
  tmp = tmp %>%
    mutate(region2 = region2, subtype2 = subtype2) %>%
    mutate(kk = reorder(region2, zscore ^ 2, FUN = max))
  tmp2 <- tmp %>% filter(abs(zscore) > z_thres)
  pp = tmp %>%  
    ggplot() + 
    geom_tile(aes(x = kk, y = subtype2, fill = zscore)) +
    coord_flip() +
    ylab('BrainXcan z-score') + 
    theme(axis.title.y = element_blank()) + 
    theme(legend.title = element_blank()) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0, na.value = 'gray')
  if (nrow(tmp2) > 0) {
    pp = pp + 
      geom_text(data = tmp2, aes(x = kk, y = subtype2, label = "*"), color = 'white', size = 4)
  }
  pp
}

load_mr_sum_stats = function(prefix) {
  kk = Sys.glob(paste0(prefix, '*.MR_sumstats.tsv'))
  if(length(kk) == 0) {
    return(NULL)
  }
  df_mr = list()
  for(k in kk) {
    mm = stringr::str_match(basename(k), '(t1_|dmri_)([A-Z-0-9a-z]+).MR_sumstats.tsv')[, 3]
    tmp = read.delim2(k)
    tmp$IDP = mm
    df_mr[[length(df_mr) + 1]] = tmp
  }
  df_mr = do.call(rbind, df_mr)
  df_mr$b = as.numeric(df_mr$b)
  df_mr$pval = as.numeric(df_mr$pval)
  df_mr
}

gen_mr_sign = function(direction, pval) {
  df_code = data.frame(
    direction = c(1, -1),
    sign = c('\u25B2', '\u25BC')
  )
  kk = left_join(
    data.frame(direction = sign(direction), pval = pval),
    df_code, by = 'direction'
  )
  paste0(signif(kk$pval, digits = 3), ' (', kk$sign, ')')
}
