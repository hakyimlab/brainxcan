p2z = function(p, b) {
  sign(b) * abs(qnorm(p / 2, lower.tail = T))
}

gen_five_scales = function(sig_cutoff, nominal_cutoff) {
  res = c(
    paste0('z > ', sig_cutoff),
    paste0(nominal_cutoff, ' < z <= ', sig_cutoff),
    paste0(-nominal_cutoff, ' <= z <= ', nominal_cutoff),
    paste0(-sig_cutoff, ' <= z < ', -nominal_cutoff),
    paste0('z < ', - sig_cutoff)
  )
  return(res)
}

zscore2category = function(zval, sig_cutoff = 4, nominal_cutoff = 2, colors = c('#cb0f0f', '#ff844c', '#FFFFFF', '#7fd1ff', '#0043ad')) {
  sig_cutoff = abs(sig_cutoff)
  nominal_cutoff = abs(nominal_cutoff)
  scales = gen_five_scales(sig_cutoff, nominal_cutoff)
  oo = rep(scales[3], length(zval))
  oo[zval > sig_cutoff] = scales[1]
  oo[zval < - sig_cutoff] = scales[5]
  oo[zval <= sig_cutoff & zval > nominal_cutoff] = scales[2]
  oo[zval >= - sig_cutoff & zval < - nominal_cutoff] = scales[4]
  color_scale = colors
  names(color_scale) = scales
  oo[is.na(zval)] = NA
  oo = factor(oo, levels = scales)
  return(list(category = oo, color_code = color_scale))
}

prep_brain_img_and_region_sf = function(region_img, bg_img, vis_meta) {
  tt = reshape2::melt(region_img) %>% filter(value != 0)
  bb = reshape2::melt(bg_img) %>% filter(value != 0) %>% rename(bg = value)
  bb = left_join(bb, tt, by = c('Var1', 'Var2', 'Var3'))
  # mid1 = bg_img@dim_[2]
  # mid2 = bg_img@dim_[3]
  # mid3 = bg_img@dim_[4]
  df_slices = vis_meta$slide_position
  res = list()
  for(i in 1 : nrow(df_slices)) {
    dim_ = df_slices$dim[i]
    idx_ = df_slices$idx[i]
    if(dim_ == 1) {
      tmp = bb %>% filter(Var1 == idx_) %>% mutate(direction = paste0('slice', i)) %>% rename(x = Var2, y = Var3, ref = Var1)
    } else if(dim_ == 2) {
      tmp = bb %>% filter(Var2 == idx_) %>% mutate(direction = paste0('slice', i)) %>% rename(x = Var1, y = Var3, ref = Var2)
    } else if(dim_ == 3) {
      tmp = bb %>% filter(Var3 == idx_) %>% mutate(direction = paste0('slice', i)) %>% rename(x = Var1, y = Var2, ref = Var3)
    }
    res[[length(res) + 1]] = tmp
  }
  # d1 = vis_meta$slide_position[1]
  # d2 = vis_meta$slide_position[2]
  # d3 = vis_meta$slide_position[3]
  # tmp = rbind(
  #   bb %>% filter(Var1 == floor(d1)) %>% mutate(direction = 'a1') %>% rename(x = Var2, y = Var3, ref = Var1),
  #   bb %>% filter(Var2 == floor(d2)) %>% mutate(direction = 'a2') %>% rename(x = Var1, y = Var3, ref = Var2),
  #   bb %>% filter(Var3 == floor(d3)) %>% mutate(direction = 'a3') %>% rename(x = Var1, y = Var2, ref = Var3)
  # )
  tmp = do.call(rbind, res)
  
  outline = list()
  for(kk in unique(tmp$direction)) {
    outline[[length(outline) + 1]] = convert_to_sf(
      tmp %>% filter(direction == kk) %>% 
        filter(!is.na(value)) %>% select(x, y, value)
    ) %>% mutate(direction = kk)
  }
  outline = do.call(rbind, outline)
  list(bg = tmp, region_sf = outline)
}

vis_by_tag = function(datadir, tag, df, score) {
  vis_bg = readRDS(paste0(datadir, '/bg_img.rds'))
  vis_meta = readRDS(paste0(datadir, '/meta_plot.rds'))
  if(!tag %in% names(vis_meta)) {
    message('Wrong tag.')
    return(NULL)
  } else {
    vis_meta = vis_meta[[tag]]
  }
  df_sub = df %>% filter(IDP %in% vis_meta$IDP)
  vis = readRDS(paste0(datadir, '/', vis_meta$vis_rds))
  
  sub = vis$df_color %>% filter(IDP %in% df_sub$IDP)
  if(sum(duplicated(sub$color_code)) > 0) {
    message('Wrong IDP set.')
    return(NULL)
  }
  
  df_color = left_join(sub, df_sub, by = 'IDP')
  
  res = prep_brain_img_and_region_sf(vis$img, vis_bg, vis_meta)
  
  if(score == 'zscore') {
    kk = zscore2category(df_color[[score]])
    df_color$value_category = kk$category
    res$region_sf = res$region_sf %>% 
      left_join(df_color, by = c('value' = 'color_code'))
    p = res$bg %>% ggplot() +
      geom_raster(aes(x, y, fill = -bg)) +
      scale_fill_gradient(na.value = "transparent", low = 'white', high = 'gray20', guide = 'none') +
      ggnewscale::new_scale_fill() + 
      geom_sf(data = res$region_sf, aes(fill = value_category)) +
      scale_fill_manual(values = kk$color_code, na.value = 'transparent', na.translate = FALSE) +
      coord_sf() + 
      facet_wrap(.~direction, labeller = label_both, ncol = 3)
      
  } else {
    df_color$value = df_color[[score]]
    res$region_sf = res$region_sf %>% 
      left_join(df_color, by = c('value' = 'color_code'))
    p = res$bg %>% ggplot() + 
      geom_raster(aes(x, y, fill = -bg)) +
      scale_fill_gradient(na.value = "transparent", low = 'white', high = 'gray20', guide = 'none') +
      ggnewscale::new_scale_fill() + 
      geom_sf(data = res$region_sf, aes(fill = value_category)) +
      scale_fill_manual(values = kk$color_code, na.value = 'transparent', na.translate = FALSE) +
      coord_sf() + 
      facet_wrap(~direction, labeller = label_both, ncol = 3)
  }
  
  p = p + 
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )
  p = p + theme(strip.background = element_blank(), strip.text.x = element_blank())
  
  p
}

vis_region = function(datadir, tag, df) {
  vis_bg = readRDS(paste0(datadir, '/bg_img.rds'))
  vis_meta = readRDS(paste0(datadir, '/meta_plot.rds'))
  if(!tag %in% names(vis_meta)) {
    message('Wrong tag.')
    return(NULL)
  } else {
    vis_meta = vis_meta[[tag]]
  }
  df_sub = df %>% filter(IDP %in% vis_meta$IDP)
  vis = readRDS(paste0(datadir, '/', vis_meta$vis_rds))
  
  sub = vis$df_color %>% filter(IDP %in% df_sub$IDP)
  if(sum(duplicated(sub$color_code)) > 0) {
    message('Wrong IDP set.')
    return(NULL)
  }
  
  df_color = left_join(sub, df_sub, by = 'IDP')
  
  res = prep_brain_img_and_region_sf(vis$img, vis_bg, vis_meta)
  res$region_sf = res$region_sf %>% 
    left_join(df_color, by = c('value' = 'color_code'))
  
  p = res$bg %>% ggplot() +
    geom_raster(aes(x, y, fill = -bg)) +
    scale_fill_gradient(na.value = "transparent", low = 'white', high = 'gray20', guide = 'none') + 
    geom_sf(
      data = res$region_sf %>% filter(!is.na(region)), 
      aes(color = region)
    ) +
    coord_sf() + 
    facet_wrap(~ direction, labeller = label_both, ncol = 3)
  
  p = p + 
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )
  p = p + theme(strip.background = element_blank(), strip.text.x = element_blank())
  
  p
}

add_title = function(p, label, tt = 10, bb = -20, ...) {
  p = p + ggtitle(label) + 
  theme(
    plot.title = element_text(
      margin = margin(t = tt, b = bb), 
      ...
    )
  )
  p
  # df_title = data.frame(x = x, y = y, direction = 'a1')
  # geom_text(data = df_title, aes(x = x, y = y), label = label, hjust = 0, ...)
}

# convert heatmap data.frame to sf object
# require packages: raster, sf, stars
# modified from: https://stackoverflow.com/questions/34756755/plot-outline-around-raster-cells
convert_to_sf = function(dd) {
  ## Convert your data.frame to a raster object
  r <- raster::rasterFromXYZ(dd)
 
  ## some ugly check but have to ..
  if(names(r) != 'value') {
    names(r) = 'value'
  } 
  ## Extract polygons
  pp <- raster::rasterToPolygons(r, dissolve=TRUE)
  
  mm = stars::st_as_stars(pp)
  mm = sf::st_as_sf(mm) %>% sf::st_cast("MULTIPOLYGON")
  mm
}

acat_vanilla = function(p_vec) {
  T_acat = sum(tan((0.5 - p_vec) * pi))
  # T_acat
  1 / 2 - atan(T_acat / length(p_vec)) / pi
}

# DEPREACATED since we don't want to use p / 2
# acat_signed = function(p_vec, stat, reference_direction) {
#   sign_stat = sign(stat) * reference_direction
#   pp = p_vec / 2
#   pp[sign_stat < 0] = 1 - p_vec[sign_stat < 0] / 2
#   acat_vanilla(pp)
# }

acat_signed = function(p_vec, stat, reference_direction) {
  sign_stat = sign(stat) * reference_direction
  pp = p_vec 
  pp[sign_stat < 0] = 1 - p_vec[sign_stat < 0] 
  acat_vanilla(pp)
}

# DEPRECATED since we refactored to visualize region using geom_raster
vis_by_tag_DEPRECATED = function(datadir, tag, df, score) {
  vis_bg = readRDS(paste0(datadir, '/bg_img.rds'))
  vis_meta = readRDS(paste0(datadir, '/meta_plot.rds'))
  if(!tag %in% names(vis_meta)) {
    message('Wrong tag.')
    return(NULL)
  } else {
    vis_meta = vis_meta[[tag]]
  }
  df_sub = df %>% filter(IDP %in% vis_meta$IDP)
  vis = readRDS(paste0(datadir, '/', vis_meta$vis_rds))
  
  sub = vis$df_color %>% filter(IDP %in% df_sub$IDP)
  if(sum(duplicated(sub$color_code)) > 0) {
    message('Wrong IDP set.')
    return(NULL)
  }
  
  img = vis$img
  img2 = img
  img2[] = NA
  for(i in 1 : nrow(sub)) {
    img2[img == sub$color_code[i]] = df_sub[[score]][ df_sub$IDP == sub$IDP[i] ] 
  }
  tt = reshape2::melt(img2) %>% filter(value != 0)
  bb = reshape2::melt(vis_bg) %>% filter(value != 0) %>% rename(bg = value)
  bb = left_join(bb, tt, by = c('Var1', 'Var2', 'Var3'))
  # mid1 = vis_bg@dim_[2]
  # mid2 = vis_bg@dim_[3]
  # mid3 = vis_bg@dim_[4]
  d1 = vis_meta$slide_position[1]
  d2 = vis_meta$slide_position[2]
  d3 = vis_meta$slide_position[3]
  tmp = rbind(
    bb %>% filter(Var1 == floor(d1)) %>% mutate(direction = 'a1') %>% rename(x = Var2, y = Var3, ref = Var1),
    bb %>% filter(Var2 == floor(d2)) %>% mutate(direction = 'a2') %>% rename(x = Var1, y = Var3, ref = Var2),
    bb %>% filter(Var3 == floor(d3)) %>% mutate(direction = 'a3') %>% rename(x = Var1, y = Var2, ref = Var3)
  )
  
  if(score == 'zscore') {
    kk = zscore2category(tmp$value)
    tmp$value_category = kk$category
    p = tmp %>% ggplot() +
      scale_fill_manual(values = kk$color_code, na.value = 'transparent', na.translate = FALSE) + 
      # scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0, na.value = 'transparent') +
      geom_tile(aes(x, y, alpha = -bg), fill = "grey20") +
      scale_alpha(range = c(0.2, 0.8)) + 
      geom_raster(aes(x, y, fill = value_category)) +
      coord_equal() + facet_grid(.~direction, labeller = label_both) +
      guides(color = guide_legend("-bg"), alpha = FALSE)
  } else {
    p = tmp %>% ggplot() + 
      geom_raster(aes(x, y, fill = value)) +
      scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0, na.value = 'transparent') +
      geom_tile(aes(x, y, alpha = -bg), fill = "grey20") +
      scale_alpha(range = c(0.2, 0.8)) +
      coord_equal() + facet_grid(.~direction, labeller = label_both) +
      guides(color = guide_legend("-bg"), alpha = FALSE)
  }
  
  p + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.title = element_blank())
}

