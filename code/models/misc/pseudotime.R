create_princurve_trace <- function(df, stretch_factor){
  df_tmp_init <- df %>% 
    group_by(concentration) %>% 
    summarize(v1 = mean(v1),
              v2 = mean(v2)) %>% 
    ungroup()
  
  fit <- principal_curve(x = df %>% dplyr::select(v1, v2) %>% as.matrix(),
                         start = df_tmp_init %>% dplyr::select(v1, v2) %>% as.matrix(),
                         approx_points = FALSE, 
                         trace = TRUE, 
                         plot_iterations = FALSE, 
                         maxit = 100,
                         stretch = stretch_factor)
  
  plot(fit)
  points(df_tmp_init%>% dplyr::select(v1, v2) %>% as.matrix())
  
  result <- list()
  
  result$fit <- fit$s[fit$ord,] %>% as_tibble()
  result$center <- df_tmp_init
  result$input <- df
  
  return(result)
}


create_plot_df_bulk <- function(loi, drug_order, n_sample, v1_cut = 99999, stretch = 2){
  set.seed(123)
  # I am pulling in all drug observations and set concentration to 1
  drug_df <- umap_tidy %>%
    filter(partition %in% c(1,2)) %>%
    filter(drug %in% drug_order) %>%
    filter(line %in% loi) %>%
    mutate(concentration = ifelse(concentration == "nan", 1, concentration))
  # I identify DMSO treated observations and set a concentration to 0
  dmso_df <- umap_tidy %>%
    filter(partition %in% c(1,2)) %>%
    filter(drug == "DMSO") %>%
    filter(line %in% loi) %>%
    mutate(concentration = 0) 
  # I reassign every DMSO concentration to a drug of interest
  assigned = rep(drug_order, times = nrow(dmso_df)/length(drug_order))
  # I join back both tables having created artificial drug treatments with concetration 0 from DMSO data
  center_df <- dmso_df %>% 
    head(length(assigned)) %>%
    cbind(assigned) %>% 
    mutate(drug = assigned) %>% 
    dplyr::select(-assigned) %>% 
    rbind(drug_df) %>%
    group_by(drug, line, concentration) %>%
    filter(v1 <= v1_cut) %>%
    sample_n(n_sample, replace = TRUE) %>% ungroup() # adjust based on pan-line
  
  plot_df <- center_df %>%
    nest(-line, -drug) %>%
    mutate(plot = purrr::map(data, ~ create_princurve_trace(.x, stretch_factor = stretch)))
  
  return(plot_df)
}

draw_trajectory_bulk <- function(df){
  gg <- ggplot() + 
    geom_point_rast(aes(v1, v2), data = umap_tidy %>% sample_frac(0.01) %>% filter(partition %in% c(1,2)) %>% dplyr::select(-line, -drug), alpha = 1, size = 0.35, color = "#f1f1f1") + 
    geom_point(data = df %>% unnest(data),
               aes(color = drug, v1, v2), alpha = 0.5, size = 0.35) +
    geom_path(data = df %>% unnest(plot_trace), 
              aes(v1, v2, group = paste0(drug, line), color = drug), size = 1.5, 
              arrow = arrow(angle = 10, ends = "last", type = "closed", length = unit(0.15, "inches")))+ 
    scale_color_brewer(type = "qual", palette = "Set2") +
    theme_cowplot() +
    labs(x = "UMAP 1",
         y = "UMAP 2")+
    theme(legend.position = "bottom") + 
    facet_wrap(~ line)
  
  return(gg)
}


create_plot_df <- function(loi, drug_order, n_sample, v1_cut = 99999, stretch = 2){
  set.seed(123)
  
  center_df <- umap_tidy %>%
    filter(partition %in% c(1,2)) %>%
    filter(drug %in% drug_order) %>%
    filter(line %in% loi) %>% #, "D027T01", "D019T01"
    filter(concentration != "nan") %>%
    mutate(drug = factor(drug, levels =drug_order)) %>% 
    arrange((drug)) %>% 
    filter(v1 <= v1_cut) %>%
    group_by(drug, line, concentration) %>%
    sample_n(n_sample, replace = TRUE) %>% ungroup() # adjust based on pan-line vs. 
  
  plot_df <- center_df %>%
    nest(-line, -drug) %>%
    mutate(plot = purrr::map(data, ~ create_princurve_trace(.x, stretch_factor = stretch)))
  
  return(plot_df)
}

draw_trajectory <- function(df){
  gg <- ggplot() + 
    geom_point_rast(aes(v1, v2), data = umap_tidy %>% sample_frac(0.01) %>% filter(partition %in% c(1,2)) %>% dplyr::select(-drug), alpha = 1, size = 0.35, color = "#f1f1f1") + 
    geom_point(data = df %>% unnest(data),
               aes(color = line, v1, v2), alpha = 0.5, size = 0.35) +
    geom_path(data = df %>% unnest(plot_trace), 
              aes(v1, v2, group = paste0(drug, line), color = line), size = 1.5, 
              arrow = arrow(angle = 10, ends = "last", type = "closed", length = unit(0.15, "inches")))+ 
    scale_color_brewer(type = "qual", palette = "Set2") +
    theme_cowplot() +
    labs(x = "UMAP 1",
         y = "UMAP 2")+
    theme(legend.position = "bottom") + 
    facet_wrap(~ drug)
  
  return(gg)
}