## ---- results='hide', warning=F, message=F---------------------------------
library(perm)
library(limma)
library(pheatmap)
library(PharmacoGx)
library(reshape2)
library(patchwork)
library(ggraph)
library(tidygraph)
library(ggrepel)
library(biobroom)
library(tidyverse)

## ---- results='hide', warning=F, message=F---------------------------------
data('ctg_data_raw', package='SCOPEAnalysis')

## ---- results='hide', warning=F, message=F---------------------------------
## read plate annotation from R data file
data('plate_anno', package='SCOPEAnalysis')

## join with ctg data
ctg_data <- ctg_data_raw %>% inner_join(plate_anno) %>% 
  extract(well, c('row', 'col'), regex='([A-P])(\\d{2})', remove=F) %>% 
  mutate(row=match(row, LETTERS), col=as.integer(col))

## ---- results='hide', warning=F, message=F---------------------------------
ctg_data <- ctg_data %>% 
  mutate(pcount = ifelse(screen == '170502_NR_M2_D004T01P006L08' & 
                         (row %in% seq(1, 15, by=2)) & 
                         (col %in% 1:13), NA, pcount))

## ---- results='hide', warning=F, message=F---------------------------------
## colour palette for drawing
colpal <- colorRampPalette(rev(c('#b2182b','#d6604d','#f4a582',
                                  '#fddbc7','#f7f7f7','#d1e5f0',
                                  '#92c5de','#4393c3','#2166ac')))(150)

## box plot of DMSO photon counts
ctg_data %>% filter(drug == 'DMSO') %>% 
  ggplot(aes(rack.concentration, pcount)) + 
  geom_jitter(width=0.2) + 
  facet_wrap(~screen) + geom_boxplot(alpha=0) + 
  theme_bw() + theme(panel.grid=element_blank())

## heat map of DMSO photon counts, ignoring other values
ctg_data %>% mutate(pcount = ifelse(drug != 'DMSO', -1000, pcount)) %>% 
  ggplot(aes(row, col, fill = pcount)) + geom_raster() +
  facet_wrap(~screen) + theme_bw() +
  theme(panel.grid=element_blank()) +
  scale_fill_gradientn(colors=colpal)

## heat map of all photon counts
ctg_data %>% ggplot(aes(row, col, fill = pcount)) + 
  geom_raster() +
  facet_wrap(~screen) + theme_bw() +
  theme(panel.grid=element_blank()) +
  scale_fill_gradientn(colors=colpal)

## ---- results='hide', warning=F, message=F---------------------------------
## scatter plot with robust loess fit for rows.
ctg_data %>% filter(drug == 'DMSO') %>% 
  ggplot(aes(row, pcount)) + geom_point() + 
  geom_smooth() +
  facet_wrap(~screen) + theme_classic()

## scatter plot with robust loess fit for columns
ctg_data %>% filter(drug == 'DMSO') %>%
  ggplot(aes(col, pcount)) + geom_point() + 
  geom_smooth() +
  facet_wrap(~screen) + theme_classic()

## ---- results='hide', warning=F, message=F---------------------------------
## split data by plate and apply loess normalization
ctg_loess <- ctg_data %>% split(.$screen) %>% lapply(function(s){
  ## loess fit. include row and column effects, fit only on DMSO
  fit <- loess(pcount ~ 0 + row + col, data=filter(s, drug=='DMSO'),
               surface = 'direct')
  ## apply normalization
  s %>% mutate(norm_fac = predict(fit, data.frame(row=row, col=col)),
               pcount_norm = pcount - (norm_fac - median(norm_fac)))
}) %>% bind_rows()

## ---- results='hide', warning=F, message=F---------------------------------
ctg_loess <- ctg_loess %>% 
  mutate(pcount_norm = ifelse(pcount_norm < 0, 0, pcount_norm),
         screen = gsub('HCO_', 'HCO-', screen)) %>%
  separate(screen, c('date', 'tag1', 'tag2', 'plate'), sep='_', remove=F) %>% 
  mutate(date = as.Date(paste0('20', date), format='%Y%m%d'))

## ---- results='hide', warning=F, message=F---------------------------------
## median of dmso controls
ctrls_med <- ctg_loess %>% filter(drug == 'DMSO') %>%
  group_by(screen) %>% 
  summarise(med_ctrl = median(pcount_norm, na.rm=T)) %>% 
  ungroup() 

## normalize to median of DMSO
ctg_norm <- ctg_loess %>% inner_join(ctrls_med) %>%
  group_by(screen) %>% 
  mutate(viability = pcount_norm/med_ctrl) %>% ungroup() %>%
  extract(screen, 'line', regex = '.+(D0.+T01)', remove=F)

## ---- results='hide', warning=F, message=F---------------------------------
ctg_norm %>% filter(drug == 'DMSO') %>% 
  ggplot(aes(plate, viability, fill=line)) + geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ctg_norm %>% filter(drug == 'Bortezomib') %>% 
  ggplot(aes(plate, viability, fill=line)) + geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ctg_norm %>% mutate(type = ifelse(drug == 'DMSO', 'control', 'other')) %>% 
  ggplot(aes(pcount_norm, fill=type)) + geom_histogram(bins=40) + 
  theme_classic() + scale_fill_manual(values=c('#4285f4', '#dddddd'))

## ---- results='hide', warning=F, message=F---------------------------------
ctg_norm %>%
  mutate(type = ifelse(drug == 'DMSO', 'pos', 
                ifelse((drug == 'Bortezomib') & 
                       (rack.concentration %in% c('0.04', '1')), 'neg', 'other'))) %>% 
  filter(type != 'other', !is.na(viability)) %>% 
  group_by(plate, type) %>% 
  summarise(med_val = median(viability), 
            se = mad(viability)/sqrt(n()), 
            ymin = med_val - 2*se, 
            ymax = med_val + 2*se) %>% ungroup() %>% 
  ggplot(aes(plate, med_val, colour=type)) + geom_point() + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width=0) + 
  theme_classic() +
  theme(axis.text.x = element_blank()) + 
  ylab('Viability') + 
  scale_colour_manual(values=c('#4285f4', '#111111')) + 
  xlab('Plate')

## ---- results='hide', warning=F, message=F---------------------------------
ctg_norm %>% ggplot(aes(row, col, fill=viability)) +
  geom_raster(hjust=0, vjust=0) + 
  facet_wrap(~plate) + theme_bw() + 
  theme(panel.grid=element_blank()) + 
  scale_fill_gradientn(colors=colpal)

## ---- results='hide', warning=F, message=F---------------------------------
## we don't plot the acutal dates but a pseudo time instead
pseudo_time <- ctg_norm %>% distinct(date) %>% 
  arrange(date) %>% mutate(pseudo_time = 1:n())

## only mean
ctg_norm %>% inner_join(pseudo_time) %>%
  group_by(line, plate, pseudo_time) %>% 
  summarise(viability=mean(viability, na.rm=T)) %>% ungroup() %>% 
  ggplot(aes(pseudo_time, viability)) + 
  geom_point(size=2) + 
  geom_smooth(method='lm', colour='red') +
  theme_classic() + ylim(c(0,1)) + 
  xlab('Time') + ylab('Viability')

## ---- results='hide', warning=F, message=F---------------------------------
pca_res <- ctg_norm %>%
  acast(well ~ screen, value.var='viability', fun.aggregate = mean) %>% 
  .[,-1] %>% t() %>% .[,apply(., 2, function(x) sum(x) > 0)] %>%
  prcomp(scale=T, center=T)

## visualize
ctg_norm %>% distinct(screen, line, batch) %>% 
  inner_join(pca_res$x %>% as.data.frame() %>% 
               rownames_to_column('screen') %>% tbl_df %>% 
               dplyr::select(screen, PC1, PC2) %>% 
               mutate(screen = gsub('^X', '', gsub('\\.', '-', screen)))) %>%
  ggplot(aes(PC1, PC2)) + 
  geom_point(aes(colour=line), size=3) + 
  theme_classic()

## ---- results='hide', warning=F, message=F---------------------------------
## replicate annotation
reps <- ctg_norm %>%
  distinct(screen, line) %>% 
  arrange(screen, line) %>% 
  group_by(line) %>% mutate(rep=paste0('rep', 1:n())) %>% ungroup()

ctg_norm <- ctg_norm %>% inner_join(reps)

cor_scatter <- ctg_norm %>% distinct() %>%
  unite(barcode, line, well, sep='_') %>% 
  dplyr::select(barcode, rep, viability) %>% 
  spread(rep, viability) %>% drop_na() 

cor_scatter %>% summarise(PCC = cor(rep1, rep2, method='pearson'), 
                          SCC = cor(rep1, rep2, method='spearman'))

cor_scatter %>% 
  ggplot(aes(rep1, rep2)) + geom_point(alpha=0.4, stroke=0) + 
  geom_hline(yintercept = 1, linetype=2) + 
  geom_vline(xintercept = 1, linetype=2) + 
  geom_abline(slope = 1) + xlim(c(0,1.5)) + ylim(c(0,1.5)) +
  annotate('text', x=0.1, y=1.4, label='PCC = 0.87\nSCC = 0.74') +
  theme_classic() + xlab('Replicate 1') + ylab('Replicate 2')

## for individual oragnoid lines
ctg_norm %>%
  unite(barcode, line, well, sep='_', remove=F) %>% 
  dplyr::select(line, barcode, rep, viability) %>% 
  spread(rep, viability) %>% drop_na() %>% 
  ggplot(aes(rep1, rep2)) + geom_point(alpha=0.2, stroke=0) + 
  geom_hline(yintercept = 1, linetype=2) + 
  geom_vline(xintercept = 1, linetype=2) + 
  geom_abline(slope = 1) + xlim(c(0,1.5)) + ylim(c(0,1.5)) +
  facet_wrap(~line) + theme_bw() + theme(panel.grid=element_blank())

## correlation coefficients
ctg_norm %>% 
  unite(barcode, line, well, sep='_', remove=F) %>% 
  dplyr::select(line, barcode, rep, viability) %>% 
  spread(rep, viability) %>% drop_na() %>% 
  group_by(line) %>% summarise(cor = cor(rep1, rep2)) %>% ungroup() %>% 
  arrange(desc(cor)) %>% mutate(line = factor(line, levels=line)) %>% 
  ggplot(aes(line, cor)) + geom_bar(stat='identity') + 
  theme_classic() + ylim(c(0,1)) +
  ylab('Pearson correlation') + xlab('Organoid line') + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

## ---- results='hide', warning=F, message=F, eval=F-------------------------
#  auc_pgx <- ctg_norm %>%
#    mutate(pc_norm = viability) %>%
#    #filter(drug == coi) %>%
#    group_by(screen, drug) %>% arrange(rack.concentration) %>%
#    ## calculate auc both for fitted curve and actual values
#    mutate(s_AUC_fit = computeAUC(rack.concentration, pc_norm, verbose=T,
#                                viability_as_pct = F, area.type = 'Fitted'),
#           s_AUC_actual = computeAUC(rack.concentration, pc_norm, verbose=T,
#                                   viability_as_pct = F, area.type='Actual')) %>%
#    group_by(line, drug) %>% arrange(rack.concentration) %>%
#    ## calculate auc both for fitted curve and actual values
#    mutate(l_AUC_fit = computeAUC(rack.concentration, pc_norm, verbose=F,
#                                viability_as_pct = F, area.type = 'Fitted'),
#           l_AUC_actual = computeAUC(rack.concentration, pc_norm, verbose=F,
#                                   viability_as_pct = F, area.type='Actual')) %>%
#    ungroup()

## ---- results='hide', warning=F, message=F, eval=T, echo=F-----------------
data('auc_pgx_0318', package='SCOPEAnalysis')

## ---- results='hide', warning=F, message=F---------------------------------
## based on AUCs
auc_pgx %>% filter(drug %in% c('DMSO', 'Bortezomib')) %>% 
  distinct(line, drug, l_AUC_fit) %>% 
  ggplot(aes(line, 1-l_AUC_fit, colour=drug)) + geom_point() + theme_classic()

## ---- results='hide', warning=F, message=F---------------------------------
plot_drc <- function(dr, lines){
  ctg_norm %>% mutate(concentration = factor(rack.concentration,
                                             levels=c('0.0016', '0.008', '0.04', '0.2', '1'))) %>% 
    filter(drug == dr, line %in% lines) %>%
    ggplot(aes(concentration, viability)) + 
    geom_point() + 
    stat_summary(aes(group=line), fun.y=mean, 
                 geom="line", colour="#4285f4") +
    facet_wrap(~line) + ylim(c(0, 1.3)) + 
    theme_bw() + theme(panel.grid=element_blank())
}

## ---- results='hide', warning=F, message=F---------------------------------
auc_pgx %>% distinct(drug, line, l_AUC_fit) %>% 
  filter(!grepl('with|mit|DMSO', drug)) %>% 
  group_by(drug) %>% 
  mutate(l_AUC_fit = 1 - l_AUC_fit,
         AUC_centered = l_AUC_fit - median(l_AUC_fit)) %>% 
  ungroup() %>%
  acast(line ~ drug, value.var = 'AUC_centered') %>%
  pheatmap(border_color = '#ffffff')

## ---- results='hide', warning=F, message=F---------------------------------
data('aucs_image', package='SCOPEAnalysis')

## ---- results='hide', warning=F, message=F---------------------------------
aucs_image %>% distinct(Line, Product.Name, l_AUC_fit) %>% tbl_df %>% 
  distinct() %>% 
  filter(!grepl('mit|with', Product.Name)) %>%
  group_by(Product.Name) %>% 
  mutate(AUC_centered = l_AUC_fit - median(l_AUC_fit)) %>%
  # ungroup() %>%
  acast(Line ~ Product.Name, value.var = 'AUC_centered') %>% 
  pheatmap(border_color = '#ffffff')

## ---- results='hide', warning=F, message=F---------------------------------
## load mutations
data('mutations', package='SCOPEAnalysis')

## binary mutation profile
mut_profile <- mutations %>% dplyr::select(line, SYMBOL) %>% 
  mutate(mut = as.factor(1))  %>% distinct() %>% 
  spread(SYMBOL, mut) %>% mutate_all(funs(ifelse(is.na(.), 0, .))) %>%
  rbind(c('D018T01', rep(0, 19)))

## ---- results='hide', warning=F, message=F---------------------------------
## annotate mutation status
drugmut <- distinct(auc_pgx, line, drug, l_AUC_fit) %>% 
  left_join(mut_profile) %>% mutate_at(vars(AKT1:TP53), 'as.factor') %>%
  mutate(AUC = 1-l_AUC_fit)

## add clinical data that might influence the results
data('clindat', package='SCOPEAnalysis')

drugmut <- drugmut %>% inner_join(clindat)

## ---- results='hide', warning=F, message=F---------------------------------
## form combinations
dm_combis <- expand.grid(
    unique(drugmut$drug[!drugmut$drug %in% c('Staurosporine_500nM', 'DMSO')]),
    colnames(mut_profile)[-1]
) %>% tbl_df %>% `colnames<-`(c('drug', 'genotype'))
dm_combis <- dm_combis %>% filter(!grepl('with|mit', drug))

## only genes that are mutated at least twice
dm_combis <- dm_combis %>% 
  left_join(mut_profile %>% dplyr::select(-line) %>%
              summarise_all(funs(sum(as.integer(.)))) %>%
              gather(genotype, mut_count)) %>% 
  filter(mut_count > 2) %>% dplyr::select(-mut_count)

## function for modelling/testing
test_druggene <- function(df){
  d <- as.character(df$drug[1])
  m <- as.character(df$genotype[1])
  
  ## get the right data
  data <- drugmut %>% filter(drug == d) %>% 
    dplyr::select(line, AUC, matches(paste0('^', m, '$')), everything())
  colnames(data)[3] <- 'mut'
  
  
  if(sum(as.numeric(as.character(data$mut))) == 0){
    return(NA)
  } else {
    ## models
    permts <- permTS(AUC ~ as.factor(mut), data = data,
                     method="exact.mc",
                     control=permControl(nmc=10^4-1))
    
    return(tibble(
      drug = d,
      mutation = m,
      p.value = permts$p.value,
      dAUC = permts$estimate,
      variance = var(data$AUC)
    ))
  }
}

dgi <- dm_combis %>% rowwise() %>% do(results=test_druggene(.)) %>%
  unnest(results) %>% mutate(padj = p.adjust(p.value, method='BH')) %>%
  arrange(p.value)

## ---- results='hide', warning=F, message=F---------------------------------
## blue/yellow colour-pallette
dgi_palette <- colorRampPalette(c('#2b8cbe', '#f7f7f7', '#f7f7f7', '#fed98e'))(7)

dgi %>% mutate(visval = ifelse(dAUC < 0, (-1) * (-log10(p.value)), -log10(p.value))) %>%
  dplyr::select(drug, mutation, visval) %>% 
  spread(drug, visval) %>% data.frame() %>% `rownames<-`(NULL) %>% 
  column_to_rownames('mutation') %>% 
  pheatmap(color=dgi_palette, border_color = '#ffffff')

## ---- results='hide', warning=F, message=F---------------------------------
int_nutlin <- drugmut %>% filter(drug == 'Nutlin3a') %>%
  mutate(TP53 = ifelse(TP53 == '0', 'WT', 'Mut.')) %>%
  ggplot(aes(TP53, AUC)) + geom_jitter(width=0.1) +
  stat_summary(fun.y='mean', fun.ymax = 'mean', 
               fun.ymin = 'mean', geom='crossbar', colour='red',
               width=0.5) + 
  theme_classic() + 
  labs(title='Treatment with Nutlin3a',
       subtitle = 'P = 0.002 (21% FDR)')

int_gefitinib <- drugmut %>% filter(drug == 'Gefitinib') %>%
  mutate(NRAS = ifelse(NRAS == '0', 'WT', 'Mut.')) %>%
  ggplot(aes(NRAS, AUC)) + geom_jitter(width=0.1) +
  stat_summary(fun.y='mean', fun.ymax = 'mean', 
               fun.ymin = 'mean', geom='crossbar', colour='red',
               width=0.5) + 
  theme_classic() + 
  labs(title='Treatment with Gefitinib',
       subtitle = 'P = 0.0008 (21% FDR)')

int_nutlin + int_gefitinib

## ---- results='hide', warning=F, message=F---------------------------------
aucs <- auc_pgx %>% group_by(drug, line, rep) %>% 
  summarise(AUC = 1-mean(l_AUC_fit)) %>% ungroup() %>% 
  group_by(drug) %>% mutate(AUC_centered = AUC - median(AUC)) %>% ungroup()

## ---- results='hide', warning=F, message=F---------------------------------
data('promise_expr', package='SCOPEAnalysis')

## ---- results='hide', warning=F, message=F---------------------------------
## dynamic phenotype range for each drug
drugs_dr <- aucs %>% group_by(drug, line) %>% 
  filter(!grepl('with|mit|Staurosporine_500nM', drug)) %>%
  summarise(AUC = mean(AUC)) %>% ungroup() %>% 
  group_by(drug) %>% 
  mutate(dynamic_range = range(AUC) %>% reduce(`-`) %>% abs()) %>% 
  ungroup() %>% distinct(drug, dynamic_range)

## visualize
drugs_dr %>% arrange(dynamic_range) %>% 
  mutate(rank = 1:n(), 
         label = ifelse(drug %in% c('DMSO', 'Gefitinib', '5-FU', 'Bortezomib'), drug, '')) %>% 
  ggplot(aes(rank, dynamic_range)) + 
  geom_point() + geom_line() +
  geom_text_repel(aes(label=label), nudge_y = 0.05) + 
  geom_hline(yintercept = 0.21, colour = 'red') + 
  theme_classic() + 
  ylab('Dynamic range of drug effect (AUCmax - AUCmin)') + 
  xlab('Drug rank')

## select drugs to include in analysis
drugs_include <- drugs_dr %>% filter(dynamic_range > 0.21) %>% .$drug

## ---- results='hide', warning=F, message=F---------------------------------
## expr matrix
expr_mat <- promise_expr %>% 
  reshape2::acast(probe ~ id, value.var = 'expr', fun.aggregate =mean) 

## drug response
pb <- progress_estimated(length(drugs_include))
expr_interactions <- aucs %>% filter(drug %in% drugs_include) %>%
  group_by(drug, line) %>% 
  summarise(AUC = mean(AUC)) %>% ungroup() %>% 
  arrange(line) %>% nest(-drug) %>% 
  mutate(res = map(data, ~{
    pb$tick()$print()    
    if(!identical(.x$line, colnames(expr_mat))){
      mm <- model.matrix(~AUC, data = bind_rows(.x, .x) %>% arrange(line))
      tidy(eBayes(lmFit(expr_mat, mm), robust=T))
    } else {
      stop('Sample frame does not fit to expression matrix')
    }
  })) %>% unnest(res) %>% arrange(p.value) %>% 
  dplyr::select(probe=gene, everything()) %>% 
  inner_join(promise_expr %>% distinct(probe, symbol)) %>%
  mutate(fdr=p.adjust(p.value, method='BH'))

## ---- results='hide', warning=F, message=F---------------------------------
drug_expr_volcano <- function(d, highlight){
  df <- expr_interactions %>% filter(drug == d) %>% 
    group_by(drug, symbol) %>% top_n(1, desc(p.value)) %>% 
    ungroup()
  
  ## fdr 0.05
  co <- df %>% filter(fdr < 0.05) %>% arrange(desc(fdr)) %>% 
    .$p.value %>% .[1] %>% log10() %>% `*`(-1)
  
  df %>% mutate(label = ifelse(symbol %in% highlight, symbol, ''), 
         type = ifelse((estimate < 0) & (fdr < 0.05), 'sensitizer', 
                ifelse((estimate > 0) & (fdr < 0.05), 'resistance', 'none'))) %>% 
  ggplot(aes(estimate, -log10(p.value))) + 
  geom_hex(aes(fill = type, colour=type), bins=100) + theme_classic() + 
  scale_fill_manual(values = c('#dddddd', '#FED98E', '#5087C8')) + 
  scale_colour_manual(values = c('#dddddd', '#FED98E', '#5087C8')) + 
  geom_text_repel(aes(label = label)) + 
  geom_hline(yintercept = co, linetype = 'dashed') + 
  ggtitle(paste('Expression markers of', d, 'response')) +
  xlab('Model coefficient')
}

## EGFR Inhibitor Pathway; q-value = 0.002; pathway source = PharmGKB; tool = ConsensusPathDB
gefitinib_hits <- c('EGF', 'CAMK2D', 'CAMK1D', 'SHC1', 'STAT3',
                    'PRKCA', 'GRB2', 'CAMK1', 'MAPK1', 'STAT5A')
drug_expr_volcano('Gefitinib', gefitinib_hits)

# TP53 Regulates Transcription of Death Receptors and Ligands; 
## pathway source = Reactome; q-value = 0.004; tool = ConsensusPathDB
nutlin_hits <- c('LIF', 'PML', 'FAS', 'RPS27L', 'DDB2', 'SERPINE1', 
                 'TIGAR', 'TNFRSF10A', 'TNFRSF10B', 'TNFRSF10C', 
                 'TNFRSF10D', 'PLK3', 'MDM2')
drug_expr_volcano('Nutlin3a', nutlin_hits)

## ---- results='hide', warning=F, message=F---------------------------------
data('aucroc', package = 'SCOPEAnalysis')

## ---- results='hide', warning=F, message=F---------------------------------
drug_morph_interaction <- mut_profile %>% 
  gather(mutation, status, -line) %>% 
  group_by(mutation) %>% filter(sum(as.integer(status)) >=3) %>%
  nest(-mutation) %>%
  mutate(tres = map(data, ~{
    .x %>% inner_join(aucroc) %>% nest(-drug) %>% 
      mutate(res = map(data, 
                       ~ broom::tidy(t.test(auroc ~ status, data = .x)))) %>% 
      unnest(res)
  })) %>% unnest(tres)

drug_morph_interaction <- drug_morph_interaction %>% 
  arrange(p.value) %>% mutate(FDR = p.adjust(p.value, method='BH'))

## ---- results='hide', warning=F, message=F, eval=F-------------------------
#  ## load feature angles which are the edges of the network
#  data('feature_angles', package='SCOPEAnalysis')
#  
#  ## add drug-genotype interactions
#  feature_angles <- drug_morph_interaction %>%
#    dplyr::select(source=mutation, target=drug, p.value) %>%
#    bind_rows(feature_angles)
#  
#  ## select drugs with at least one interaction
#  good_drugs <- feature_angles %>% filter((angle < 45) | (p.value < 0.05)) %>%
#    distinct(source, target)
#  good_drugs <- unique(c(good_drugs$source, good_drugs$target))
#  
#  ## select mutations/drugs to be labeled
#  for_labeling <- c('APC', 'PIK3CA', 'KRAS', 'TP53',
#                    'Tyrphostin AG 879', 'K02288',
#                    'WIKI4', 'AT7867')
#  
#  ## make a tidy graph
#  graph <- feature_angles %>%
#    filter((angle < 45) | (p.value < 0.05),
#           source != target,
#           source != 'NRAS') %>%
#    as_tbl_graph() %>%
#    activate(nodes) %>%
#    mutate(type = ifelse(name %in% drug_morph_interaction$mutation, 'mutation', 'drug')) %>%
#    activate(nodes) %>%
#    filter(name %in% good_drugs) %>%
#    mutate(label = ifelse(name %in% for_labeling, name, ''),
#           index = 1:n(),
#           has_edge = ifelse((index %in% .E()$from) |
#                             (index %in% .E()$to), T, F),
#           size = ifelse(type == 'mutation', 2, 1))
#  
#  ## draw the graph
#  graph %>% ggraph(layout = 'fr') +
#    geom_edge_fan(colour = '#dddddd') +
#    geom_node_point(aes(shape = type, colour = type, size = size)) +
#    geom_node_text(aes(label = label), nudge_y = 0.2) +
#    theme_graph() +
#    scale_colour_manual(values = c('#aaaaaa', '#4285f4'))

## --------------------------------------------------------------------------
sessionInfo()

