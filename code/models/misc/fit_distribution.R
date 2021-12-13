library(tidyverse)
library(here)
library(monocle3)
library(ggrastr)
library(cowplot)
library(princurve)
library(scico)

library(fitdistrplus)
set.seed(123)

# I wish I could solve my path problems with the here package. 
PATH = "/dkfz/groups/shared/OE0049/B110-Isilon2/promise/"

# Loading data and formatting for easier access
umap_tidy <- read_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_tidy.Rds"))

df <- umap_tidy %>% filter(drug == "DMSO")

plotdist(df$size, demp = TRUE)
descdist(df$size, boot = 1000)

fw <- fitdistrplus::fitdist(df$size, "weibull")
#fg <- fitdist(df$size, "gamma")
fln <- fitdistrplus::fitdist(df$size, "lnorm")
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal")

pdf(paste0(PATH, "reports/figures/all_lines_dist.pdf"))
denscomp(list(fw, fln), legendtext = plot.legend)
qqcomp(list(fw, fln), legendtext = plot.legend)
cdfcomp(list(fw, fln), legendtext = plot.legend)
ppcomp(list(fw, fln), legendtext = plot.legend)
dev.off()

# I run the same distribution fit on each line independently
# function definition
fitdist <- function(data, loi){
  df <- data %>% dplyr::filter(line == loi)
  
  pdf(paste0(PATH, "reports/figures/", loi, "_dist.pdf"))
  plotdist(df$size, demp = TRUE)
  descdist(df$size, boot = 1000)
  
  fw <- fitdistrplus::fitdist(df$size, "weibull")
  #fg <- fitdist(df$size, "gamma")
  fln <- fitdistrplus::fitdist(df$size, "lnorm")
  par(mfrow = c(2, 2))
  plot.legend <- c("Weibull", "lognormal")
  
  denscomp(list(fw, fln), legendtext = plot.legend)
  qqcomp(list(fw, fln), legendtext = plot.legend)
  cdfcomp(list(fw, fln), legendtext = plot.legend)
  ppcomp(list(fw, fln), legendtext = plot.legend)
  dev.off()
}

# running function
df <- colData(obj) %>% as_tibble() %>% filter(drug == "DMSO")
lapply(df$line %>% unique(), fitdist, data = df)


# count distribution
df <- colData(obj) %>% as_tibble() %>% filter(drug == "DMSO") %>%
  dplyr::count(concentration, line, replicate, well)

plotdist(df$n, demp = TRUE)
descdist(df$n, boot = 1000)

fw <- fitdistrplus::fitdist(df$n, "weibull")
#fg <- fitdist(df$size, "gamma")
fln <- fitdistrplus::fitdist(df$n, "lnorm")
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal")

pdf(paste0(PATH, "reports/figures/count_all_lines_dist.pdf"))
denscomp(list(fw, fln), legendtext = plot.legend)
qqcomp(list(fw, fln), legendtext = plot.legend)
cdfcomp(list(fw, fln), legendtext = plot.legend)
ppcomp(list(fw, fln), legendtext = plot.legend)
dev.off()
