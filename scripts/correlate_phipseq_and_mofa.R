suppressPackageStartupMessages({
    library(tidyverse)
    library(here)
    library(glue)
    library(patchwork)
    library(Hmisc)
    library(nlme)
    library(emmeans)
    library(PMCMRplus)
    library(princurve)
    library(ggpubr)
    library(ggbeeswarm)
    library(vegan)
    library(ggnewscale)
})

PATHS <- list(
    latent_factors = "dataCell/urban-rural-study-manuscript/outputs/mofa/final_mofa_model_factors.rds",
    meta = "dataCell/urban-rural-study-manuscript/data/PMA_metadata.rds",
    mofa_model = "dataCell/urban-rural-study-manuscript/outputs/mofa/final_mofa_model.rds"
)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(stats::predict)

color_groups <- c("#E69F00", "#56B4E9", "#009E73")
color_groups <- purrr::set_names(color_groups, c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))
theme_set(theme_pubr(10, legend="right"))



# load data ---------------------------------------------------------------

setwd("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq")

meta <- readRDS(PATHS$meta) %>%
    rename(sample = patient_id) %>%
    as.data.frame() %>%
    mutate(group3 = factor(group3,
                           levels=c("RUR", "URB", "EUR"),
                           labels=c("Rural Senegalese", "Urban Senegalese", "Urban Dutch")))
rownames(meta) <- meta$sample

# load mofa latent factors
latent.factors <- readRDS(PATHS$latent_factors)



# correlate with correspondence analysis coordinates ----------------------

path <- "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/dataCell/processed/PCoA_coords.csv"
pc.coords <- data.table::fread(path)


pc.coords <- pc.coords |>
    right_join(latent.factors, by=c("sampleID"="sample")) |>
    drop_na() |>
    mutate(Residence = factor(Residence, levels=c("Rural Senegalese", "Urban Senegalese", "Urban Dutch")))

centroids <- pc.coords %>%
    group_by(Residence) %>%
    dplyr::summarize(center_x = mean(Axis.1), center_y = mean(Axis.2))

pc.coords |>
    ggplot(aes(Axis.1, Axis.2, color=Residence)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    ggpubr::stat_stars() +
    geom_point(aes(shape=Residence, fill=Residence)) +
    geom_point(data = centroids, aes(center_x, center_y, fill=Residence, shape=Residence), color="black", cex=3) +
    scale_color_manual(values=color_groups) +
    scale_fill_manual(values=color_groups) +
    cowplot::theme_minimal_grid(font_size=10) +
    scale_shape_manual(values=c(21:23))
ggsave("plotsCell/PCoA_coords_cytof_samples.pdf", width=5, height=3)

# principal curve ---------------------------------------------------------

latent.factors <- readRDS("dataCell/urban-rural-study-manuscript/outputs/mofa/final_mofa_model_factors.rds")
latent.factors.coords <- latent.factors[, grepl("Factor", colnames(latent.factors))]

mds.res <- cmdscale(dist(latent.factors.coords))
latent.factors$mds1 <- mds.res[,1]
latent.factors$mds2 <- mds.res[,2]

set.seed(42)
pcurve <- principal_curve(as.matrix(latent.factors[,c("mds1", "mds2")]), maxit=100, thresh = 0.0005)
latent.factors$lambda <- pcurve$lambda

lf.reord <- latent.factors
ix <- order(lf.reord$lambda)
pcurve.coords <- pcurve$s[ix,]
lambda <- lf.reord$lambda[ix]

pc.coords |>
    right_join(latent.factors, by=c("sampleID"="sample")) |>
    drop_na() |>

    ggplot(aes(Axis.2, lambda)) +
    geom_point(aes(color=Residence, shape=Residence, fill=Residence))  +
    # ggpubr::stat_stars(aes(color=Residence)) +
    stat_smooth(method="lm", se=FALSE) +
    ggpubr::stat_cor(method="spearman") +
    scale_color_manual(values=color_groups) +
    scale_fill_manual(values=color_groups)  +
    scale_shape_manual(values=c(21:23)) +
    labs(x="PhIP-seq MDS2", y="Rural-Urban Trajectory Metric") +
    cowplot::theme_half_open(10)
ggsave("plotsCell/phipseq_lambda_correlations.pdf", dpi=600, width=4.2, height=3)

