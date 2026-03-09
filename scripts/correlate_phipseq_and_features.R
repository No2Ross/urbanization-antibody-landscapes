suppressPackageStartupMessages({
    #library(MOFA2)
    library(tidyverse)
    library(here)
    library(glue)
    library(ggpubr)
    library(cowplot)
  library(patchwork)
})
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")



PATHS <- list(
    data = list(
        EXV = "dataCell/urban-rural-study-manuscript/EXV_percentLineage.rds",
        PMA = "dataCell/urban-rural-study-manuscript/data/PMA_percentLineage.rds",
        MPL = "dataCell/urban-rural-study-manuscript/data/MPL_percentLineage.rds",
        MET = "dataCell/urban-rural-study-manuscript/data/MET_80ptile.rds",
        GLY = "dataCell/urban-rural-study-manuscript/data/GLY_igg.rds"
    ),
    meta = "dataCell/urban-rural-study-manuscript/data/PMA_metadata.rds",
    results = "dataCell/urban-rural-study-manuscript/outputs/",
    models = "dataCell/urban-rural-study-manuscript/outputs/mofa/trained_models/"
)

color_groups <- c("#E69F00", "#56B4E9", "#009E73")
color_groups <- purrr::set_names(color_groups, c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))
theme_set(theme_pubr(10, legend="right"))

# for scaling
tst <- function(matrix) t(scale(t(matrix)))

## PARAMETERS
prop_top_features <- 1
mofa_iter <- 50

# load and prepare data ------------------------------------------------------------

setwd("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq")

mofa <- readRDS("dataCell/urban-rural-study-manuscript/outputs/mofa/final_mofa_model.rds")

mofa.data <- MOFA2::impute(mofa)
mofa.data <- MOFA2::get_data(mofa.data, as.data.frame=TRUE)
mofa.data <- mofa.data |>
    pivot_wider(names_from=c(view, feature), values_from=value)

luminex <- readRDS("dataCell/urban-rural-study-manuscript/data/data_luminex.rds") |>
    select(donor, ccl24_eotaxin2_mpif2:cd40ligand_tnfsf5)
colnames(luminex)[2:10] <- paste0("luminex_", colnames(luminex)[2:10])

# get PCA coordinates ---------------------------------

path <- "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/dataCell/processed/PCoA_coords.csv"
pc.coords <- data.table::fread(path)


combined.dfs <- right_join(mofa.data, pc.coords, by=c("sample"="sampleID")) |>
    #drop_na(Axis.1, exvivo_gdT_DN_CD45RApos_CXCR3pos) |>
    left_join(luminex, by=c("sample"="donor")) #|>
    #drop_na(ccl24_eotaxin2_mpif2)

dim2 <- combined.dfs$Axis.2
# features <- combined.dfs[,c(3:513)]
features <- combined.dfs[,str_detect(colnames(combined.dfs), "metflow", negate = TRUE)] # exclude metflow data
# metflows <- combined.dfs[,170:505]

features <- as.data.frame(features)
row.names(features) <- features$sample
features$sample <- NULL
features$group <- NULL

#Remove individuals with absolutely no information or barely any information (i.e. one person is missing 176 values)
features <- features[rowSums(is.na(features) == TRUE) != 185 | rowSums(is.na(features) == TRUE) != 176,]

#Remove non-numeric columns
features <- features %>%
                select(where(is.numeric))

results <- list()
for (i in 1:ncol(features)) {
  
    results[[i]] <- cor.test(dim2, features[[i]], method="spearman")
}

estimates <- sapply(results, function(x) x$estimate)
pvalues <- sapply(results, function(x) x$p.value)

df.results <- data.frame(
    feature = colnames(features),
    estimate = estimates,
    pvalue = pvalues
)
df.results$fdr <- p.adjust(df.results$pvalue, method="BH")
df.results <- df.results |>
    mutate(data = str_extract(feature, "exvivo|pma|mpl|metflow|glyco|luminex")) |>
    mutate(fdr_perdata = p.adjust(pvalue, method="BH"), .by=data)

#Remove the MDS axis as features
df.results <- df.results[str_detect(df.results$feature, "exvivo|pma|mpl|metflow|glyco|luminex"),]

p1 <- df.results |>
    ggplot(aes(estimate, -log10(fdr))) +
    geom_point() +
    geom_point(data = df.results |>
                   filter(fdr<.1) |>
                   mutate(Direction = ifelse(estimate>0, "Urban Dutch", "Rural Senegalese")),
               aes(color=Direction)) +
    geom_hline(yintercept=-log10(0.1), lty=2) +
    labs(x="Correlation with urban-rural metric (Spearman's rho)", y="-log10(FDR)") +
    scale_color_manual(values=color_groups[c(1,3)]) +
    cowplot::theme_cowplot(10)
p1
ggsave("plotsCell/correlation_volcano.pdf", p1, dpi=600, width=5, height=3)

p2 <- df.results |>
    filter(fdr<0.1) |>
    mutate(direction = sign(estimate)) |>
    mutate(Direction = ifelse(estimate<0, "Urban Dutch", "Rural Senegalese")) |>
    mutate(feature = fct_reorder(feature, estimate)) |>
    ggplot(aes(estimate, feature, fill=Direction)) +
    geom_vline(xintercept = 0, lty=1) +
    geom_col() +
    labs(x="Correlation with MDS2", y=NULL) +
    scale_fill_manual(values=color_groups[c(1,3)]) +
    cowplot::theme_minimal_grid(10) +
    scale_y_discrete(position="right") +
    theme(legend.position = "none")
p2
ggsave("plotsCell/correlation_features.pdf", p2, dpi=600, width=6, height=3)


# scatter plots -----------------------------------------------------------

sigs <- df.results |>
    filter(fdr<0.1) |>
    arrange(estimate) |>
    pull(feature)

plots <- list()
for (feature in sigs){
    rho <- round(df.results$estimate[df.results$feature==feature],3)
    fdr <- round(df.results$fdr[df.results$feature==feature], 3)
    plots[[feature]] <- combined.dfs |>
        ggplot(aes(Axis.2, .data[[feature]])) +
        geom_point(aes(color=Residence, shape=Residence, fill=Residence))  +
        geom_smooth(method="lm", se=F) +
        labs(x="PhIP-Seq MDS2", y="Feature", title=feature,
             subtitle=glue::glue("Spearman's rho={rho}, FDR={fdr}")) +
        scale_color_manual(values=color_groups) +
        scale_fill_manual(values=color_groups)  +
        scale_shape_manual(values=c(21:23)) +
        cowplot::theme_half_open(10) +
        theme(legend.position = "none")
}
p3 <- wrap_plots(plots, nrow=3)
p3
ggsave("plotsCell/correlation_scatterplots.pdf", p3, dpi=600, width=10, height=7)


apply(combined.dfs[,sigs], 2, function(x) sum(!is.na(x)))
