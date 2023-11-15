## Saunders, Srivatsan, et al. Nature (2023) - Embryo-scale reverse genetics at single-cell resolution
## This script contains code to generate plots for Figure 5
## from processed data files available in the GEO repository GSE202639 or available via Github.

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(viridis)
  library(data.table)
  library(tidyr)
  library(monocle3)
  library(plotly)
  library(htmlwidgets)
  library(VGAM)
  library(forcats)
  
})

source("saunders_srivatsan_2023_utils.R")

# Load data ---------------------------------------------------------------

# notochord cells from ctrls, tbxta, noto and smo
notochord_cds <- readRDS("data/ctrl_tbxta-noto-smo_notochord_20k_cds.RDS")

# export umap coords to coldata
notochord_cds <- append_umap_coordinates(notochord_cds, 
                                        umap_3D = TRUE)

# extract cell metadata
notochord_coldata <- colData(notochord_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(timepoint))

# collapse all ctrl genos into a new column
notochord_coldata <- notochord_coldata %>%
  mutate(gene_target_test = ifelse(grepl("ctrl", gene_target), "ctrl", gene_target))

# add the test col back to the cds
colData(notochord_cds)$gene_target_test <- notochord_coldata$gene_target_test


# Figure 5a ---------------------------------------------------------------

## plot celltype abundance heatmap for tbxta and noto relative to controls

# load abundance test data
res_df <- fread("data/zperturb_beta-bin_results_broadtype_clean.csv", 
               sep = ",", data.table = F, stringsAsFactors = F, na.strings = "")

# filter results
qval_thresh <- 0.01
sig_df <- res_df %>% 
  filter(geno.v.wt %in% c("tbxta", "noto")) %>% 
  dplyr::mutate(sig_fc = case_when(qval < qval_thresh ~ TRUE, 
                                   TRUE ~ FALSE))

# assign meso celltypes
grps <- c("notochord", "mesodermal progenitor cells (contains PSM)", 
         "mature slow muscle", "slow-committed myocyte",
         "mature fast muscle", "floor plate", 
         "fast-committed myocyte", "hatching gland")

# summarize fc in a matrix
hm_wide <- sig_df %>%
  filter(cell_group %in% grps) %>% # use for all cell groups
  select(geno.v.wt, timepoint, cell_group, abund_log2fc) %>%
  pivot_wider(
    names_from = cell_group,
    values_from = abund_log2fc,
    values_fill = c(0)
  )

hm_plot <- hm_wide %>%
  pivot_longer(!geno.v.wt:timepoint,
               names_to = "cell_group",
               values_to = "plot") %>%
  left_join(
    sig_df %>%
      select(geno.v.wt, timepoint, cell_group, sig_fc),
    by = c("geno.v.wt", "timepoint", "cell_group"),
    
  ) %>%
  mutate(
    sig_fc = replace_na(sig_fc, FALSE),
    geno_timepoint = paste0(geno.v.wt, "_", timepoint)
  ) %>%
  arrange(timepoint, geno.v.wt) %>%
  mutate(timepoint = as.factor(timepoint),
         cell_group = as.factor(cell_group))

new_levels <- hm_plot %>% distinct(geno_timepoint) %>% pull %>% as.character
hm_plot$geno_timepoint <- factor(hm_plot$geno_timepoint, levels = new_levels)

# get clustering results
hm_wider <- hm_wide %>%
  mutate(geno_timepoint = paste0(geno.v.wt, "_", timepoint)) %>%
  select(geno_timepoint, everything(),-geno.v.wt,-timepoint)
hm_mat <- as.matrix(hm_wider[,-1])
rownames(hm_mat) <- hm_wider$geno_timepoint
hm_mat <- t(hm_mat)
hm.dendro <- as.dendrogram(hclust(d = dist(x = hm_mat)))
hm.order <- order.dendrogram(hm.dendro)

# reorder cell type factor
hm_plot$cell_group <- factor(as.character(hm_plot$cell_group),
                             rownames(hm_mat[hm.order, ]))

# make df for plotting significance boxes
fms <- hm_plot[hm_plot$sig_fc, c("geno_timepoint", "cell_group")]
fms$geno_timepoint <- as.integer(fms$geno_timepoint)
fms$cell_group <- as.integer(fms$cell_group)

# plot heatmap
ggplot(hm_plot,
       aes(x = geno_timepoint, y = cell_group, fill = plot)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#3234a9",
    mid = "white",
    high = "darkred",
    na.value = "black",
    name = ""
  ) +
  geom_rect(
    data = fms,
    size = 0.5,
    fill = NA,
    colour = "black",
    aes(
      xmin = geno_timepoint - 0.5,
      xmax = geno_timepoint + 0.5,
      ymin = cell_group - 0.5,
      ymax = cell_group + 0.5
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    text = element_text(size = 10),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.75),
    plot.margin = unit(c(.5, .5, .5, 2.5), "cm"),
    legend.key.size = unit(0.5, 'cm'),
    axis.line = element_line(colour = "black", linewidth = 1),
    axis.ticks = element_line(colour = "black", linewidth = 0.25)
  ) +
  coord_equal() +
  monocle3:::monocle_theme_opts()

ggsave("fig5a_tbxta-noto_abund_heatmap.pdf",
       units = "in",
       width = 6,
       height = 4)


# Figure 5b ---------------------------------------------------------------

# load normalized counts/embryo
summary_df <- fread("data/zperturb-all_broadtype_summary_merge.csv", 
                      sep = ",", data.table = F, stringsAsFactors = F) # cell_type_broad

# Calculate cells per thousand
summary_df <- summary_df %>% 
  mutate(cpt = round(cells / total_cells * 1000), 
         timepoint = as.character(timepoint))

# select genotypes, timepoints and tissues 
types <- c("notochord")
genos <- c("ctrl-inj", "tbxta", "noto")
times <- c("18", "24", "36")
clrs <- c("#009999", "#F16667", "#DE281A")

summary_df %>% 
  filter(genotype %in% genos) %>% 
  filter(cell_type %in% types) %>%
  filter(timepoint %in% times) %>% 
  ggplot(aes(x = timepoint, 
             y = cpt, 
             fill = genotype)) +
  geom_boxplot(color = "black", size = 0.3, outlier.size = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.5, color = "black", width = 0.1) +
  theme(axis.title.x = element_blank(), legend.position = "none", 
        axis.ticks = element_line(color = "black", size = .3),
        legend.title = element_blank(), 
        axis.text.x = element_text(angle= 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black")) +
  facet_wrap(~genotype, nrow = 1) +
  scale_fill_manual(values = clrs) +
  monocle3:::monocle_theme_opts()

ggsave("fig5b_ctrl-noto-tbxta_notochord_counts-boxplot_time-facet.pdf",
       dpi = 500, width = 4, height = 3)

# Figure 5c ---------------------------------------------------------------

## Plot notochord umap colored by timepoint
# get timepoint colors
rainbow_timepoint_colors = 
  c("18h" = "#DF4828", 
    "24h" = "#E78C35", 
    "36h" = "#F6C141", 
    "48h" = "#4EB265", 
    "72h" = "#1965B0")

num_colors <- notochord_coldata %>%
  pull(timepoint) %>%
  unique() %>%
  sort() %>%
  length()

full_spectrum_timepoint <- colorRampPalette(rainbow_timepoint_colors)(num_colors) 

names(full_spectrum_timepoint) <- notochord_coldata %>%
  filter(!is.na(timepoint)) %>%
  pull(timepoint) %>%
  unique() %>%
  sort()

# plot timepoint in 2D
ggplot(notochord_coldata %>%
         sample_n(size = dim(notochord_coldata)[1])) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_3),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_3,
                 color = timepoint),
             stroke = 0,
             size = 0.6) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = full_spectrum_timepoint) 
ggsave("fig5c_all-notochord_ctrl-tbxta-noto-smo_timepoint_umap.png",
         dpi = 750,
         height = 2,
         width = 2.5,
         bg = "transparent")

# Figure 5d ----------------------------------------------------

# plot tbxta cells in red
ctrl_tbxta_coldata <- notochord_coldata %>%
  filter(gene_target_test %in% c("ctrl", "tbxta"))

# plot timepoint in 2D
ggplot(ctrl_tbxta_coldata) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_3),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_3,
                 color = gene_target_test),
             stroke = 0,
             size = 0.6) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey90", "tomato"))
ggsave("fig5d_all-notochord_ctrl-tbxta_geno_umap.png",
       dpi = 750,
       height = 2,
       width = 2.5,
       bg = "transparent")

# Figure 5e ----------------------------------------------

# genes to plot
nt_genes <- c("epyc", "col2a1a", "sox5", "sox6", "shha")

rowData(notochord_cds) %>% 
  as.data.frame() %>% 
  filter(grepl("sox", gene_short_name)) %>% 
  arrange(-num_cells_expressed) %>% 
  head(20)

# plot with legend
plot_genes_by_group(notochord_cds[,colData(notochord_cds)$gene_target_test %in% c("ctrl", "tbxta")], 
                    group_cells_by = "gene_target_test", 
                    markers = nt_genes, 
                    axis_order = "marker_group", max.size = 6) +
        ggplot2::coord_flip() +
        theme(axis.text = element_text(size = 6, color = "black"), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6))

ggsave("fig5e_ctrl_tbxta_notochord_expr_dotplot_hm_wlegend.pdf",
       width = 3, height = 2.5)

 # EDF20 Parachordal cartilage expr ----------------------------------------------------

plot_cells(notochord_cds[,colData(notochord_cds)$gene_target_test %in% c("ctrl", "tbxta")], x = 1, y = 3,  genes = c("epyc"), label_cell_groups = F) +
  facet_wrap(~gene_target_test) +
  scale_color_viridis_c() +
  theme_void() +
  theme(legend.position = "none")
ggsave("fig_plots/all-notochord_ctrl-tbxta_epyc-expr_umap_fig4.png",
       dpi = 750,
       height = 2,
       width = 5,
       bg = "transparent")

# plot tgm2l expression
plot_cells(notochord_cds[,colData(notochord_cds)$gene_target_test %in% c("ctrl")], x = 1, y = 3,  genes = c("tgm2l"), label_cell_groups = F) +
  scale_color_viridis_c() +
  theme_void() +
  theme(legend.position = "none")
ggsave("plots/all-notochord_ctrl_tgm2l-expr_umap_fig4.png",
       dpi = 750,
       height = 2,
       width = 3,
       bg = "transparent")

# EDF20a Parachordal Cartilage differential gene expression -----------------------------------

# subset data for a pairwise comparison of the proper timepoint, cells and genotypes
comp_cds <- notochord_cds[,colData(notochord_cds)$timepoint %in% c("36") &
                           colData(notochord_cds)$gene_target %in% c("ctrl-inj", "tbxta") &
                           !(clusters(notochord_cds) %in% c(2, 5, 7, 8))]

# perform DE testing
comp_cds <- comp_cds[Matrix::rowSums(SingleCellExperiment::counts(comp_cds) > 0) > 10,] # filter for expressed genes
fits <- fit_models(comp_cds, model_formula_str = "~gene_target", cores = 4)
mod_coef <- coefficient_table(fits)

# clean up results
celltype.degs <- mod_coef %>% 
  select(-model_summary, -model) %>% 
  filter(term != "(Intercept)") %>%
  mutate(up_in = case_when(
    estimate < 0 ~ "control",
    estimate > 0 ~ "tbxta"))

sig.celltype.degs <- celltype.degs %>% 
  filter(q_value < 0.01 ) %>%  
  arrange(up_in, q_value) %>%
  select(up_in, id, gene_short_name, q_value, 
         estimate)

# save
fwrite(sig.celltype.degs, "notochord_vs_PC_36hpf_DEGs.csv", 
       sep = ",", na = "NA")

## volcano plot

celltype.degs %>%
  mutate(sig = case_when(q_value < 0.01 ~ "yes",
                         TRUE ~ "no")) %>% 
  ggplot(aes(x= normalized_effect, y = -log10(q_value), color = sig)) +
  geom_point() +
  #gghighlight::gghighlight(gene_short_name == "epyc", label_key = gene_short_name) +
  #gghighlight::gghighlight(q_value < 1e-25, label_key = gene_short_name) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(file = "nc-v-pc_volcano_plot2.pdf", 
       width = 3, height = 3)

# EDF20b-e - plot gene expression for PC -------------------------------

# tidy genotypes and subset data
colData(notochord_cds)$geno_summary = colData(notochord_cds) %>% 
  as.data.frame() %>% 
  mutate(tmp = case_when(grepl("ctrl", gene_target) ~ "control",
                         TRUE ~ gene_target)) %>% 
  pull(tmp)

sub_cds <- notochord_cds[,colData(notochord_cds)$geno_summary %in% c("control", "tbxta")]

# filter for 36 hpf cells and downsample ctrls to compare better
unique(colData(sub_cds)$timepoint)
sub_cds <- notochord_cds[,colData(notochord_cds)$geno_summary %in% c("control", "tbxta")
                        & colData(notochord_cds)$timepoint == "36"]

coldata_df <- colData(sub_cds) %>% 
  as.data.frame()

coldata_df$cluster <- as.numeric(clusters(sub_cds))
coldata_df$cell_barcode <- rownames(coldata_df)

downsamp_df <- coldata_df %>% 
  filter(geno_summary == "control" &
           !(cluster %in% c(4, 8))) %>% 
  group_by(cluster) %>% 
  sample_frac(size = 0.14)

samp_wt_cells <- downsamp_df %>% 
  pull(cell_barcode)

tbxta_cells <- coldata_df %>% 
  filter(gene_target == "tbxta") %>% 
  pull(cell_barcode)

keep_cells <- c(tbxta_cells, samp_wt_cells)

# get the downsampled cds
samp_cds <- sub_cds[,keep_cells]
ncol(samp_cds)

# plot PC genes in controls vs. tbxta

pc_genes <- c("cnmd", "mvp", "matn4", "tgm2l", "epyc")

# plot gene expression
for (gene in pc_genes) {
  plot_cells(samp_cds, x = 1, y = 3,  genes = c(gene), 
             label_cell_groups = F, cell_size = 0.6) +
    scale_color_viridis_c() +
    theme_void() +
    theme(legend.position = "none", strip.text.x = element_blank()) +
    facet_wrap(~geno_summary)
  
  ggsave("EDF20_all-notochord_", gene, "-expr_facet_36hpf_downsamp.png",
         dpi = 750,
         height = 2,
         width = 5,
         bg = "transparent")
}

# get legend info
plot_cells(sub_cds, x = 1, y = 3,  genes = c("epyc"), 
           label_cell_groups = F, cell_size = 0.5) +
  scale_color_viridis_c() +
  theme_void()

