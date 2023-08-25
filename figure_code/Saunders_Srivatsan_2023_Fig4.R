## Saunders, Srivatsan, et al. (2023)
## This script contains code to generate plots for Figure 4
## from processed data files available in the GEO repository GSE202639 or available via Github.

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(viridis)
  library(data.table)
  library(devtools)
  library(tidyr)
  library(plotly)
  library(htmlwidgets)
  library(VGAM)
  library(forcats)
})

source("saunders_srivatsan_2023_utils.R")

# Load data ---------------------------------------------------------------

# sensory neurons from all genotypes
cg_cds <- readRDS("data/all-geno_sensory-cranial-ganglion_neuron_29k_cds.RDS")

# append 3D umap coordinates
colData(cg_cds)$sub_umap1 <- reducedDim(x = cg_cds,
                                        type = "UMAP")[,1]
colData(cg_cds)$sub_umap2 <- reducedDim(x = cg_cds,
                                        type = "UMAP")[,2]
colData(cg_cds)$sub_umap3 <- reducedDim(x = cg_cds,
                                        type = "UMAP")[,3]

colData(cg_cds)$pseudotime <- pseudotime(cg_cds)

# get coldata
cg_coldata <- colData(cg_cds) %>% 
  as.data.frame()

# get cell type colors
ganglia_colors = 
  c("cranial ganglion progenitor" = "#A29ADE", 
    "trigeminal ganglion" = "#B4286C", 
    "epibranchial ganglion" = "#2E74BC", 
    "lateral line ganglion" = "#5ABC6B",
    "statoacoustic ganglion" = "#5ED3F3",
    "unknown sensory ganglion" = "#ADADAD",
    "rohon-beard neuron" = "#C9AE56")


# Figure 4b  --------------------------

## reference inset

# load full reference metadata
global_coldata <- fread("data/final_final_ref/R_objects_ref/full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_coldata.csv", 
                       sep = ",", stringsAsFactors = F, data.table = F)

cg_cells <- cg_all_coldata %>% 
  pull(cell)

cg_cells <- intersect(cg_cells, global_coldata$cell)

ggplot() +
  geom_point(data = global_coldata, 
             aes(x = umap3d_1,
                 y = umap3d_2),
             color = "grey90",
             stroke = 0,
             size = 0.1) +
  geom_point(data = global_coldata %>% 
               filter(cell %in% cg_cells),
             aes(x = umap3d_1,
                 y = umap3d_2),
             color = "black",
             stroke = 0,
             size = 0.1) +
  theme_void() +
  theme(legend.position = "none")
ggsave("fig4b-1_ref_cells_ganglia-highlight_umap.png",
         dpi = 750,
         height = 1.5,
         width = 1.5,
         bg = "transparent")

## Cranial ganglia UMAP by timepoint

# get timepoint colors
rainbow_timepoint_colors <- 
  c("18h" = "#DF4828", 
    "24h" = "#E78C35", 
    "36h" = "#F6C141", 
    "48h" = "#4EB265", 
    "72h" = "#1965B0")

num_colors <- cg_coldata %>%
  pull(timepoint) %>%
  unique() %>%
  sort() %>%
  length()

full_spectrum_timepoint <- colorRampPalette(rainbow_timepoint_colors)(num_colors) 

names(full_spectrum_timepoint) <- cg_coldata %>%
  filter(!is.na(timepoint)) %>%
  pull(timepoint) %>%
  unique() %>%
  sort()

# plot timepoint in 2D
ggplot(cg_coldata %>%
         sample_n(size = dim(cg_coldata)[1])) +
  geom_point(aes(x = sub_umap1,
                 y = sub_umap3),
             color = "black",
             stroke = 0,
             size = 0.5) +
  geom_point(aes(x = sub_umap1,
                 y = sub_umap3,
                 color = timepoint),
             stroke = 0,
             size = 0.4) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = full_spectrum_timepoint) 
ggsave("fig4b-2_all_cranial-ganglia_timepoint_umap.png",
       dpi = 750,
       height = 2,
       width = 2.2,
       bg = "transparent")

# Figure 4c ---------------------------------------------

# get cell type colors
ganglia_colors = 
  c("cranial ganglion progenitor" = "#A29ADE", 
    "trigeminal ganglion" = "#B4286C", 
    "epibranchial ganglion" = "#2E74BC", 
    "lateral line ganglion" = "#5ABC6B",
    "statoacoustic ganglion" = "#5ED3F3",
    "unknown sensory ganglion" = "#ADADAD",
    "rohon-beard neuron" = "#C9AE56")

# plot in 2D
ggplot(cg_coldata %>%
         sample_n(size = dim(cg_coldata)[1])) +
  geom_point(aes(x = sub_umap1,
                 y = sub_umap3),
             color = "black",
             stroke = 0,
             size = 0.5) +
  geom_point(aes(x = sub_umap1,
                 y = sub_umap3,
                 color = cell_type_sub),
             stroke = 0,
             size = 0.4) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = ganglia_colors) 
ggsave("fig4c_all_cranial-ganglia_subtype_umap.png",
       dpi = 750,
       height = 2,
       width = 2.2,
       bg = "transparent")

# Figure 4d -----------------------------------------------

## Load TFs and model their expression along each trajectory arm ##

ref_cds <- readRDS("data/all-ref_sensory_neuron_cds.RDS")

# load genes to test
genes_to_test <- fread("data/cranial-ganglia_TF_top-markers.csv",
                      sep = ",", stringsAsFactors = F, data.table = F)[,1]
test_ids <- get.gene.ids(ref_cds, names = genes_to_test)

## Fit a linear model using pseudotime as the predictor 

# 1. Epibranchial pred
epi_cds <- ref_cds[test_ids, colData(ref_cds)$cell_type_sub %in% c("cranial ganglion progenitor", 
                                                                  "epibranchial")]
epi_trajectory_degs <- fit_models(epi_cds,
             model_formula_str = "~splines::ns(pseudotime, df=3)")

new_data_epi <- data.frame(
    Size_Factor = 1 ,
    pseudotime = seq(min(pseudotime(epi_cds)),
                     max(pseudotime(epi_cds)), 
                     length.out= 100))

# Use the model to predict the gene expression along this trajectory
predictions_epi <- model_predictions(model_tbl = epi_trajectory_degs, 
                    new_data=new_data_epi)

# 2. Statoacoustic pred
sa_cds <- ref_cds[test_ids, colData(ref_cds)$cell_type_sub %in% c("cranial ganglion progenitor", 
                                                                 "statoacoustic")]
sa_trajectory_degs <- fit_models(sa_cds,
             model_formula_str = "~splines::ns(pseudotime, df=3)")

new_data_sa <- data.frame(Size_Factor = 1,
    pseudotime = seq(min(pseudotime(epi_cds)),
                     max(pseudotime(epi_cds)), 
                     length.out= 100))

predictions_sa <- model_predictions(model_tbl = sa_trajectory_degs, 
                    new_data=new_data_sa)

# 3. Trigeminal pred
tg_cds <- ref_cds[test_ids, colData(ref_cds)$cell_type_sub %in% c("cranial ganglion progenitor", 
                                                                 "trigeminal")]
tg_trajectory_degs <- fit_models(tg_cds,
             model_formula_str = "~splines::ns(pseudotime, df=3)")

new_data_tg <- data.frame(
    Size_Factor = 1 ,
    pseudotime = seq(min(pseudotime(tg_cds)),
                     max(pseudotime(tg_cds)), 
                     length.out= 100))

predictions_tg <- model_predictions(model_tbl = tg_trajectory_degs, 
                    new_data=new_data_tg)

# 4. Lateral line pred
ll_cds <- ref_cds[test_ids, colData(ref_cds)$cell_type_sub %in% c("cranial ganglion progenitor", 
                                                                 "lateral line ganglion")]
ll_trajectory_degs <- fit_models(ll_cds,
             model_formula_str = "~splines::ns(pseudotime, df=3)")

new_data_ll <- data.frame(Size_Factor = 1 ,
    pseudotime = seq(min(pseudotime(ll_cds)),
                     max(pseudotime(ll_cds)), 
                     length.out= 100))

# Use the model to predict the gene expression along this trajectory
predictions_ll <- model_predictions(model_tbl = ll_trajectory_degs, 
                    new_data=new_data_ll)

## Now plot these genes in a pseudotime heatmap

row_center = function(m){
  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  return(m)
}

predictions_intersection = 
  cbind(row_center(predictions_epi),
        row_center(predictions_ll),
        row_center(predictions_sa),
        row_center(predictions_tg))

predictions_intersection[predictions_intersection > 3] <- 3
predictions_intersection[predictions_intersection < -1] <- -1

bks <- seq(-1,3, by = 0.1)
hmcols <- viridis::viridis(length(bks) - 1,option = "D")

col_gaps_ind = c(100, 200, 300)

rownames(predictions_intersection) = rowData(ref_cds)[rownames(predictions_intersection),"gene_short_name"]
rowSums(is.na(predictions_intersection)) %>% sort()

col.annotations = 
  data.frame(row.names = cutree(ph$tree_row,k = 4) %>% names(),
             cut = cutree(ph$tree_row,k = 4) %>% as.character())
annotation_colors = 
  list(cut = c("1" = RColorBrewer::brewer.pal(n = 4,name = "Set1")[1],
               "2" = RColorBrewer::brewer.pal(n = 4,name = "Set1")[2],
               "4" = RColorBrewer::brewer.pal(n = 4,name = "Set1")[4],
               "3" = RColorBrewer::brewer.pal(n = 4,name = "Set1")[3]))                               

# save plots
png("fig4d_pseudotime_heatmap.png", width = 4, height = 3, units = "in", res = 700)
pheatmap(predictions_intersection,
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = col_gaps_ind, 
         border_color = NULL,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         show_colnames = F,
         breaks=bks,
         legend = F,
         color=hmcols, width = 6, height = 5)
dev.off()

png("fig4d_pseudotime_heatmap_w-legend.png", width = 4, height = 3, units = "in", res = 700)
pheatmap(predictions_intersection,
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = col_gaps_ind,
         border_color = NULL,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         show_colnames = F,
         breaks=bks,
         legend = T,
         color=hmcols)
dev.off()

png("fig4d_pseudotime_heatmap_rownames.png", width = 4, height = 6, units = "in", res = 700)
pheatmap(predictions_intersection,
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = col_gaps_ind,
         border_color = NULL,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = T,
         show_colnames = F,
         breaks=bks,
         legend = F,
         color=hmcols)
dev.off()

# Figure 4e and EDF 18 --------------------------------------------------

# list in-situ hybridization genes
ish_genes <- c("syt9b", "kcnq2b", "ndnfl", "cpne7", "cpne4a", 
               "tspan2b", "scn5lab", "hs6st3a", "ntrk2a", 
               "irx4b", "irx5a", "p2rx3b")

for (g in ish_genes){
plot_cells(cg_cds, x = 1, y = 3,  genes = g, 
           label_cell_groups = F, show_trajectory_graph = F,
           cell_size = 0.5, cell_stroke = 0, alpha = 0.8) +
  scale_color_viridis_c() +
  theme_void() +
  theme(legend.position = "right")

  ggsave(paste0("fig4e_all-ganglia_", g, "-expr_umap.png"),
       dpi = 750,
       height = 2.5,
       width = 3,
       bg = "transparent")
}

# Figure 4f ----------------------

# load data
summary_merge <- fread("data/cranial-ganglia_all-geno_summary.csv",
                       sep = ",", data.table = F, stringsAsFactors = F, na.strings = "")

# plot by celltype facet
types <- c("cranial ganglion progenitor", "epibranchial ganglion", "trigeminal ganglion", "statoacoustic ganglion", "lateral line ganglion")
genos <- c("ctrl-inj", "phox2a", "foxi1", "tfap2a-foxd3")
times <- c("48")

# separate by genotype
cell_cols <- c('#BBBBBB', "#4477AA", '#3BD4F5', "#AA3377")

plot_df <- summary_merge %>% 
  filter(genotype %in% genos) %>% 
  filter(cell_type %in% types) %>%
  filter(timepoint %in% times)

# plot boxplots of cell counts for each ganglion across selected genotypes
ggplot(plot_df, aes(x = factor(genotype, levels = c("ctrl-inj", "phox2a", "foxi1", "tfap2a-foxd3")), 
                    y = cells, 
                    fill = cell_type)) +
  geom_boxplot(color = "black", size = 0.3, outlier.size = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.5, color = "black", width = .2) +
  theme(axis.title.x = element_blank(), legend.position = "none", 
        axis.ticks = element_line(color = "black", size = .3),
        legend.title = element_blank(), 
        axis.text.x = element_text(angle= 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black")) +
  #facet_wrap(~genotype, nrow = 1) +
  scale_fill_manual(values = ganglia_colors) +
  facet_wrap(~factor(cell_type, levels = types), 
             nrow = 1, scales = "free_y") +
  monocle3:::monocle_theme_opts()
ggsave("fig4f_ctrl-phox2a-foxi1-tfap-foxd3_cg-all_counts-boxplot_geno-facet.pdf", width = 7, height = 3)

# EDF 18 Subtype abundance heatmap -----------------------------------------------

# just genos for this fig
res_df <- fread("/Volumes/GoogleDrive/Shared drives/Trapnell Lab/Projects/SDG/GAPFISH/data/gap16/analysis/cranial_ganglia/cranial-ganglia_bb_res_list_36-48-72hpf_foxi1-phox2a-tfap-foxd3_df.csv", 
               sep = ",", data.table = F, stringsAsFactors = F, na.strings = "")

# filter results
qval_thresh <- 0.05
sig_df <- res_df %>%
  dplyr::mutate(sig_fc = case_when(qval < qval_thresh ~ TRUE, 
                                   TRUE ~ FALSE))

# filter on p-val
pval_thresh = 0.05
sig_df = res_df %>%
  filter(geno.v.wt %in% c("phox2a", "foxi1")) %>% 
  dplyr::mutate(sig_fc = case_when(p.value < pval_thresh ~ TRUE, 
                                   TRUE ~ FALSE))

# summarize fc in a matrix
hm_wide = sig_df %>%
  select(geno.v.wt, timepoint, cell_group, estimate) %>%
  pivot_wider(
    names_from = cell_group,
    values_from = estimate,
    values_fill = c(0)
  )

hm_plot = hm_wide %>%
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

new_levels = hm_plot %>% distinct(geno_timepoint) %>% pull %>% as.character %>% sort
hm_plot$geno_timepoint = factor(hm_plot$geno_timepoint, levels = new_levels)

# get clustering results
hm_wider = hm_wide %>%
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
fms = hm_plot[hm_plot$sig_fc, c("geno_timepoint", "cell_group")]
fms$geno_timepoint = as.integer(fms$geno_timepoint)
fms$cell_group = as.integer(fms$cell_group)

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
    axis.line = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 0.25)
  ) +
  coord_equal() +
  monocle3:::monocle_theme_opts()

ggsave("EDF18_phox2a-foxi1-tfap2-foxd3_est-heatmap_p05.pdf",
       units = "in",
       width = 6,
       height = 4)
