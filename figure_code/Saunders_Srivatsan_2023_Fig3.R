## Saunders, Srivatsan, et al. Nature (2023) - Embryo-scale reverse genetics at single-cell resolution
## This script contains code to generate plots for Figure 3 and EDF 15 
## from processed data files available in this github repository ("data" folder), at https://cole-trapnell-lab.github.io/zscape/
## or in the GEO repository GSE202639.

# startup ---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(stringr)
  library(pheatmap)
  library(viridisLite)
  library(data.table)
  library(gghighlight)
  library(gprofiler2)
})

source("saunders_srivatsan_2023_utils.R")

# Figure 3a ------------------------------------------------

# load DEG results
deg_res <- fread("data/Saunders_Srivatsan_DEG-list_q05filt.csv", 
                 sep = ",", data.table = F, 
                 stringsAsFactors = F)

# load DEG summary based on a q-value cutoff of 0.05
deg_summary <- fread("data/all_broad_type_pseudobulked_DEG_q05_summary.csv", 
                    sep = ",", data.table = F, 
                    stringsAsFactors = F)

# per type DEGs
cns.results.mat <- deg_res %>%
  filter(!gene_target %in% c("noto-mut")) %>%
  group_by(gene_target,
           cell.type) %>%
  summarise(n.degs = n()) %>%
  left_join(deg_summary %>% 
              dplyr::rename(deg_total = "n.degs"), 
            by = "cell.type", relationship = "many-to-many") %>%
  dplyr::filter(tissue == "Central Nervous System") %>% 
  #dplyr::filter(deg_total > 50) %>%
  distinct() %>%
  dplyr::select(-deg_total, -tissue) %>% 
  tidyr::pivot_wider(names_from = cell.type, 
               values_from = n.degs, 
               values_fill = 0) %>%
  as.data.frame()
rownames(cns.results.mat) <- cns.results.mat$gene_target
cns.results.mat <- cns.results.mat %>% 
  dplyr::select(-gene_target)
cns.results.mat <- as.matrix(cns.results.mat)

# plot heatmap with log scale
pheatmap(log10(cns.results.mat+1),
         color = colorRampPalette(colors = c("white","white","#f08080", "#ce4257", "#720026", "#4f000b"))(100),
         clustering_method = "ward.D2",
         cluster_cols = F,
         border_color = NULL,
         filename = "Fig3a_cns.DEG.heatmap.q05.log10.pdf",
         cellheight = 8,
         cellwidth = 8)


# Figure 3b ---------------------------------------------------------------

# load abundance data
dact_df <- fread("data/zperturb_beta-bin_results_broadtype_clean.csv", 
                 sep = ",", data.table = F) # DACT results

# get mutant x cell type DEG summary
deg_perturb_ct_df <- deg_res %>% 
  dplyr::group_by(gene_target, cell.type) %>% 
  dplyr::mutate(n.degs = sum(n())) %>% 
  select(gene_target, cell.type, n.degs, tissue) %>% 
  distinct()

# filter results
qval_thresh <- 0.01
hb_targets <- c("egr2b", "cdx4", "cdx4-cdx1a", "hoxb1a", "mafba", "epha4a", "wnt3a-wnt8", "smo")
timept <- 24 

cns_abund_dea_df <- dact_df %>% 
  filter(timepoint == timept) %>%
  mutate(cell.type = gsub('\\/','-', cell_group)) %>%
  select(gene_target = geno.v.wt, cell.type, abund_qval = qval, abund_log2fc) %>% 
  left_join(deg_perturb_ct_df, by = c("gene_target", "cell.type")) %>%
  filter(tissue == "Central Nervous System") %>%
  mutate(n.degs = replace_na(n.degs, replace = 0)) %>%
  mutate(abund_dir = case_when(abund_log2fc >= 0 ~ "up",
                               abund_log2fc < 0 ~ "down"))

# plot abundance vs DEG scatter with and without legend
cns_abund_dea_df %>% 
  filter(gene_target %in% hb_targets) %>% 
  ggplot(aes(x = abs(abund_log2fc), y = n.degs, color = abund_dir)) +
  geom_point(size = 0.75) +
  scale_color_manual(values = c("tomato", "steelblue")) +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(size = 8))

ggsave("Fig3b_cns_deg_vs_abundFC_scatter_fig3.pdf", width = 2.3, height = 2.3, units = "in")
# ggsave("Fig3b_cns_deg_vs_abundFC_scatter_w-legend.pdf", width = 3, height = 3, units = "in", dpi = 700) # legend.position = "right"


# Figure 3c  ----------------------------------

# load HB DEG data
hb_dea_res_df <- fread("data/neural-progenitor_all-genos_deg-heatmap_df.csv", 
                       sep = ",", data.table = F, 
                       stringsAsFactors = F)

# select tissues and celltypes
hb_targets <- c("egr2b", "cdx4", "cdx4-cdx1a", "hoxb1a", "mafba", "epha4a", "ctrl-inj", "wnt3a-wnt8", "smo")
hb_types <- c("neural progenitor (hindbrain)", "neural progenitor (hindbrain R7/8)", 
              "differentiating neuron (hindbrain)", 
              "neural progenitor (telencephalon/diencephalon)")

# make matrix of the DEG test betas
hb_beta_wide <- hb_dea_res_df %>% 
  filter(cell.type == "neural progenitor (hindbrain)") %>% 
  filter(mutant %in% hb_targets) %>% 
  filter(id %in% deg_ids) %>%
  distinct() %>%
  select(id, normalized_effect, mutant) %>% 
  pivot_wider(
    names_from = mutant,
    values_from = normalized_effect,
    values_fill = c(0))

# format matrix and set bounds
hb_beta_mat <- as.matrix(hb_beta_wide[,-1])
rownames(hb_beta_mat) <- hb_beta_wide$id
hb_beta_mat[hb_beta_mat > 3] <- 3
hb_beta_mat[hb_beta_mat < -3] <- -3

# set color palette
BuRd <- c('#2166AC', '#4393C3', '#92C5DE', '#D1E5F0', '#F7F7F7', "#fad2d2", "#f08080", "#ce4257", "#720026")

# no cluster breaks
pdf("fig3c_hindbrain-np_8-geno_no-clust_expr-heatmap.pdf", width = 3, height = 4, bg = "transparent", useDingbats = F)           
pheatmap(hb_beta_mat,
         color = colorRampPalette(BuRd)(100),
         cluster_rows = T,
         cluster_cols = T, border_color = "NA",
         show_rownames = F, show_colnames = T, clustering_method = "ward.D2")
dev.off()

# plot heatmap with clusters
pdf("fig3c_hindbrain-np_8-geno_expr-heatmap.pdf", width = 3, height = 4, bg = "transparent", useDingbats = F)           
pheatmap(hb_beta_mat,
         color = colorRampPalette(BuRd)(100),
         cluster_rows = T,
         cluster_cols = T, border_color = "NA", cutree_rows = 5,
         show_rownames = F, show_colnames = T, clustering_method = "ward.D2")
dev.off()

# save with gene names
pdf("fig3c_hindbrain-np_8-geno_expr-heatmap_names.pdf", width = 5, height = 22, bg = "transparent", useDingbats = F)           
pheatmap(hb_beta_mat,
         color = colorRampPalette(BuRd)(100),
         cluster_rows = T,
         cluster_cols = T, border_color = "NA", cutree_rows = 5,
         show_rownames = T, show_colnames = T, clustering_method = "ward.D2", fontsize = 4)
dev.off()

### Gene ontology analysis for HM clusters and by target group ###

# 1. heatmap clusters
# extract genes from heatmap (above) clusters
hm_out <- pheatmap(hb_beta_mat,
                   color = colorRampPalette(BuRd)(100),
                   cluster_rows = T,
                   cluster_cols = T, border_color = "NA", cutree_rows = 5,
                   show_rownames = F, show_colnames = T, clustering_method = "ward.D2")

clust_df = data.frame(cluster = sort(cutree(hm_out$tree_row, k=5))) %>% 
  tibble::rownames_to_column("id")

# which genes are in which cluster?
clust_df %>% 
  filter(cluster == 2) %>% 
  head()

# format for GO testing
mods_for_go = tapply(clust_df$id, clust_df$cluster, list)

go_res = gost(mods_for_go, 
              organism = "drerio",
              multi_query = FALSE, evcodes = T) #, custom_bg = hb_expr_genes

go_res_df = go_res$result
go_res_df = go_res_df %>% 
  select(deg_clust = query, term_id, term_name, 
         p_value, term_size, intersection_size)

dim(go_res_df)

go_res_df %>% 
  filter(deg_clust == 4) %>% 
  head(10)

# 2. By Target group
grp_1 = c("egr2b", "mafba")
grp_2 = c("cdx4", "cdx4-cdx1a", "wnt3a-wnt8", "smo")

# import and filter HB DEGs
degs_by_grp <- hb_dea_res_df %>% 
  filter(cell.type == "neural progenitor (hindbrain)") %>% 
  filter(q_value < 0.05) %>% 
  mutate(grp = case_when(mutant %in% grp_1 ~ "egr2_grp",
                         mutant %in% grp_2 ~ "cdx_grp",
                         TRUE ~ NA_character_)) %>% 
  filter(!is.na(grp)) %>% 
  select(grp, id)

# format for GO testing
mods_for_go <- tapply(degs_by_grp$id, degs_by_grp$grp, list)

go_res = gost(mods_for_go, 
              organism = "drerio",
              multi_query = FALSE, evcodes = T)
go_res_df_grp = go_res$result
go_res_df_grp = go_res_df_grp %>% 
  select(grp = query, term_id, term_name, 
         p_value, term_size, intersection_size)

go_res_df_grp %>% 
  filter(grp == "cdx_grp") %>% 
  head(10)

# Figure 3d ---------------------------------------------------------------

# load data
hb_cds <- readRDS("data/hindbrain-targets_hb-mn-cells_18-36h_60k_cds.RDS")

# add umap coords to coldata
hb_coldata <- colData(hb_cds) %>% 
  as.data.frame() 

umap_coords <- as.data.frame(reducedDim(hb_cds, "UMAP"))
colnames(umap_coords) = c("subumap_1", "subumap_2")
umap_coords$cell = rownames(umap_coords)

hb_coldata <- hb_coldata %>% 
  left_join(umap_coords, by = "cell")

# assign colors for celltypes
my.colors = c("neural progenitor (hindbrain)" = "#66CCEE",
              "neural progenitor (hindbrain R7/8)" = "#EE6677",
              "differentiating neuron (hindbrain)" = "#4477AA",
              "neural progenitor (telencephalon/diencephalon)" = "#56A85C")

ggplot(hb_coldata %>% 
         filter(cell_type_broad != "motor neuron") %>% 
         filter(subumap_1 < 10)) +
  geom_point(aes(x = subumap_1,
                 y = subumap_2),
             color = "black",
             stroke = 0,
             size = 0.6) +
  geom_point(aes(x = subumap_1,
                 y = subumap_2,
                 color = cell_type_broad),
             stroke = 0,
             size = 0.5,
             alpha = 0.8) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = my.colors) 

ggsave("Fig3d_hindbrain_hb-targets_all-timepoints_celltypebroad_umap.png",
       dpi = 750,
       height = 2,
       width = 2,
       bg = "transparent")

# Figure 3e ----------------------------------------------

## load data ##
# gene expression
hb_cds <- readRDS("data/hindbrain-targets2_hb-mn-cells_18-36h_66k_cds.RDS")

# coldspot analysis results
all_bind_hs_df <- fread("data/hindbrain-subset_all-types-together_coldspot-results_220121.csv", 
                        stringsAsFactors = F, sep = ",", data.table = F)

coldspot_summary <- all_bind_hs_df %>% 
  filter(gene_target == "ctrl-inj") %>%
  mutate(sig = adjusted_pval < 0.05) %>% 
  filter(sig) %>% 
  group_by(comparison, timepoint, cell_type_broad) %>% 
  tally()

# plot areas where mutant cells are depleted relative to controls

targets_only <- unique(all_bind_hs_df$gene_target)[unique(all_bind_hs_df$gene_target) != "ctrl-inj"]

for (gene in targets_only){
  
  filt_df = all_bind_hs_df %>% 
    filter(gene_target == "ctrl-inj") %>% 
    filter(comparison == gene)
  LG_thresh = mean(filt_df$localG + 1*sd(filt_df$localG))
  
  plot_df = all_bind_hs_df %>% 
    filter(cell_type_broad != "motor neuron") %>% 
    filter(gene_target == "ctrl-inj") %>%
    filter(comparison == gene) %>% 
    mutate(hs_call = localG > LG_thresh)
  
  ggplot() +
    geom_point(data = plot_df %>%
                 arrange(hs_call),
               aes(subumap_1, subumap_2), color = "gray90", alpha = 1, size = 0.5, stroke = 0) +
    geom_point(data = plot_df %>% 
                 filter(hs_call),
               aes(subumap_1, subumap_2), color = "#0466c8", alpha = 0.6, size = 0.6, stroke = 0) +
    theme_void() +
    theme(legend.position = "none")
  ggsave(filename = paste0("fig3e_hindbrain-coldspot_", gene, "_combined-tmpts_simple_umap.png"), 
         dpi = 750, width = 2, height = 2, bg = "transparent", units = "in")
  
}

# Figure 3f ---------------------------------

# plot epha4a expression in controls and egr2b-crispant hindbrain cells
my_gene <- c("epha4a")
my_genos <- c("ctrl-inj", "egr2b")

plot_cells(hb_cds[,colData(hb_cds)$gene_target %in% my_genos & !(colData(hb_cds)$cell_type_broad == "motor neuron")], 
           genes = my_gene, label_cell_groups = F, cell_size = 0.4, alpha = 0.8) +
  scale_color_viridis_c() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~gene_target)

ggsave(paste0("fig3f_hindbrain-cds_ctrl_hb-genes_", my_gene, "_DEG_expr_umap-facet_legend.png"),
       dpi = 750,
       height = 2,
       width = 4,
       bg = "transparent")


# Figure 3g ---------------------------------------------------------------

# plot hox gene expression between ctrl and cdx4;cdx1a-crispant HB cells
my_genes <- c("hoxb3a", "hoxc3a", "hoxc6b")
my_genos <- c("ctrl-inj", "cdx4-cdx1a")

for (g in my_genes){
plot_cells(hb_cds[,colData(hb_cds)$gene_target %in% my_genos & !(colData(hb_cds)$cell_type_broad == "motor neuron")], 
           genes = g, label_cell_groups = F, cell_size = 0.4, alpha = 0.8) +
  scale_color_viridis_c() +
  theme_void() +
  facet_wrap(~forcats::fct_rev(gene_target), dir = "v") +
  theme(legend.position = "bottom", strip.text.x = element_blank())
ggsave(paste0("hindbrain-cds_ctrl_hb-genes_", g, "_cdx4-cdx1-DEG_expr_umap-facet_vert_legend_fig3.png"),
       dpi = 750,
       height = 4,
       width = 2,
       bg = "transparent")
}

# Extended data figure 15 -------------------------------------------------

# load DEG summary based on a q-value cutoff of 0.05
deg_summary <- fread("data/broad_type_pseudobulked_DEG_summary.csv", 
                     sep = ",", data.table = F, 
                     stringsAsFactors = F)

# summarize deg results into matrix
results.mat =
  deg_summary %>%
  filter(tissue %in% tissue_sel) %>% 
  pivot_wider(key = cell.type,
              value = (n),
              fill = 0) %>%
  as.data.frame()
rownames(results.mat) = results.mat$mutant
results.mat = results.mat %>% 
  dplyr::select(-mutant)
results.mat = as.matrix(results.mat)

# group matrix by tissue
ord.df = results %>% 
  select(tissue, cell.type) %>% 
  distinct() %>% 
  arrange(tissue)
type.order = ord.df$cell.type

order.int = intersect(type.order, colnames(results.mat))
ord.results.mat = results.mat[,order.int]

#plot and save (this is for supp fig 15)
pheatmap(log10(ord.results.mat+1),
         color = colorRampPalette(colors = c("white", "white","#f08080", "#ce4257", "#720026", "#4f000b"))(100),
         clustering_method = "ward.D2",
         border_color = NA,
         filename = "DEG.heatmap.q05.emb-filt.log10.tissue-order.pdf",
         cluster_rows = T, cluster_cols = F,
         cellheight = 8,
         cellwidth = 8, angle_col = 90)
