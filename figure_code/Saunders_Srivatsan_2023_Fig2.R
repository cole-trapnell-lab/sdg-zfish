## Saunders, Srivatsan, et al. Nature (2023) - Embryo-scale reverse genetics at single-cell resolution
## This script contains code to generate plots for Figure 2 
## from processed data files available in this github repository ("data" folder), at https://cole-trapnell-lab.github.io/zscape/
## or in the GEO repository GSE202639.
## Included here is also code for Extended data figures 12 and 14 

# startup ---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(plotly)
  library(viridis)
  library(data.table)
  library(DelayedArray)
  library(devtools)
  library(monocle3)
  library(htmlwidgets)
  library(VGAM)
  library(forcats)
 
})

source("saunders_srivatsan_2023_utils.R")

# Load data -------------------------------------------------------

ref_emb_cds <- readRDS("data/all-ref_cell_count_cds.RDS") # reference emb cell data set
mut_emb_cds <- readRDS("data/zperturb_cell_count_cds.RDS") # zperturb emb cds
emb_cds_36h <- readRDS("data/mut_emb_cds_36hpf_with_coords.RDS") # embryo cds
dact_df <- fread("data/zperturb_beta-bin_results_broadtype_clean.csv", 
                 sep = ",", data.table = F) # DACT results

# Fig 2c -------------------------------------------

# summarize DACT number per perturbation
dact_summary <- dact_df %>% 
  select(gene_target = geno.v.wt, cell_group, timepoint, qval) %>% 
  filter(qval < 0.05) %>% 
  distinct() %>% 
  group_by(gene_target, timepoint) %>% 
  tally(name = "dact_n")

dact_summary$timepoint <- as.character(dact_summary$timepoint)

# add control label
colData(emb_cds_36h)$is_ctrl = colData(emb_cds_36h) %>% 
  as.data.frame() %>% 
  mutate(is_ctrl = case_when(grepl("ctrl", gene_target) ~ TRUE,
                             TRUE ~ FALSE)) %>% 
  pull(is_ctrl)

# extract data for plotting
coldat_filt <- colData(emb_cds_36h) %>% 
  as_tibble() %>%
  left_join(dact_summary, by = c("gene_target", "timepoint")) %>% 
  distinct()

# plot embryos colored by genotype and sized by number of dacts
mut_colors = c("ctrl-inj" = "grey80", 
               "cdx4" = "#F12CBA", 
               "cdx4-cdx1a" = "#A63587", 
               "tbxta" = "#F47534", 
               "wnt3a-wnt8" = "#E4A341",
               "tbx16" = "#A12926", 
               "tbx16-mut" = "#D63732", 
               "tbx16-tbx16l" = "#D63732", 
               "tbx16-msgn1" = "#E68784", 
               "smo" = "#510027", 
               "noto" = "#892E58", 
               "noto-mut" = "#892E58", 
               "egr2b" = "#056596", 
               "epha4a" = "#2D358D",
               "mafba" = "#29686D", 
               "mafba-mut" = "#29686D", 
               "hoxb1a" = "#785BDA", 
               "zc4h2" = "#50849A",
               "tbx1" = "#5640F1", 
               "hand2" = "#DA99B8", 
               "phox2a" = "#AC8BF0", 
               "met" = "#33B8D8",
               "met-mut" = "#33B8D8", 
               "foxi1" = "#2E91FF",
               "hgfa" = "#F5BE66", 
               "hgfa-mut" = "#F5BE66", 
               "tfap2a" = "#4AB633", 
               "foxd3" = "#3F7933", 
               "tfap2a-foxd3" = "#157B00")

ggplot() +
  geom_point(data = coldat_filt %>%
               filter(gene_target_test == "ctrl"),
             mapping = aes(x = umap1, y = umap2), 
             color = "grey70", size = 2, stroke = 0) +
  geom_point(data = coldat_filt %>%
               filter(gene_target_test != "ctrl"),
             mapping = aes(x = umap1, y = umap2,
                           color = gene_target_test, size = dact_n, alpha = 0.5, stroke = 0)) +
  scale_color_manual(values = mut_colors) +
  scale_size_binned(range = c(1, 4)) +
  theme_void() #+ 
  #theme(legend.position = "none")

# ggsave("fig2c_emb_36hpf_umap_DACT-size.pdf", width = 4, height = 4)
ggsave("emb_36hpf_umap_DACT-size_w-legend.pdf", width = 4, height = 4)

# Plot control embryos only
ggplot() +
  geom_point(data = coldat_filt %>% 
               filter(is_ctrl == T),
             mapping = aes(x = umap1, y = umap2,
                           label = gene_target), color = "grey70", 
             size = 1, stroke = 0, alpha = 0.6) +
  theme_void() +
  theme(legend.position = "none")

ggsave("fig2c_ctrl-only_36h_emb_umap.pdf", 
       width = 0.6, height = 0.6)


# Figure 2d ---------------------------------------------------------------

sig_df = fread("data/sig_broadtype_count_per-genotype-timepoint_q05_fc-pt5.csv", sep = ",", stringsAsFactors = F, 
                 data.table = F)

geno_order = c("cdx4",
                 "cdx4-cdx1a",
                 "tbxta",
                 "wnt3a-wnt8",
                 "tbx16",
                 "tbx16-msgn1",
                 "tbx16-tbx16l",
                 "smo",
                 "noto",
                 "egr2b",
                 "epha4a",
                 "mafba",
                 "hoxb1a",
                 "zc4h2",
                 "tbx1",
                 "hand2",
                 "phox2a",
                 "tfap2a",
                 "foxd3",
                 "tfap2a-foxd3",
                 "met",
                 "hgfa-mut",
                 "foxi1")

num_ct_df = dact_summary %>% 
  ungroup() %>% 
  filter(gene_target %in% geno_order) %>%
  mutate(gene_target = factor(gene_target, levels = geno_order),
         timepoint = factor(timepoint))

ggplot(num_ct_df, aes(x = timepoint, 
                      y = fct_rev(gene_target), 
                      fill = dact_n)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  scale_fill_viridis_c(option = "G") +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_blank(), 
        text = element_text(size = 8, color = "black"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", size = 1), 
        axis.ticks = element_line(colour = "black", size = 0.25)) +
  monocle3:::monocle_theme_opts()

ggsave("geno_timepoint_num-types_dotplot.pdf", width = 2.3, height = 4)
  
# Figure 2f ---------------------------------------------------------------

# tbx16 group boxplots

# load data
summary_merge = fread("data/zperturb-all_broadtype_summary_merge.csv",
                      sep = ",", data.table = F, stringsAsFactors = F, na.strings = "")

# plot by cell type facet
genos <- c("ctrl-inj", "tbx16", "tbx16-msgn1", "tbx16-tbx16l")
types <- c("posterior spinal cord progenitors",
          "fast-committed myocyte",
          "mesodermal progenitor cells (contains PSM)")
times <- c("24")

cell_cols <- c("grey80", "#ffd5c2", "#f28f3b", "#c8553d")

summary_merge %>% 
  filter(genotype %in% genos) %>% 
  filter(cell_type %in% types) %>%
  filter(timepoint %in% times) %>%
  mutate(cells_per_1000 = round(cells / total_cells * 1000)) %>%
  ggplot(aes(x = factor(genotype, levels = genos), 
             y = cells_per_1000, 
             fill = genotype)) +
  geom_boxplot(color = "black", size = 0.3, outlier.size = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.5, color = "black", width = .2) +
  theme(axis.title.x = element_blank(), legend.position = "none", 
        axis.ticks = element_line(color = "black", size = .3),
        legend.title = element_blank(), 
        axis.text.x = element_text(angle= 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black")) +
  scale_fill_manual(values = cell_cols) +
  facet_wrap(~cell_type, nrow = 1, scales = "free_y") +
  monocle3:::monocle_theme_opts()

ggsave("fig2f_tbx16-grp_counts-boxplot_CPT_celltype-facet.pdf", width = 6, height = 3)

# Related extended data figures -------------------------------------------

# Extended data figure 12 -------------------------------------------------

# plot embryo UMAP with all time points
emb_cds = readRDS("zperturb_cell_count_cds.RDS")

timept_clrs = 
  c("18" = "#DF4828", 
    "24" = "#E78C35", 
    "36" = "#F6C141", 
    "48" = "#4EB265", 
    "72" = "#1965B0")

# plot
emb_coldata <- colData(ref_emb_cds) %>% 
  as.data.frame() 

ggplot() +
  geom_point(data = emb_coldata %>% 
               filter(!grepl("ctrl", gene_target)),
             aes(x = umap1,
                 y = umap2,
                 color = timepoint),
             stroke = 0,
             size = .8,
             inherit.aes = FALSE) +
  geom_point(data = emb_coldata %>%
               filter(grepl("ctrl", gene_target)),
             aes(x = umap1,
                 y = umap2,
                 fill = timepoint),
             shape = 21,
             color = "black",
             size = .8,
             stroke = 0.5,
             inherit.aes = FALSE) +
  scale_color_manual(values = timept_clrs) +
  scale_fill_manual(values = timept_clrs) +
  theme_void() +
  theme(legend.position = "none")


ggsave("EDF12_gap16-all_embryo_wt-highlight_timept_umap.pdf", 
       width = 3, height = 3,
       bg = "transparent")

# Extended data figure 14 -------------------------------------------------

res_df <- fread("data/zperturb_beta-bin_results_broadtype_clean.csv", 
               sep = ",", data.table = F, stringsAsFactors = F, na.strings = "")

# filter for cell type and genotypes of interest
genos = c("tbx16", "tbx16-msgn1", "tbx16-tbx16l", "ctrl-inj", "tbx16-mut")
celltypes = c("posterior spinal cord progenitors",
              "mesodermal progenitor cells (contains PSM)",
              "dorsal spinal cord neuron",
              "mature slow muscle", 
              "mature fast muscle",
              "fast-committed myocyte", "slow-committed myocyte",
              "dorsal spinal cord neuron",
              "myoblast",
              "neuron (+ spinal cord)",
              "red blood cell")

qval_thresh = 0.01
sig_df = res_df %>% 
  filter(geno.v.wt %in% genos) %>% 
  filter(cell_group %in% celltypes) %>% 
  dplyr::mutate(sig_fc = case_when(qval < qval_thresh ~ TRUE, 
                                   TRUE ~ FALSE))

# summarize fc in a matrix
hm_wide = sig_df %>%
  select(geno.v.wt, timepoint, cell_group, abund_log2fc) %>%
  pivot_wider(
    names_from = cell_group,
    values_from = abund_log2fc,
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

new_levels = hm_plot %>% distinct(geno_timepoint) %>% pull %>% as.character
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
       aes(y = cell_group, x = geno_timepoint, fill = plot)) +
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
    size = 0.3,
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
  coord_flip() +
  monocle3:::monocle_theme_opts()

ggsave("EDF14_tbx16-grp_sc-neuron_meso_abund-heatmap.pdf",
       units = "in",
       width = 5,
       height = 5)

# plot a representation of mean cells per embryo for each cell type above
# use control mean values

ctrl_coldata <- fread("reference_cell_metadata.csv") # download from ZSCAPE

coldat_filt <- ctrl_coldata %>% 
  filter(gene_target == "ctrl-inj") %>%
  filter(timepoint == 24) %>% 
  group_by(embryo, cell_type_broad) %>% 
  summarize(ct_tot = sum(n())) %>%
  ungroup()

emb_tot_mean <- ctrl_coldata %>% 
  filter(gene_target == "ctrl-inj") %>% 
  filter(timepoint == 24) %>% 
  group_by(embryo) %>% 
  tally(name = "emb_tot")

coldat_filt <- coldat_filt %>% 
  left_join(emb_tot_mean, 
            by = "embryo") %>% 
  mutate(ct_frac = ct_tot/emb_tot)

mean_fracs_df <- coldat_filt %>% 
  group_by(cell_type_broad) %>% 
  summarize(mean_frac = round(100*mean(ct_frac), 2))

mean_fracs_df %>% 
  arrange(-mean_frac) %>% 
  head(20)

sum(mean_fracs_df$mean_frac)

# plot
ggplot(mean_fracs_df %>% 
         filter(cell_type_broad %in% celltypes),
       aes(y = 1, x = cell_type_broad, fill = mean_frac)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("#fdf8e1","#ffea00", "#ffb700", "#ff9500", "#ff8800")) +
  coord_equal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank()) +
  monocle3:::monocle_theme_opts()

ggsave("EDF14_tbx16-grp_sc-neuron_meso_abund-heatmap_mean-percent.pdf",
       units = "in",
       width = 5,
       height = 5)

