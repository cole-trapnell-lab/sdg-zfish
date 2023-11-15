## Saunders, Srivatsan, et al. Nature (2023) - Embryo-scale reverse genetics at single-cell resolution
## This script contains code to generate plots for Figure 1 
## from processed data files available in this github repository ("data" folder), at https://cole-trapnell-lab.github.io/zscape/
## or in the GEO repository GSE202639.
## Included here is also code for Extended data figure 6

# Startup -------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)
})

source("saunders_srivatsan_2023_utils.R")
 
# Load data --------------------------------------------------------

ref_cds <- readRDS("reference_cds.RDS") # download from ZSCAPE

ref_coldata <-  
  ref_cds %>%
  colData() %>%
  as.data.frame()

# Define colors for time points -------------------------------------------

rainbow_timepoint_colors <-  
  c("18h" = "#DF4828", 
    "24h" = "#E78C35", 
    "36h" = "#F6C141", 
    "48h" = "#4EB265", 
    "72h" = "#1965B0")

num_colors <- ref_coldata %>%
  filter(!is.na(timepoint)) %>%
  pull(timepoint) %>%
  unique() %>%
  sort() %>%
  length()

full_spectrum_timepoint <- colorRampPalette(rainbow_timepoint_colors)(num_colors) 

names(full_spectrum_timepoint) <- ref_coldata %>%
  filter(!is.na(timepoint)) %>%
  pull(timepoint) %>%
  unique() %>%
  sort()


# Make plot for time-series -----------------------------------------------

# Figure 1a ---------------------------------------------------------------

ggplot(ref_coldata %>%
         group_by(timepoint) %>%
         summarise(n = n()) %>%
         drop_na()) +
  geom_bar(aes(x = reorder(as.character(timepoint),-timepoint),
               y = n,
               fill = as.character(timepoint)),
           stat = "identity",
           color = "black",
           linewidth = 0.2) +
  theme_classic() +
  theme(legend.position = "none",
        #axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.line.y = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = full_spectrum_timepoint) +
  coord_flip()

ggsave("fig1a_reference_cells_barplot.png",
         dpi = 750,
         height = 3,
         width = 1,
         bg = "transparent")

# Figure 1b ---------------------------------------------------------------

ggplot(ref_coldata %>%
         group_by(timepoint,Oligo) %>%
         summarise(n = n()) %>%
         drop_na()) +
  geom_boxplot(aes(x = reorder(as.character(timepoint),-timepoint),
                   y = log10(n),
                   fill = as.character(timepoint)),
               color = "black",
               size = 0.1,
               outlier.size = 0.2,
               outlier.stroke = 0) +
  theme_classic() +
  scale_y_continuous(breaks = c(2,3,4), limits = c(1.5,4.5)) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = full_spectrum_timepoint) +
  coord_flip()

ggsave("fig1b_reference_cells_boxplot.png",
         dpi = 750,
         height = 3,
         width = 0.5,
         bg = "transparent")

# Figure 1c ---------------------------------------------------------------

# Individual UMAP 24hours 
ggplot(ref_coldata %>%
         sample_n(size = dim(ref_coldata)[1]) %>%
         filter(!is.na(timepoint))) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_2),
             color = "grey80",
             stroke = 0,
             size = 0.1) +
  geom_point(data = ref_coldata %>%
               filter(Oligo == "24h_ctrl-inj_P10_H7"),
             aes(x = umap3d_1,
                 y = umap3d_2,
                 color = as.character(timepoint)),
             stroke = 0,
             size = 0.125) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = full_spectrum_timepoint)

ggsave("fig1c_reference_umap_indiv_24.png",
         dpi = 600,
         height = 1,
         width = 1,
         bg = "transparent")


# Individual UMAP 48 hours 
ggplot(ref_coldata %>%
         sample_n(size = dim(ref_coldata)[1]) %>%
         filter(!is.na(timepoint))) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_2),
             color = "grey80",
             stroke = 0,
             size = 0.1) +
  geom_point(data = ref_coldata %>%
               filter(Oligo == "48h_ctrl-inj_P2_D11"),
             aes(x = umap3d_1,
                 y = umap3d_2,
                 color = as.character(timepoint)),
             stroke = 0,
             size = 0.125) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = full_spectrum_timepoint)

ggsave("fig1c_reference_umap_indiv_48.png",
         dpi = 600,
         height = 1,
         width = 1,
         bg = "transparent")


# Figure 1d ---------------------------------------------------------------

vibrant.colors <-  
  c('#EE7733', '#0077BB', '#228833', '#33BBEE', '#EE3377', '#CC3311',
    '#AA3377', '#009988', '#004488', '#DDAA33', '#99CC66','#D590DD')

bright.colors = <-  
  c('#4477AA', 
    '#EE6677', 
    '#228833', 
    '#CCBB44', 
    '#66CCEE',
    '#AA3377', 
    '#BBBBBB')


num.colors.tissue <-  
  ref_coldata %>%
  filter(!is.na(tissue)) %>%
  pull(tissue) %>%
  unique() %>%
  sort() %>%
  length()

tissue.colors <- colorRampPalette(c(vibrant.colors,bright.colors))(num.colors.tissue) 

# plot cells by time
ggplot(ref_coldata %>%
         sample_n(size = dim(ref_coldata)[1]) %>%
         filter(!is.na(timepoint))) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_2),
             color = "black",
             stroke = 0,
             size = 0.15) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_2,
                 color = as.character(timepoint)),
             stroke = 0,
             size = 0.11) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = full_spectrum_timepoint)

ggsave("fig1d_reference_umap_timepoint.png",
         dpi = 750,
         height = 1.5,
         width = 1.5,
         bg = "transparent")

# plot cells by tissue
ggplot(ref_coldata %>%
         sample_n(size = dim(ref_coldata)[1]) %>%
         filter(!is.na(timepoint)) %>%
         group_by(tissue) %>%
         add_tally() %>%
         arrange(-n)) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_2),
             color = "black",
             stroke = 0,
             size = 0.2) +
  geom_point(aes(x = umap3d_1,
                 y = umap3d_2,
                 color = as.character(tissue)),
             stroke = 0,
             size = 0.15) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = tissue.colors)

ggsave("fig1d_reference_umap_tissue.png",
         dpi = 1200,
         height = 3,
         width = 3,
         bg = "transparent")


# Figure 1e -------------------------------------------------------

# load broad cell type variance tests
celltype_disps <- fread("data/ref_broad_celltype_disps.csv", 
                        sep = ",", stringsAsFactors = F, 
                        data.table = F)

# set thresholds
sig_thresh = 0.05
times = c(20, 22, 24, 42)

plot_df = celltype_disps %>% 
  filter(cells_per_embryo > 3) %>% 
  filter(timepoint %in% times)

# plot variance faceted by time
ggplot(data=plot_df) + 
  geom_line(aes(cells_per_embryo, model_fit), color="black", data=plot_df) + 
  geom_ribbon(aes(cells_per_embryo, ymin=model_fit_lower, ymax=model_fit_upper), alpha=0.2) +
  geom_point(aes(cells_per_embryo, cells_per_embryo_cv, color=cv_z_stat_p_val < sig_thresh), size = 0.75) + 
  geom_linerange(aes(cells_per_embryo, 
                     ymin=cells_per_embryo_cv - cells_per_embryo_cv.stderr, 
                     ymax=cells_per_embryo_cv + cells_per_embryo_cv.stderr),
                 data=plot_df %>% filter(cv_z_stat_p_val < sig_thresh)) +
  facet_wrap(~timepoint, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  scale_y_log10() +
  scale_x_log10() +
  theme_classic() +
  theme(legend.position = "none")

ggsave("fig1e_reference_variance-sig_timepoint-select.pdf",
       width = 5,
       height = 2.2)

# Related extended data figures -------------------------------------------

# Extended data figure 6  ---------------------------------------------------------------

# load broad cell type variance tests
celltype_disps <- fread("code_to_submit/ref_broad_celltype_disps.csv", 
                        sep = ",", stringsAsFactors = F, 
                        data.table = F)

ggplot(celltype_disps %>% 
         filter(cells_per_embryo > 3), 
       aes(x = reorder(cell_group, cv_test_stat, FUN = median), y = cv_test_stat, fill = cell_group)) + 
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  theme_classic() +
  xlab("Cell type (broad)") +
  ylab("Cell-count variance (CV test statistic)") +
  theme(legend.position = "none", 
        axis.text = element_text(size = 6, color = "black"), 
        axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"))

ggsave(filename = "EDF6_ref_variance_ranked_boxplot.pdf", width = 8, height = 8)

