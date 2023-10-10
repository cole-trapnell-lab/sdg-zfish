# Overview 

## **Single cell RNA-seq of wild type and perturbed zebrafish development**

We explored the individual embryo and cell type-specific responses of zebrafish embryos to genetic and temperature perturbations using a scalable combinatorial-indexing scRNA-seq (sci-RNA-seq) approach. By barcoding hundreds-thousands of individual embryos, we statistically assess the changes in cell type abundances in response to each perturbation (3 temperatures, or 23 genetic loss of function experiments). Totalling ~3 million single cells across almost 2000 individual embryos, we put forward an analytical framework for high resolution phenotyping of diverse developmental traits in whole organisms. 

Our findings are published in [Saunders, Srivatsan, et al. _Nature, in press_ (2023)](https://doi.org/10.1101/2022.08.04.502764) and [Dorrity, et al. _Cell, in press_ (2023)](https://www.biorxiv.org/content/10.1101/2022.08.04.502669v1). 

For code and processed tables, please visit our [GitHub repository](https://github.com/cole-trapnell-lab/sdg-zfish).

## **Study Design**

Our studies produced three main datasets: 

**1. Reference atlas** - A hierarchically annotated, individually-resolved reference atlas, comprising 1,167 individuals and 1.2 million cells over 19 timepoints (from 18-96 hpf).

**2. Genetic perturbation atlas** - Two million cells from fish recieving 23 different genetic perturbations over multiple timepoints, totaling 645 individually barcoded animals (we collected 8 or more embryos per condition/timepoint).

Tissue | gene targets
----|----
Mesoderm | _cdx4, cdx4;cdx1a, tbxta, tbx16, tbx16;tbx16l, tbx16;msgn1, wnt3a;wnt8a, noto, smo, tbx1, hand2_
Central or peripheral nervous system | _egr2b, epha4a, hoxb1a, mafba, zc4h2, phox2a, foxi1, hgfa, met_ 
Neural crest lineages | _tfap2a, tfap2a;foxd3_  

**3. Temperature perturbation atlas** - 400,000 cells from fish raised at three different temperatures over three timepoints, totaling 288 individually barcoded animals.

## **Cell metadata breakdown**

Our cell metadata contains lots of information about time, perturbation, statistical metrics, and annotation. Here is a breakdown of those attributes according to our column names: 

* `timepoint`: The developmental stage in hours post fertilization (hpf) of embryos. Embryos were staged according to key landmarks according to [(Kimmel, et al (1995))](https://zfin.org/zf_info/zfbook/stages/index.html). Options are 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 72, and 96 hpf.
* `expt`: Denotes unique library preparation instances.
* `cell_type_sub`: The finest cell type annotation level.
* `cell_type_broad`: A broader level of cell type annotation than `cell_type_sub`, but still capturing all uniquely identified cell types.
* `tissue`: The tissue type annotation that contains the cell types.
* `germ_layer`: The germ layer of origin for each cell type, if known.
* `gene_target`: The genetic perturbation target. Controls include injected (scrambled) and unjected. Sibling controls are listed as `ctrl-<target>`, for null mutants included in the study.
* `mean_nn_time`: The mean time points of the nearest 15 neighbor cells
* `embryo`: The individual embryo barcode.
* `temp`: The growth temperature of the embryos (standard conditions are 28C).
* `pseudostage`: Embryo-level staging prediction by cell composition.

