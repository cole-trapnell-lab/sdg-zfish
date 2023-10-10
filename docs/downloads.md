# Downloads

Here, you'll find data and code to reproduce analyses in the papers or explore the data for yourself.

For code and processed tables, please visit our [GitHub repository](https://github.com/cole-trapnell-lab/sdg-zfish).

### **Reference Atlas**

A processed and annotated monocle3 `cell_data_set` object (.RDS file) containing ~1.25M cells of wild type zebrafish development. This reference includes both the uninjected controls (this dataset only) and the injected-control embryos (also in the perturbation atlas).

[Reference atlas](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_reference_cds.RDS.gz){ .md-button .md-button--primary }

##### For individual parts of the dataset:
[Raw count matrix (.mtx)](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_reference_raw_counts.mtx){ .md-button }
[Cell metadata (.csv)](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_reference_cell_metadata.csv.gz){ .md-button }
[Gene metadata (.csv)](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_reference_gene_metadata.csv.gz){ .md-button }

### **Genetic Perturbation Atlas**

A processed and annotated monocle3 `cell_data_set` object (.RDS file) containing 2M cells from control-injected and perturbed zebrafish embryos. Timepoints: 18, 24, 36, 48, 72 hpf.

[Perturbation atlas](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_perturb_cds.RDS){ .md-button .md-button--primary }

##### For individual parts of the dataset:
[Raw count matrix (.mtx)](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_perturb_full_raw_counts.mtx){ .md-button }
[Cell metadata (.csv)](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_perturb_full_cell_metadata.csv.gz){ .md-button }
[Gene metadata (.csv)](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_perturb_full_gene_metadata.csv.gz){ .md-button }

##### Processed data subsets:

* Subclustered and finely annotated sensory cranial ganglia cells (plus Rohon-Beard neurons). This dataset contains just wild type and control-injected cells.  

[Cranial ganglia subset (.RDS)](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_reference_sensory_neuron_cds.RDS.gz){ .md-button }

* Subclustered and finely annotated notochord and parachordal cartilate cells (control-injected and wild type only).  

[Notochord subset (.RDS)](http://trapnell-lab-s3-zscape.s3-website-us-west-2.amazonaws.com/zscape_reference_notochord_cds.RDS.gz){ .md-button }

### **Temperature Atlas**

A processed and annotated monocle3 `cell_data_set` object (.RDS file) containing 400,000 cells from embryos raised at increased temperature (28&deg;C, 32&deg;C, and 34&deg;C).

[Temperature atlas](Shared drives/Trapnell Lab/ZSCAPE_site/data/zscape_hotfish_cds.rds){ .md-button .md-button--primary }

##### For individual parts of the dataset:
[Count matrix (.mtx)](Shared drives/Trapnell Lab/ZSCAPE_site/data/zscape_temp_full_raw_counts.mtx){ .md-button }
[Cell metadata (.csv)](Shared drives/Trapnell Lab/ZSCAPE_site/data/zscape_temp_full_cell_metadata.csv.gz){ .md-button }
[Gene metadata (.csv)](Shared drives/Trapnell Lab/ZSCAPE_site/data/zscape_temp_full_gene_metadata.csv.gz){ .md-button }