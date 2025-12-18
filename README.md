# RNA-seq Gene Expression Analysis

This repository contains scripts, data, and models for the RNA-seq gene expression analysis presented in the paper titled "Immune signals dominate transcriptomic responses to developmental and life-history transitions in Antarctic fur seals", including filtering, differential expression, and WBC-adjusted analysis. The workflow is designed for reproducible analyses in R and Quarto. The raw data is available at *link*.

## Repository structure

### Data Files

| File | Description |
|-----------------------|------------------------------------------------|
| `rnaseq_filtered_Liv.Rdata` | Filtered RNA-seq dataset for all samples |
| `sample_sheet_rna_Liv.csv` | Sample metadata including IDs and library type. |
| `WBCcounts_GEP.xlsx` | White blood cell counts to control for leucocyte composition |

### Scripts

| File | Description |
|------------------------------------|------------------------------------|
| 0.RNAseq.filtering.R | Filters raw RNA-seq data |
| 1.Gene.expression.analysis.qmd | Main gene expression analysis |
| 2.Gene.expression.analysis.WBCadjusted.qmd | WBC-adjusted analysis (presented in publication) |

### Folders

| Folder | Description |
|----------------------------|--------------------------------------------|
| `GenesGO/` | Gene Ontology results and gene lists (all, upregulated, downregulated) |
| `Models/` | Saved model objects (`vfit_tp.rds`, etc.) for reproducibility |
| `Plots/` | Output plots from analyses (volcano plots, MA plots, etc.) |

#### Contact details

This script was run by Bernice Sepers (bioinformatics pipeline) and Ane Liv Berthelsen (analysis). If you have any questions, please contact: ane_liv.berthelsen[at]uni-bielefeld.de
