## BF591-Final

The final project is an R Shiny application that implements multiple bioinformatics processes in R, including Sample Information Exploration, Counts Matrix Exploration, Differential Expression, and Gene Set Enrichment Analysis, each in different tabs.

### Shiny App Features

Each tab has the following features:

#### Sample Information Exploration

- Summary of the sample 
- Sample information 
- Histograms of continuous variables 

#### Counts Matrix Exploration (filtered)

- Table summarizing the effect of filtering 
- Diagnostic scatter plots 
- Clustered heatmap 
- Plot of PCA 

#### Differential Expression

- Table of differential expression genes
- Volcano plots 

#### Gene Set Enrichment Analysis

- Barplot of fgsea NES for top pathways 
- Downloadable data table of all, positive, or negative NES pathways 
- Scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis 

### Data

The data used in the app is from Labadorf et al. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4670106/), which utilized RNA Sequence Analysis of Human Huntington Disease Brain in the prefrontal cortex of 20 HD patients and 49 neuropathologically normal controls.

The count matrix and DEGs results were downloaded from here (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810). The summary data was downloaded using SRA run selector. fGSEA analysis was generated with DEGs results and analyzed with the fGSEA.R file, which filtered padj<0.05 and resulted in 5,480 genes. The gene set used is the Canonical pathways 2023 version (c2.cp.v2023.1.Hs.symbols.gmt), resulting in 688 pathways.