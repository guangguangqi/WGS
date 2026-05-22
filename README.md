Running Python Multi-Omics Bioinformatics Pipeline

This study establishes an end-to-end, high-performance Python workflow designed to process raw The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA: portal.gdc.cancer.gov/projects/TCGA-BRCA ) RNA-Seq datasets and isolate robust downstream biological pathways without the computational overhead of large-scale matrix operations.

Input:  The raw master matrix was parsed using TCGA barcodes to differentiate sample source types (tumor vs. healthy controls). To optimize local computational memory footprint and enhance script reproducibility for GitHub integration, a rigorous 6 vs. 6 sample subsetting strategy was executed (6 Tumor vs. 6 Healthy Normal uploaded in the input data) and   4 Tumor Vs 4 Healthy database were used for the analysis.
Raw RNA-Seq expression matrices downloaded directly from the TCGA-BRCA public repository. Individual sample files were programmatically merged, aligned, and converted locally into an integrated master count matrix containing all available clinical cases.  

Task: . Downstream data integration, formatting, and differential gene expression (DGE) statistical testing were handled natively in a Python environment using anndata(v0.9.2) and pydeseq2(v0.4.1). Functional enrichment and pathway analysis were subsequently evaluated using gseapy(v1.1.4).

Outputs:  A structured, lightweight subsetted count matrix matching the exact structure of the parent TCGA file format, along with a comprehensive table of statistically significant Differentially Expressed Genes (DEGs). Pathway analysis revealed a profound upregulation of core cell-cycle checkpoints and chromosome segregation machinery in tumor tissue. 
The top enriched biological processes include Spindle Assembly Checkpoint Signaling (GO:0071173), Mitotic Spindle Assembly Checkpoint Signaling (GO:0007094), and Negative Regulation of Mitotic Metaphase/Anaphase Transition (GO:0045841), all demonstrating high Combined Scores (>200) and strong statistical significance score -log_10 (FDR) >5.


Verification Metric: pytest_verify = 1.0. All custom matrix operations, local concatenation modules, and environment dependencies passed validation tests cleanly.
