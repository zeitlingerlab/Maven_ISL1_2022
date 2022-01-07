# Introduction

This is the analysis code done by Melanie Weilert in the Zeitlinger Lab as part of the collaboration headed by Bonnie Maven in the Srivastava Lab for the submission "The multi-lineage transcription factor ISL1 controls cardiomyocyte cell fate through interaction with NKX2.5" to Cell Stem Cell. Below is a brief overview of the order of analysis, sources of the code, and links to the figures generated with posted code.

# Order of analysis

+ `Snakefile`: Processed the .BAM files in order to create merged Input and ChIP .bw files for BPNet training.
+ `analysis/1_model_training_and_grid_search.ipynb`: Test different model parameters to optimize training.
+ `analysis/2_map_motifs.ipynb`: After training a model, run TF-MoDISco on highly contribution regions, then map motifs to the genome
+ `analysis/3_motif_curation.Rmd`: Curate motifs from `2_*` in order to obtain a high confidence, denoised set for analysis.
+ `analysis/4_motif_summary.Rmd`: Summarize mapped motifs and actual signal/contribution across them.
+ `analysis/5_collect_genomic_perturbs.ipynb`: Generate predictions of ISL1 response upon motif mutation.
+ `analysis/6_pairwise_genomic_perturbs.Rmd`: Summarize mutation predictions in pairwise fashion.
+ `analysis/7_enhancer_predictions.Rmd`: Generate plots showing motif mappings across annotated enhancers.
+ `analysis/8_nko_analysis.Rmd`: Perform differential enrichment analysis across ISL1 WT peaks to see which peaks increase/decrease in signal upon NKX2.5 knockdown. After ISL1 WT peaks are classified with their mutation response, compare the enrichment of each motif across these different mutation responses.

# Figure links

+ figure2C-[MYL4 enhancer motifs](figures/7_enhancer_predictions/MYL4%20upstream.pdf)
+ figure2C-[PHOX2b enhancer motifs](figures/7_enhancer_predictions/PHOX2B%20downstream%202.pdf)
+ figure2C-[SETD5 enhancer motifs](figures/7_enhancer_predictions/SETD5%20gene%20body%201.pdf)
+ figure2D-[TF-MoDISco motif summary](figures/4_motif_summary/figure2D-curated_motif_metaplot.pdf)
+ figure2D-[reverse complement of NKX-alt](figures/2_map_motifs/figure2d-nkxalt_logos_revcomp.pdf)
+ figure3C-[In silico genomic perturbation summaries](figures/6_pairwise_genomic_perturbs/fc_med_cat_perturbs_500bp.pdf)
+ figure5A-[DESeq2 logFC of mutation response in ISL1 ChIP peaks](figures/8_nko_analysis/figure5A-peak_enrichments.pdf)
+ figure5B-[EXOSC2 enhancer motifs](figures/7_enhancer_predictions/EXOSC2%20downstream.pdf)
+ figure5C-[Motif enrichment across different ISL1 responses upon NKX mutation](figures/8_nko_analysis/figure5C-motif_enrichments_across_peak_states.pdf)
+ figureS3C-[Validation metrics for optimized BPNet model]() #TODO

# Intermediate data

Intermediate data files that are substituted as symlinks in the GitHub repo can be found in Zenodo at the following link (). #TODO upon publication

Below is a description of the intermediate files that can be found here.

+ data/bw: `.bw` coverage files generated via the `Snakefile`
+ bed: `.bed` files generated via analysis performed by Angelo
+ models: model `.h5` and associated files for the hyperoptimization analysis as well as the final model
+ preds: contribution (`.h5` and `.bw`) and prediction scores (`.bw`) generated across regions the model was trained on
+ modisco: TF-MoDISco results and mapped motifs from CWM-scanning generated on the contribution scores originating from the trained model
