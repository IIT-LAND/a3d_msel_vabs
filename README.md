# A 3D approach to understanding heterogeneity in early developing autisms

This repository has all the code for the analyses in **Mandelli, V., Severino, I., Eyler, L., Pierce, K., Courchesne, E., & Lombardo, M. V. A 3D approach to understanding heterogeneity in early developing autisms.**

## Requirements

> Python = 3.8.1

Activate the environment in the folder and then Install all the libraries with:

> pip install -r requirements.txt

The stratification analysis also requires the `reval` Python library.

  + **reval** (https://github.com/IIT-LAND/reval_clustering) - Used for stability-based relative clustering validation analyses.


To run the R code, please use:

> R >=4.0


## Analysis pipeline

0) Download NDA data from NDA.

1) **_01_data_cleaning.Rmd**. Clean and merge the VABS and MSEL dataset from NDA. 

2) **_02_trts_split.ipynb**. Split the cleaned dataset into the train and test set before using reval.

3) **_03_batch_correction.Rmd**. Run batch correction for VABS.

4) **_04_run_reval.ipynb**. Run reval to identify subtypes.

5) **_05_plot_reval_output.Rmd**. Visualize reval results and graphically describe the subtypes.

6) **_06_reval_longitudinal_plots.Rmd**. Run longitidudinal analysis on NDA data.

7) **_07_clean_UCSD_data_classifier.Rmd**. Preprocess the UCSD datset before running the classification.

8) Use the application (**https://landiit.shinyapps.io/Autisms3D/**) to classify UCSD ACE subjects into subtypes.

9) **_08_ucsd_long.Rmd**. Run longitidudinal analysis on UCSD ACE data.

10) **_09_roi_analysis.Rmd**. Run ROI analysis on fMRI data.

11) **_10_pls.Rmd**. Run PLS analyses.

12) **_11_ari_enrichment_analysis.Rmd**. Run enrichment analysis with ARI genes.

