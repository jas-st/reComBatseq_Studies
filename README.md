# reComBat-seq Studies
Summary of the pipelines and experiments associated with `reComBat-seq`

## File descriptions

### Simulation Studies

#### Hyperparameter Study

+ **results_pipeline.R**
    + Runs a hyperparameter grid search varying `alpha` and `lambda` for different mean (`mbfoldX`) and dispersion (`dbfold`) batch effects.
    + Outputs CSV files (`mbfoldX_dbfoldY.csv`) storing precision values for each parameter combination (results stored in `result_csvs.zip`).
+ **visualisations.R**
    + Creates the heatmap based on the result dataframes.
    + Results visualised in `hyperparams_heatmap2.png`.

#### Gene Count Study

+ **results_pipeline_step1.R**
    + Simulates RNA-seq datasets with and without batch effects - 1024 samples and genes varying between 128 and 20 171.
    + Applies `ComBat-seq` and `reComBat-seq` for batch correction, runs DE analysis (edgeR) on all data, computes precision/TPR/FPR.
    + Saves raw counts, normalized counts (voom/TMM), confounded covariate matrices and corrected count matrices as CSV files.
+ **results_pipeline_step2_pycombat.ipynb, results_pipeline_step2_recombat.ipynb**
    + Applies `reComBat` and `pyComBat-seq` for batch correction on the simulated data from the previous step.
    + For `reComBat` it uses the transformed normalized counts instead of the raw matrix.
    + A modified version of `pyComBat-seq` is used, that allows confounded matrices.
    + Runs LDA analysis on all data (raw, `ComBat-seq`, `reComBat-seq`, `reComBat`, `pyComBat-seq`).
+ **results_pipeline_step3**
    + Runs DE analysis on the `reComBat` (limma) and `pyComBat-seq` (edgeR) corrected data.
+ **visualisations_pca.R**
    + Creates PCA plots for the simulated and corrected data.
    + All plots are stored in the **PCA Plots** directory.
+ **visualisations_summary.R**
    + Creates the summary plot for all the methods.
    + Precision/LDA/Time are illustrated in **gene_plot.png**, TPR/FPR in **genes_plot_2.png**.


 
#### Sample Count Study

+ **results_pipeline_step1.R**
    + Simulates RNA-seq datasets with and without batch effects - 1024 genes and samples varying between 16 and 16 384.
    + Applies `ComBat-seq` and `reComBat-seq` for batch correction, runs DE analysis (edgeR) on all data, computes precision/TPR/FPR.
    + Saves raw counts, normalized counts (voom/TMM), confounded covariate matrices and corrected count matrices as CSV files.
+ **results_pipeline_step2_pycombat.ipynb, results_pipeline_step2_recombat.ipynb**
    + Applies `reComBat` and `pyComBat-seq` for batch correction on the simulated data from the previous step.
    + For `reComBat` it uses the transformed normalized counts instead of the raw matrix.
    + A modified version of `pyComBat-seq` is used, that allows confounded matrices.
    + Runs LDA analysis on all data (raw, `ComBat-seq`, `reComBat-seq`, `reComBat`, `pyComBat-seq`).
+ **results_pipeline_step3**
    + Runs DE analysis on the `reComBat` (limma) and `pyComBat-seq` (edgeR) corrected data.
+ **visualisations_pca.R**
    + Creates PCA plots for the simulated and corrected data.
    + All plots are stored in the **PCA Plots** directory.
+ **visualisations_summary.R**
    + Creates the summary plot for all the methods.
    + Precision/LDA/Time are illustrated in **sample_plot.png** and **sample_plot_3.png**, TPR/FPR in **genes_plot_2.png**.
 
#### Batch Count Study

+ **results_pipeline_step1.R**
    + Simulates RNA-seq datasets with and without batch effects - 6000 samples, 1024 genes and batches varying between 2 and 40.
    + Applies `reComBat-seq` for batch correction, runs DE analysis (edgeR) on all data, computes precision/TPR/FPR.
    + Saves raw counts, normalized counts (voom/TMM), confounded covariate matrices and corrected count matrices as CSV files.
+ **results_pipeline_step2.ipynb**
    + Applies `reComBat` for batch correction on the simulated data from the previous step.
    + Uses the transformed normalized counts instead of the raw matrix.
    + Runs LDA analysis on all data (raw, `reComBat-seq`, `reComBat`).
+ **results_pipeline_step3**
    + Runs DE analysis on the `reComBat` (limma) corrected data.
+ **visualisations_umap.ipynb**
    + Creates UMAP plots for the simulated and corrected data.
    + Example plots are **batch_comparison.png** and **batch_iterations_plot.png**
+ **visualisations_summary.R**
    + Creates the summary plot for all the methods.
    + Precision/LDA/Time are illustrated in **batch_plot.png**, TPR/FPR in **batch_plot_2.png**.
 
### Real Data

+ **real_data_pipeline_recombatseq.R**
    + Applies `reComBat-seq` for batch correction on the dataset (tissue-wise).
    + Saves the corrected count matrices as CSV files.
+ **real_data_pipeline_recombatseq.R**
    + Creates the UMAP plots for each tissue - raw and `reCombat`/`reComBat-seq` corrected.
    + The plots can be found in the `UMAP Plots` directory
