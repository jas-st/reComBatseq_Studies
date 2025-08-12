# reComBat-seq Studies
Summary of the pipelines and experiments associated with `reComBat-seq`

## File descriptions

### Simulation Studies

#### Hyperparameter Study

+ **sim_DEpipe.R, sim_DEpipe_helpers.R**
    + Pipeline and helper functions for simulations. Run sim_DEpipe.R to produce the simulation results. 
    + Usage: ``Rscript sim_DEpipe.R <mean batch effect> <dispersion batch effect> <total number of samples>``
    + True and false positive rates will be stored in CSV files. 
    + Modify the parameters to change the study design, level of biological signal, sequencing depth, etc.
+ **qsub_simDEpipe.py**
    + Script to run (qsub to cluster) multiple experiments
+ **visualize.R, visualize_helpers.R**
    + Script and helper functions to visualize the simulation results. Run visualize.R to generate the plot based on the CSV result files. 
    + Change the paths to files if necessary.

### Real data application

+ **gfrn_application.R, gfrn_DE.R, gfrn_helpers.R**
    + Script and helper functions for application example on the GFRN signature dataset. Run gfrn_application.R for the PCA analysis. Run gfrn_DE.R for differential expression analysis.
    + Change the paths to files at the top of the script, if necessary.
+ **signature_data.rds**
    + RDS object for the cleaned signature dataset, published[4] and used in our previous work[5,6].
+ **ras-pathway-gene-names.csv**
    + genes in RAS signaling pathway, obtained from [NCI website](https://www.cancer.gov/research/key-initiatives/ras/ras-central/blog/2015/ras-pathway-v2).

6. Zhang, Y., Jenkins, D. F., Manimaran, S., & Johnson, W. E. (2018). Alternative empirical Bayes models for adjusting for batch effects in genomic studies. *BMC bioinformatics*, 19(1), 262.
