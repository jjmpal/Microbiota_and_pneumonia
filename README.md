# Microbiota and pneumonia
Microbiota and pneumonia using FINRISK 2002 data. The code is based on Helius replication code by Bob Kullberg. 

**Running the the analysis**

1. Prepare the environment, if not prepared before. 
  * Install finriskR_env conda environment with R4.3.1: https://github.com/finrisk2002/finriskR/tree/main/installation
  * Install additional dependences at the atlas node: `bash scripts\run_cluster_install_dependences.bash`  
  * Install R libraries at the atlas node `bash scripts\run_cluster_install_packages.bash`
  
2. Run the analysis at a HPC
  * Creates automatically output as *md file*.
  * Requires also manually created figure directories figs_mg and data_mg.
  * Requires functions_yingtools2.R (https://github.com/ying14/yingtools2), which includes some functions from the yingtools2 package. We had earlier problems with package dependences so we sourced individual functions instead of whole package. 

 Generate data: Run command `Rscript scripts/gen_data_mg.R > data_mg/gen_data.log` at the cluster in Atlas.
```
mkdir data_mg
bash scripts\run_cluster_gen_data.bash
```

Main analysis: Run the command  `R -e "knitr::knit('Microbiota_and_pneumonia_mg2.Rmd')` at the node in the atlas
```
mkdir figs_mg
bash script\run_cluster_create_md.bash
```
3. Pull file *'microbiota_pneumonia_mg2.md* and *figs_mg directory* to local a computer with the pandoc and run:

```
R -e "rmarkdown::render('Microbiota_and_pneumonia_mg2.md')"
```

**Files**

```

├── README.md                                                # Overview
├── Microbiota_pneumonia_mg2.Rmd                             # Script to run the analysis (metagenome data, R4.3.1, finriskr_env environment)
├── scripts
  ├── functions.R                                               # Minor R functions for the analysis script
  ├── functions_yingtools2.R                                    # Functions from package Yingtools2
  ├── gen_data_mg.R                                             # Generates the data, metagenome
  ├── gen_data.R                                                # Generates the data, 16S
  ├── run_cluster_gen_data.bash                                 # Runs data generation at the cluster
  ├── run_cluster_create_md.bash                                # Runs analysis script at the cluster
  ├── run_cluster_install_dependences.bash                      # Installs required dependences to conda environment
  ├── run_cluster_install_packages.bash                         # Runs package installation at the cluster
  ├── install_packages.R                                        # Package installation script
  ├── convert_mia2phy.R                                         # Converts mia micriobiome data to phyloseq format; for Oleg, not needed

```
