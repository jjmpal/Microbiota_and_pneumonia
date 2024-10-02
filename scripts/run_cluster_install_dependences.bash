#Installs dependences required by R packages to the conda environment.
#In reality these were run one by one. Running all at once is not tested.

cmd="source /homes/akauko/.bashrc
     conda activate finriskR_env
     conda install -c conda-forge r-rcpp
     conda install -c conda-forge lxml
     conda install -c conda-forge cmake
     conda install -c conda-forge gsl
     conda install -c conda-forge gmp
     conda install -c conda-forge r-rmpfr"

grun.py -n install_dependences -q highmem.q -c "$cmd"

