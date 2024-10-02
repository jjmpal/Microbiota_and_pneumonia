# Run cluster
cmd="source $HOME/.bashrc; conda activate finriskR_env; Rscript scripts/install_packages.R"
grun.py -n install_packages -q highmem.q -c "$cmd"

