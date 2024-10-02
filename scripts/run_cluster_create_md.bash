# Run cluster
cmd="knitr::knit('Microbiota_and_pneumonia_mg2.Rmd')"
grun.py -n create_md -q highmem.q -c "source $HOME/.bashrc; conda activate finriskR_env; R -e \"$cmd\""
