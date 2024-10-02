# Run cluster
cmd="source /homes/akauko/.bashrc; conda activate finriskR_env; Rscript scripts/gen_data_mg.R > data_mg/gen_data.log"
grun.py -n run_gen_data -q highmem.q -c "$cmd"

