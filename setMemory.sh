ulimit -s unlimited
export OMP_STACKSIZE=100g
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-36}
