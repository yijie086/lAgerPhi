#!/bin/bash
#SBATCH --account=eic
#SBATCH -N1
#SBATCH -n1
#SBATCH --time=12:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=2048       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=pcsim-1
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load parallel

mkdir -p /lcrc/project/eic/data/mc/pcsim-1
cd /lcrc/project/eic/data/mc/pcsim-1

parallel -j36 srun -l \
  bash /home/whit/projects/Pcsim/containers/slurm/run_singularity.sh \
  /home/whit/projects/Pcsim/examples/run_pcsim.sh -I {1} -E {2} -c -o {2}_on_{1} ::: $(seq 50 10 160) :::  $(seq 4 10)

wait
