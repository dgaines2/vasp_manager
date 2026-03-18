#!/bin/bash

#SBATCH -N 1
#SBATCH -q regular
#SBATCH -J bBCC_Ti_s-4
#SBATCH -A m4545
#SBATCH -t 01:00:00
#SBATCH -C cpu
#SBATCH --mem=0

ulimit -s unlimited
export OMP_NUM_THREADS=1
module load vasp/6.4.3-cpu
mpitasks=$(echo "$SLURM_JOB_NUM_NODES * 128" | bc)

starttime=$(date +%s)

srun -t 00:59:00 -u -n $mpitasks --cpu_bind=cores vasp_std > stdout.txt 2> stderr.txt

stoptime=$(date +%s)
tottime=$(echo "$stoptime - $starttime" | bc -l)
echo "total time (s): $tottime"
to_hours=$(echo "scale=3; $tottime/3600" | bc -l)
echo "total time (hr): $to_hours"
