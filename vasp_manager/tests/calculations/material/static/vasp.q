#! /bin/bash

#SBATCH -N 1
#SBATCH -q regular
#SBATCH -J sNaCl
#SBATCH -A m1673
#SBATCH -t 1:00:00
#SBATCH -C cpu
#SBATCH --mem=0

#OpenMP settings:
ulimit -s unlimited
export OMP_NUM_THREADS=1

#run the application:
module load vasp/6.3.2-cpu

starttime=$(date +%s)

mpitasks=$(echo "$SLURM_JOB_NUM_NODES * 128" |bc)
srun -t 0:59:00 -u -n $mpitasks --cpu_bind=cores vasp_std > stdout.txt 2> stderr.txt

stoptime=$(date +%s)
tottime=$(echo "$stoptime - $starttime" | bc -l)
echo "total time (s): $tottime"
to_hours=$(echo "scale=3; $tottime/3600" | bc -l)
echo "total time (hr): $to_hours"
