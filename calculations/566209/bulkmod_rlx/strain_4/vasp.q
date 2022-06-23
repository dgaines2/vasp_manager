#! /bin/bash

#SBATCH -N 1
#SBATCH -p regular
#SBATCH -J 566209_4
#SBATCH -t 00:10:00
#SBATCH -A m1673
#SBATCH -L SCRATCH
#SBATCH -C knl,quad,cache

#OpenMP settings:
ulimit -s unlimited
export OMP_NUM_THREADS=1

#run the application:
module load vasp/20181030-knl

starttime=$(date +%s)
ulimit -s unlimited

mpitasks=$(echo $SLURM_JOB_NUM_NODES*64|bc)
srun -n $mpitasks -c 4 --cpu_bind=cores vasp_std > stdout.txt 2> stderr.txt

stoptime=$(date +%s)
tottime=$(echo "${stoptime} - ${starttime}" | bc -l)
echo "total time (s): $tottime"
to_hours=$(echo $tottime/3600 | bc -l)
to_charge_factor=$(echo "scale=3; $to_hours*$SLURM_JOB_NUM_NODES*80/1" | bc)
echo "SU: $to_charge_factor"
