#! /bin/bash

#SBATCH -N 2
#SBATCH -n 56
#SBATCH -p short
#SBATCH -J rcAlSb
#SBATCH -A p31151
#SBATCH -t 02:00:00
#SBATCH --mem=0

#OpenMP settings:
ulimit -s unlimited
export OMP_NUM_THREADS=1

#run the application:
module load vasp/5.4.4-vtst-openmpi-4.0.5-intel-19.0.5.281

starttime=$(date +%s)

mpi_tasks=$(echo "$SLURM_NTASKS * 24/28" | bc)
mpirun -np $mpi_tasks vasp_std > stdout.txt 2> stderr.txt

stoptime=$(date +%s)
tottime=$(echo "$stoptime - $starttime" | bc -l)
echo "total time (s): $tottime"
to_hours=$(echo $tottime/3600 | bc -l)
echo "total time (hr): $to_hours"
