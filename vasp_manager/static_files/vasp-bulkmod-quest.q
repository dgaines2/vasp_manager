#! /bin/bash

#SBATCH -N {n_nodes}
#SBATCH -n {n_procs}
#SBATCH -p {queuetype}
#SBATCH -J {jobname}
#SBATCH -A {allocation}
#SBATCH -t {walltime}
#SBATCH --mem=0

#OpenMP settings:
ulimit -s unlimited
export OMP_NUM_THREADS=1

#run the application:
module load {vasp_module}

starttime=$(date +%s)

mpi_tasks=$(echo "$SLURM_NTASKS * 24/28" | bc)
for p in strain*; do
    cd $p
    mpirun -np $mpi_tasks vasp_std > stdout.txt 2> stderr.txt
    cd ..
done

stoptime=$(date +%s)
tottime=$(echo "$stoptime - $starttime" | bc -l)
echo "total time (s): $tottime"
to_hours=$(echo "scale=3; $tottime/3600" | bc -l)
echo "total time (hr): $to_hours"
