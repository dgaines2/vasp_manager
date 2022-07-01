#! /bin/bash

#SBATCH -N {n_nodes}
#SBATCH -p {queuetype}
#SBATCH -J b{jobname}
#SBATCH -A {allocation}
#SBATCH -t {walltime}
#SBATCH -L SCRATCH
#SBATCH -C knl,quad,cache

#OpenMP settings:
ulimit -s unlimited
export OMP_NUM_THREADS=1

#run the application:
module load {vasp_module}

starttime=$(date +%s)

mpitasks=$(echo "$SLURM_JOB_NUM_NODES * {ncore_per_node}" |bc)
for p in strain*; do
    cd $p
    srun -n $mpitasks -c 4 --cpu_bind=cores vasp_std > stdout.txt 2> stderr.txt
    cd ..
done

stoptime=$(date +%s)
tottime=$(echo "$stoptime - $starttime" | bc -l)
echo "total time (s): $tottime"
to_hours=$(echo $tottime/3600 | bc -l)
echo "total time (hr): $to_hours"
