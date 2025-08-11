# Job Scripts

For each calculation, its job is handled by
[JobManager][vasp_manager.job_manager.JobManager]. This class can submit a job,
store its jobid, and monitor its progress.

The creation of the actual creation of the `vasp.q` SLURM submission script is
handled by [VaspInputCreator][vasp_manager.vasp_input_creator.VaspInputCreator]
using settings from `calc_config.json` and `computing_config.json`.

!!! note "Implementation Details"
    The [`vasp.q` file](https://github.com/dgaines2/vasp_manager/blob/main/vasp_manager/static_files/vasp.q)
    is just a template string that needs to be passed 3 separate sections:
    ```title="vasp.q"
    --8<-- "vasp_manager/static_files/vasp.q"
    ```

    1. `{sbatch_params}`: sbatch tags (e.g. `#SBATCH -N`)
    2. `{preamble}`: environment or module settings (e.g. `export OMP_NUM_THREADS=1`)
    3. `{command}`: the `mpi` or `srun` command to launch `VASP`

    The actual commands for each supercomputer and job type are determined by
    [`q_mapper.json`](https://github.com/dgaines2/vasp_manager/blob/main/vasp_manager/static_files/vasp.q)
    which lists the path of a yaml file that contains the `sbatch_params`,
    `preamble`, and `command` to use.

    You can modify the exact commands in the yaml file of your choice if needed.
    For reference, here's the yaml file for normal jobs on Perlmutter:
    ```title="vasp.yml"
    --8<-- "vasp_manager/static_files/vasp.yml"
    ```
