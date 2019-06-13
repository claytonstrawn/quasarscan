#PBS -S /bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:model=hsw
#PBS -l walltime=00:10:00
#PBS -j oe
#PBS -o my_job_output_file.txt

