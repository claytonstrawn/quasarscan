#PBS -S /bin/bash
#PBS -l select=4:ncpus=24:model=has
#PBS -l walltime=00:15:00
#PBS -o my_job_output_file.txt
#PBS -q debug

bash quasarscan/batch_scripts/run_continuously_nasa.sh 96