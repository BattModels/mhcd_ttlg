#!/bin/bash
#SBATCH -J mhcd_ttlg_Eo_0.1
#SBATCH -n 128 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH -p RM
#SBATCH --mem=253000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o out.out # File to which STDOUT will be written
#SBATCH -e out.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mbabar@andrew.cmu.edu # Email to which notifications will be sent

echo "Job started on `hostname` at `date`"

julia eta_run_script.jl
#for i in $(seq -0.4 0.1 0.4)
#do
#    echo "Starting eta = $i \n"
#    julia script_cutout.jl i
#    echo "Finished eta = $i \n"
#done
    
echo " "
echo "Job Ended at `date`"
