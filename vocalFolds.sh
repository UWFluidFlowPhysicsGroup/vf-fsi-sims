#!/bin/bash
#SBATCH --account=def-peterson
#SBATCH --ntasks-per-node=96      # number of MPI processes
#SBATCH --nodes=1
#SBATCH --mem=576G      # memory; default unit is megabytes
#SBATCH --time=2-00:00           # time (DD-HH:MM)
#SBATCH --output=./slurm-sim-%j.out  # create log file for outputs, outputs data from console

# displays outputs to terminal use for debugging, but dont need for Nibi cluster simulation
# set -x

# Script used to run simulation on Nibi cluster

# Everything present in projects/def-peterson/nhewko is accessible on computational cluster,
# compile everything beforehand through ssh terminal to Nibi, then can run without having to
# recompile.

# if short run simulations work through ssh access, then it should work with the compute cluster.
# only need to run 'mpiexec -n 10 main' line or similar for compute cluster batch job.

module load dealii/9.4.1
export LD_PRELOAD=$EBROOTMETIS/lib/libmetis.so
srun main


