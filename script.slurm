#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=64   # number of processor cores (i.e. tasks)
###SBATCH --nodes=2   # number of nodes
#SBATCH -J "sound"   # job name

srun hostname -s > hostfile

julia --machine-file hostfile script.jl
