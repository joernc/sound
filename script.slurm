#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=128   # number of processor cores (i.e. tasks)
###SBATCH --nodes=2   # number of nodes
#SBATCH -J "sound"   # job name
###SBATCH --exclusive

srun hostname -s > hostfile

julia --machine-file hostfile script.jl
