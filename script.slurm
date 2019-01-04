#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
###SBATCH --nodes=2   # number of nodes
#SBATCH -J "sound"   # job name
###SBATCH --exclusive

srun hostname -s > hostfile

julia --machine-file hostfile script.jl
