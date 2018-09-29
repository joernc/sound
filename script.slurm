#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH -J "sound_1e-3"   # job name

julia -p 32 script.jl
