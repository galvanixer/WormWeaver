#!/bin/bash
#SBATCH -N 1 --exclusive
#SBATCH --mem-per-cpu=4000
#SBATCH -p grant -A g2026a136c
#SBATCH -t 4-00:00:00
commandes=jobfile
hpc_multilauncher $commandes
