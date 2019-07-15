#!/bin/bash

#SBATCH -t 5-0
#SBATCH --mem 2G
#SBATCH -p medium

dbp_relax_b 1aay.pdb
