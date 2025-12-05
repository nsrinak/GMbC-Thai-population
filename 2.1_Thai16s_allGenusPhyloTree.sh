#!/bin/bash

#SBATCH --job-name=phylo
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500000
#SBATCH --time=48:00:00
#SBATCH --output=phylo.out
#SBATCH --error=phylo.err
#SBATCH --partition=base
#SBATCH --qos=normal


muscle -align asv_genus_sequences_short.fasta -output asv_genus_sequences_short_alignment.fasta

Gblocks asv_genus_sequences_short_alignment.fasta -t=d

FastTree -gtr -nt asv_genus_sequences_short_alignment.fasta-gb > tree_asv