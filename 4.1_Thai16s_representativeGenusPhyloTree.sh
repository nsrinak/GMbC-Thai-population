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


muscle -align asv_genus_sequences_short_nogap.fasta -output asv_genus_nogap_alignment.fasta

Gblocks asv_genus_nogap_alignment.fasta -t=d

FastTreeMP -nt -gtr -gamma -spr 4 -mlacc 2 -slownni -boot 500 asv_genus_nogap_alignment.fasta-gb > tree_asv_genus_nogap-gb_boot
