#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

#load conda environment
module load conda

source activate /mainfs/home/lag1e24/.conda/envs/rmats_4.1.2/

#rMATS
conda run rmats.py \
	--gtf /mainfs/BaralleLab/for_Luke/genome_data/gtf/gencode.v47.primary_assembly.annotation.gtf \
	--b1 "/mainfs/BaralleLab/for_Luke/puf60rmats/b1.txt"\
	--b2 "/mainfs/BaralleLab/for_Luke/puf60rmats/b2.txt" \
	--od "/mainfs/BaralleLab/for_Luke/puf60rmats/" \
	--tmp "/mainfs/BaralleLab/for_Luke/puf60rmats/tmp" \
	--libType fr-firststrand \
	--allow-clipping \
	--readLength 101 \
	--nthread 8


conda deactivate


