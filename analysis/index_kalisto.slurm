#!/bin/bash
#SBATCH --job-name=kalisto_index
#SBATCH --partition=shared
#SBATCH --time=2-00:00:00

#SBATCH --cpus-per-task=4
#SBATCH --mem=128G

#SBATCH --error=hello-%A.err
#SBATCH --output=hello-%A.out

#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=knakats@emory.edu

module load lang/Python/3.9.6-GCCcore-11.2.0

pip install kb-python

cd /home/kenmn/nfs_scratch/nb12799275/analysis_exp/

kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno /home/kenmn/nfs_scratch/nb12799275/analysis_exp/Homo_sapiens.GRCh38.cdna.all.fa /home/kenmn/nfs_scratch/nb12799275/analysis_exp/homo_sapiens/Homo_sapiens.GRCh38.96.gtf
