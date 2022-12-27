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

kb count -i index.idx -g t2g.txt -x 10xv2 --workflow lamanno --loom -c1 cdna_t2c.txt -c2 intron_t2c.txt /home/kenmn/nfs_scratch/nb12799275/fastq/SRR10156295_S1_L001_R1_001.fastq /home/kenmn/nfs_scratch/nb12799275/fastq/SRR10156295_S1_L001_R2_001.fastq