#!/bin/bash
#SBATCH --job-name=cellranger_nb
#SBATCH --partition=shared
#SBATCH --time=2-00:00:00

#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

#SBATCH --error=hello-%A.err
#SBATCH --output=hello-%A.out

#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=kenmn@emory.edu

export PATH=/home/kenmn/bioinfo_tools/cellranger-7.0.1:$PATH

bash crd_raw_config.sh

##run single sample through cellranger pipeline
if (($crd_loop==0))
then
cellranger count --id=$crd_jobid --fastqs=$crd_fastqs --sample=$crd_sample --transcriptome=$crd_transcriptome
fi


##run multiple samples through cellranger pipeline (for 2-3 samples, for bulk samples run multithread)
if (($crd_loop==1))
then
for crd_sample in $(cat crd_samples.txt)
do
cellranger count --id=$crd_jobid --fastqs=$crd_fastqs --sample=$crd_sample --transcriptome=$crd_transcriptome
done
fi

