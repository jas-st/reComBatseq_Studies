#!/bin/bash
#
#SBATCH --job-name=recombatseq
#SBATCH --cpus-per-task=4
#SBATCH --time=01-00:00:00
#SBATCH --mem=1000GB
#SBATCH --output=correctfull.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel@lbi-netmed.com
#SBATCH --no-requeue

export OPENBLAS_NUM_THREADS=1
#conda activate recombatseq
#echo 'correcting tregs'
#Rscript sccorrect.Rscript -i skin.tregs.exp.tsv -m skin.tregs.meta.tsv -b patient_id --nproc 4
#echo 'correcting filtered tregs'
#Rscript sccorrect.Rscript -i skin.tregs.exp.filtered.tsv -m skin.tregs.meta.tsv -b patient_id --nproc 4
#echo 'correcting all sample'
#Rscript sccorrect.Rscript -i skin.inflammatory.exp.sample.tsv -m skin.inflammatory.meta.sample.tsv -b patient_id --nproc 4
echo 'correcting all'
Rscript sccorrect.Rscript -i skin.inflammatory.exp.tsv -m skin.inflammatory.meta.tsv -b patient_id --nproc 4
