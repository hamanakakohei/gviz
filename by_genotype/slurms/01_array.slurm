#!/bin/bash
#SBATCH --job-name=gviz
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate gviz_2025

INPUT_FILE=$1
VCF_PREFIX=$2
SAMPLE_BAM_MAP=$3
SAMPLE_NORM_MAP=$4
TXDB_RDS=$5
MARGIN=$6

read -r VARIANT GENE CHR STA END STRAND < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $INPUT_FILE)

echo $SLURM_ARRAY_TASK_ID
echo $INPUT_FILE
echo ${VCF_PREFIX}
echo ${CHR}
echo ${VARIANT}
echo ${GENE}
echo ${MARGIN}
echo $STA
echo $END
echo $STRAND

VCF=${VCF_PREFIX}${CHR}.vcf.gz
OUT=results/${VARIANT}.${GENE}.${MARGIN}.both.png
PLOT_STA=`expr $STA - $MARGIN`
PLOT_END=`expr $END + $MARGIN`

#Rscript scripts/01_gviz_rnaseq_profile_by_gt_group.R \
#       --vcf $VCF \
#       --variant $VARIANT \
#       --bam_map $SAMPLE_BAM_MAP \
#       --norm_map $SAMPLE_NORM_MAP \
#       --plt_chr $CHR \
#       --plt_sta $PLOT_STA \
#       --plt_end $PLOT_END \
#       --txdb_db $TXDB_RDS \
#       --ylim_range_strand both \
#       --width 1000 \
#       --height 150 \
#       --out $OUT
#       #--ylim_range_strand $STRAND \


Rscript scripts/01_gviz_rnaseq_profile_of_random_samples_of_each_gt.R \
        --vcf $VCF \
        --variant $VARIANT \
        --bam_map $SAMPLE_BAM_MAP \
        --sample_n 10 \
        --plt_chr $CHR \
        --plt_sta $PLOT_STA \
        --plt_end $PLOT_END \
        --txdb_db $TXDB_RDS \
        --width 1000 \
        --height 150 \
        --out $OUT
