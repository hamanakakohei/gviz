#!/usr/bin/env bash
#
# 目的
#
# 入力
# 0.
# - gtf: 遺伝子名から遺伝子領域を取る（1.で撮影範囲とする）
# - 図示したいeqtl-egeneの組み合わせファイル
# 1.
# - 0で作ったファイル
# - vcf
# - sample - bam path
# - sample - normalization factor
# - gene model用tx db
# - 撮影時マージン
#
# 出力
# - gvizでcoverage, sashimi画像を撮る
set -euo pipefail


# 0
awk -F'\t' -v OFS='\t' '
  FNR==NR && $3=="gene" {
    if (match($9, /gene_name "([^"]+)"/, a))
      coord[a[1]] = $1 "\t" $4 "\t" $5 "\t" $7
    next
  }
  {
    if ($1 in coord)
      print $2, $1, coord[$1]
  }' \
  inputs/gm24_gencodev47.chr.scaffold.444141isoforms.geneid_corrected.geneid_refined.add_txType_geneType_geneRow.sort.gtf \
  inputs/trait-associated_new_gene_qtl.txt \
  > inputs/eqtl_egene_egeneRegion.txt


# 1
STA=1
END=$(wc -l < inputs/eqtl_egene_egeneRegion.txt)
#END=3

sbatch \
  --cpus-per-task=1 \
  --array=${STA}-${END}%28 \
  slurm/01_array.slurm \
  inputs/eqtl_egene_egeneRegion.txt \
  inputs/vcf/ALL.correctRefAlt.norm. \
  inputs/sample_bam.map \
  inputs/sample_normFactor.map \
  inputs/gm24_gencodev47.chr.scaffold.444141isoforms.geneid_corrected.geneid_refined.add_txType_geneType_geneRow.sort.rds \
  5000
  #--array=1,17,85,95,96,103,104,108,111,121,182,183,184,236,237,238,255,256,273,288%20 \
