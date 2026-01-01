#!/usr/bin/env Rscript
#
# 1. VCFから座標で抜き出す-->idで抜き出す-->ジェノタイプをdfにする
#
# 2. 各ジェノタイプからランダムに最大で10サンプル選んでBAMパスに変える
#
# 3. gviz絵を描く
library(optparse)
library(VariantAnnotation)
library(Gviz)
library(GenomicFeatures)
library(Rsamtools)


option_list <- list(
  make_option("--vcf", type="character", help="マルチサンプルVCF、指定したvariant_idの各サンプルのジェノタイプをここから抜き出す"),
  make_option("--variant", type="character", help="バリアントIDで「chr:100:A:G」みたいなフォーマットを想定して、座標とIDの両方でバリアント行をVCFから取り出す"),
  make_option("--bam_map", type="character", help="描くサンプルリストをBAMパスリストに変えるmapファイル"),
  make_option("--sample_n", type="integer", help="ジェノタイプごとに描くサンプルの数、ランダウに選ぶ"),
  make_option("--plt_chr", type="character"),
  make_option("--plt_sta", type="integer"),
  make_option("--plt_end", type="integer"),
  make_option(c("--txdb_db"), type="character", default=NULL, help="--gtfか--txdb_dbを指定する、描く遺伝子モデル"),
  make_option(c("--gtf"), type="character", help="--gtfか--txdb_dbを指定する"),
  make_option(c("--width"), type="integer", default=3000),
  make_option(c("--height"), type="integer", default=2000),
  make_option(c("--out"), type="character", default="out.gviz.png")
)
opt <- parse_args(OptionParser(option_list=option_list))

set.seed(1)

VCF = opt$vcf
VARIANT = opt$variant
BAM_MAP = opt$bam_map
SAMPLE_N = opt$sample_n
CHROM = opt$plt_chr
START = opt$plt_sta
END = opt$plt_end
TXDB_DB = opt$txdb_db
GTF = opt$gtf
WIDTH = opt$width
HEIGHT = opt$height
OUT = opt$out


# 1
parts <- strsplit(VARIANT, ":")[[1]]
chr <- parts[1]
pos <- as.integer(parts[2])

tbx <- TabixFile(VCF)
open(tbx)
param <- GRanges(seqnames = chr, ranges = IRanges(pos, pos))
vcf_sub <- readVcf(tbx, genome = "hg38", param = param)
close(tbx)

vcf_var <- vcf_sub[names(rowRanges(vcf_sub)) == VARIANT]

gt <- geno(vcf_var)$GT
gt_df <- data.frame(
  sample = colnames(gt),
  GT     = as.character(gt[1, ]),
  stringsAsFactors = FALSE
)

table(gt_df$GT, useNA = "ifany")


# 2. サンプルをジェノタイプごとにランダムに選んでジェノタイプ順にする
pick_samples <- function(df, gt, n=3) {
  sub <- df[df$GT == gt, ]
  if (nrow(sub) <= n) {
    sub
  } else {
    sub[sample(nrow(sub), n), ]
  }
}

gt_00 <- pick_samples(gt_df, "0/0", SAMPLE_N)
gt_01 <- pick_samples(gt_df, "0/1", SAMPLE_N)
gt_11 <- pick_samples(gt_df, "1/1", SAMPLE_N)
selected_samples <- rbind(gt_00, gt_01, gt_11)

bam_map <- read.table(BAM_MAP, header=TRUE, sep="\t")
selected_samples <- merge(selected_samples, bam_map, by="sample")

selected_samples$GT <- factor(
  selected_samples$GT,
  levels=c("0/0", "0/1", "1/1")
)

selected_samples <- selected_samples[
  order(selected_samples$GT), ]


# geneなど
source("scripts/igv_patch.orig.R")
libType = 'fr-firststrand'

if (!is.null(TXDB_DB)) {
  txdb <- loadDb(TXDB_DB)
} else {
  txdb <- makeTxDbFromGFF(GTF, format="gtf")
}

gene_track <- GeneRegionTrack(
  txdb,
  chromosome=CHROM,
  start=START,
  end=END,
  transcriptAnnotation="symbol",
  name="Gene Model")

axis_track <- GenomeAxisTrack()


# y軸の最大値を決める
get_max_cov <- function(bam, region) {

  param <- ScanBamParam(
    which = region,
    what  = c("pos", "qwidth"),
    flag  = scanBamFlag(isUnmappedQuery = FALSE)
  )

  aln <- scanBam(bam, param = param)[[1]]

  if (length(aln$pos) == 0) return(0L)

  # ★ region 内相対座標に変換
  rel_start <- aln$pos - start(region) + 1

  # region 外に出るものを除外
  keep <- rel_start > 0 & rel_start <= width(region)

  cov <- coverage(
    IRanges(start = rel_start[keep],
            width = aln$qwidth[keep]),
    width = width(region)
  )

  max(as.numeric(cov), na.rm = TRUE)
}

region <- GRanges(CHROM, IRanges(START, END))
cov_max <- max(
  sapply(selected_samples$bam, get_max_cov, region = region)
)

cov_max <- ceiling(cov_max * 1.2)


# トラックを作る
tracks <- list(axis_track, gene_track)

for (i in seq_len(nrow(selected_samples))) {

  bam <- selected_samples$bam[i]
  sample <- selected_samples$sample[i]
  GT <- selected_samples$GT[i]

  tracks <- c(
    tracks,
    DataTrack(
      bam,
      genome = 'hg38',
      chromosome = CHROM,
      importFunction=strandedBamImport,
      stream=TRUE, legend=TRUE,
      col=c("cornflowerblue","purple"),
      groups=c("Forward","Reverse"),
      name = paste(sample, GT),
      #ylim = c(0, cov_max),
      type = "hist"
    ),
    AlignmentsTrack(
      bam,
      isPaired = TRUE,
      name = sample,
      type = "sashimi",
      #sashimiHeight = 5,
      sashimiNumbers = TRUE,
      sashimiScore = 20
    )
  )
}

png(OUT, width = WIDTH, height = HEIGHT * (nrow(selected_samples)*2+5), res = 150)
plotTracks(
  tracks,
  chromosome = CHROM,
  from = START,
  to = END,
  ylim = c(0, cov_max),
  sizes = c(0.5, 2, rep(0.5, nrow(selected_samples)*2))
)
dev.off()
