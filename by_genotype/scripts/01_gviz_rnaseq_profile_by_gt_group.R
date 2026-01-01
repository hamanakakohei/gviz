#!/usr/bin/env Rscript
#
# 1. VCFから座標で抜き出す-->idで抜き出す-->ジェノタイプをdfにする（ついでにジェノタイプ毎のサンプル数を出す）
#
# 2. サンプルの情報（bamパス、normalization factorなど）を付ける
#
# 3. 全サンプルの平均covを計算して、gviz 用 DataTrack を作る（FとRをもつGRanges 化）
#
# 4. gvizの用意を色々して描く
library(optparse)
library(VariantAnnotation)
library(Gviz)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)


option_list <- list(
  make_option("--vcf", type="character", help="マルチサンプルVCF、指定したvariant_idの各サンプルのジェノタイプをここから抜き出す"),
  make_option("--variant", type="character", help="バリアントIDで「chr:100:A:G」みたいなフォーマットを想定して、座標とIDの両方でバリアント行をVCFから取り出す"),
  make_option("--bam_map", type="character", help="描くサンプルリストをBAMパスリストに変えるmapファイル"),
  make_option("--norm_map", type="character", help="描くサンプルリストをnormalization factorリストに変えるmapファイル"),
  make_option("--sample_n", type="integer", help="ジェノタイプごとに描くサンプルの数、ランダウに選ぶ"),
  make_option("--plt_chr", type="character"),
  make_option("--plt_sta", type="integer"),
  make_option("--plt_end", type="integer"),
  make_option(c("--txdb_db"), type="character", default=NULL, help="--gtfか--txdb_dbを指定する、描く遺伝子モデル"),
  make_option(c("--gtf"), type="character", help="--gtfか--txdb_dbを指定する"),
  make_option(c("--ylim_range_strand"), type="character", default="both", help="both, +, or -"),
  make_option(c("--width"), type="integer", default=3000),
  make_option(c("--height"), type="integer", default=2000),
  make_option(c("--out"), type="character", default="out.gviz.png")
)
opt <- parse_args(OptionParser(option_list=option_list))

VCF = opt$vcf
VARIANT = opt$variant
BAM_MAP = opt$bam_map
NORM_FACTOR_MAP = opt$norm_map
SAMPLE_N = opt$sample_n
CHROM = opt$plt_chr
START = opt$plt_sta
END = opt$plt_end
TXDB_DB = opt$txdb_db
GTF = opt$gtf
YLIM_RANGE_STRAND = opt$ylim_range_strand
WIDTH = opt$width
HEIGHT = opt$height
OUT = opt$out


# 関数
get_stranded_cov <- function(bam, region) {

  param <- ScanBamParam(
    which = region,
    what = c("flag", "pos", "cigar", "strand", "qwidth")
  )

  sb <- scanBam(bam, param = param)[[1]]
  #  リードが1本もない場合、0カバレッジを返す
  if (length(sb$pos) == 0) {
    zero <- Rle(0L, width(region))
    return(list(
      plus  = zero,
      minus = zero
    ))
  }

  # FLAG 判定
  is_r1 <- bitwAnd(sb$flag, 0x40) != 0
  is_r2 <- bitwAnd(sb$flag, 0x80) != 0

  # strand
  st <- sb$strand

  A <- (is_r1 & st == "+") | (is_r2 & st == "-")
  B <- (is_r1 & st == "-") | (is_r2 & st == "+")

  # relative coordinate
  rel_start <- sb$pos - start(region) + 1
  keep <- rel_start > 0 & rel_start <= width(region)

  covA <- coverage(
    IRanges(rel_start[A & keep], width = sb$qwidth[A & keep]),
    width = width(region)
  )

  covB <- coverage(
    IRanges(rel_start[B & keep], width = sb$qwidth[B & keep]),
    width = width(region)
  )

  list(
    plus  = as.numeric(covA),
    minus = as.numeric(covB)
  )
}

# GT ごとに「正規化平均カバレッジ」を計算
get_mean_cov_by_gt <- function(df, gt, region) {
  sub <- df[!is.na(df$GT) & df$GT == gt, ]
  n <- nrow(sub)

  # サンプル数が0の場合、本当は何も返すべきでないが今は0カバレッジということにした
  if (n == 0) return(list(
    plus  = Rle(0L, width(region)),
    minus = Rle(0L, width(region)),
    n = 0
  ))

  sum_plus  <- Rle(0L, width(region))
  sum_minus <- Rle(0L, width(region))

  for (i in seq_len(n)) {
    print(i)
    cov <- get_stranded_cov(sub$bam[i], region)
    scale <- sub$norm[i]

    sum_plus  <- sum_plus  + cov$plus  / scale
    sum_minus <- sum_minus + cov$minus / scale
  }

  list(
    plus  = sum_plus  / n,
    minus = sum_minus / n,
    n     = n
  )
}

make_cov_gr <- function(cov, region) {
  pos <- start(region):end(region)

  GRanges(
    seqnames = CHROM,
    ranges = IRanges(pos, pos),
    Forward = as.numeric(cov$plus),
    Reverse = -as.numeric(cov$minus)
  )
}

get_cov_max <- function(cov_list, strand = c("both", "+", "-")) {
  strand <- match.arg(strand)

  vals <- unlist(lapply(cov_list, function(cov) {
    if (strand == "+") {
      as.numeric(cov$plus)
    } else if (strand == "-") {
      as.numeric(cov$minus)
    } else {
      c(as.numeric(cov$plus), as.numeric(cov$minus))
    }
  }))

  max(abs(vals), na.rm = TRUE)
}


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

print(table(gt_df$GT, useNA = "ifany"))


# 2
gt_df <- gt_df[, c("sample", "GT")]

bam_map <- read.table(BAM_MAP, header=TRUE, sep="\t")
norm_df <- read.table(NORM_FACTOR_MAP, header=TRUE, sep="\t",
                      col.names=c("sample", "norm"))

sample_df <- Reduce(
  function(x, y) merge(x, y, by="sample"),
  list(gt_df, bam_map, norm_df)
)

sample_df$GT <- factor(sample_df$GT, levels=c("0/0", "0/1", "1/1"))


# 3
region <- GRanges(CHROM, IRanges(START, END))
mean_cov_00 <- get_mean_cov_by_gt(sample_df, "0/0", region)
mean_cov_01 <- get_mean_cov_by_gt(sample_df, "0/1", region)
mean_cov_11 <- get_mean_cov_by_gt(sample_df, "1/1", region)

gr_00 <- make_cov_gr(mean_cov_00, region)
gr_01 <- make_cov_gr(mean_cov_01, region)
gr_11 <- make_cov_gr(mean_cov_11, region)


# 4
# gene model作るor読む
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
  transcriptAnnotation="gene", # "symbol",
  name="Gene Model")


# ylim範囲計算
cov_all <- list(mean_cov_00, mean_cov_01, mean_cov_11)

cov_max <- get_cov_max(cov_all, strand = YLIM_RANGE_STRAND)
cov_max <- ceiling(cov_max * 1.1)

ylim_val <- c(-cov_max, cov_max)


# トラックを作る
tracks <- list(
  GenomeAxisTrack(),
  gene_track,
  DataTrack(
    gr_00,
    name = paste0("0/0 (n=", mean_cov_00$n, ")"),
    type="hist",
    groups=c("Forward","Reverse"),
    col=c("cornflowerblue","purple"),
    legend=FALSE,
    ylim = ylim_val
  ),
  DataTrack(
    gr_01,
    name = paste0("0/1 (n=", mean_cov_01$n, ")"),
    type="hist",
    groups=c("Forward","Reverse"),
    col=c("cornflowerblue","purple"),
    legend=FALSE,
    ylim = ylim_val
  ),
  DataTrack(
    gr_11,
    name = paste0("1/1 (n=", mean_cov_11$n, ")"),
    type="hist",
    groups=c("Forward","Reverse"),
    col=c("cornflowerblue","purple"),
    legend=FALSE,
    ylim = ylim_val
  )
)

# 図を保存
png(OUT, width = WIDTH, height = HEIGHT * 8, res = 200)
plotTracks(
  tracks,
  chromosome = CHROM,
  from = START,
  to = END,
  sizes = c(0.5, 2, 0.5, 0.5, 0.5) # rep(0.5, 3+2)
)
dev.off()

