# gviz

`Gviz` は、`conda`でのインストールが難しいが、以下で成功した。

## ✅ 環境構築

```bash
mamba create -n gviz -c bioconda -c conda-forge \
  bioconductor-gviz \
  bioconductor-genomicfeatures \
  bioconductor-biomart \
  bioconductor-rtracklayer \
  bioconductor-bsgenome \
  bioconductor-rsamtools \
  bioconductor-genomicalignments \
  bioconductor-annotationdbi \
  bioconductor-summarizedexperiment \
  bioconductor-variantannotation \
  bioconductor-ensembldb \
  r-base=4.3
