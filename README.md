# CNAPE
**C**opy **N**umber **A**lteration **P**rediction from gene **E**xpression in human cancers


## Ownership

[Wang Lab at HKUST](http://wang-lab.ust.hk/)

## Status
Active development


## Introduction

Copy number alterations (CNAs) are important features of cancer. While the standard methods for CNA detection (CGH arrays, SNP arrrays, DNA sequencing) rely on DNA, occasionally DNA data are not available, especially in cancer studies (e.g. biopsies, legacy data). CNAPE comes into play by predicting CNAs based on gene expression data from RNA-seq.


## How to run
Please clone this repository by
```
git clone https://github.com/Kevin-Moo/CNAPE
```
To test the environment, go to the CNAPE folder and run :
```
./run_example.sh
```

The main function of CNAPE is packaged in *cnape.R*. Get your gene expression profile prepared, and run it like this:

```
Rscript cnape.R expressionMatrix outputPrefix
```
CNAPE.R takes two arguments: the first one is your expression matrix, and the second one is the prefix of the output.
The format of the input gene expression matrix is exemplified in the example. Meanwhile, please make sure that your RNA-seq data is processed using [TCGA's RNA-seq processing pipeline](https://webshare.bioinf.unc.edu/public/mRNAseq_TCGA/UNC_mRNAseq_summary.pdf) (i.e., reads were
aligned to the human genome using [MapSplice](https://academic.oup.com/nar/article/38/18/e178/1068935) and expression was quantified/normalized using [RSEM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323) against UCSC genes).

The output contains *prefix*.chromosome_level.cna.txt and *prefix*.arm_level.cna.txt, where 1 means amplified, -1 means deleted, while 0 means no CNA change.

## Dependencies

The models are trained on the [TCGA Pancancer Atlas data](https://gdc.cancer.gov/about-data/publications/pancanatlas), using [*glmnet*](https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html) package in R. Other dependencies include *reshape2*, *ggplot2* and *pheatmap*, all for visualization. The dependency requirements are automatically solved while running the program.

## Contact
For technical issues please email professor Jiguang Wang: jgwang@ust.hk.
