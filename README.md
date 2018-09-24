# CNAPE
__C__opy __N__umber __A__lteration __P__rediction from gene __E__xpression in human cancers


## Ownership

[Wang Lab at HKUST](http://wang-lab.ust.hk/)

## Status
Active development


## Introduction

Copy number alterations (CNAs) are important features of cancer. While the standard methods for CNA detection (CGH arrays, SNP arrrays, DNA sequencing) rely on DNA, occasionally DNA data are not available, especially in cancer studies (e.g. biopsies). CNAPE comes into play by predicting CNAs based on gene expression data from RNA-seq.


## How to run

An example run is shown in a shell script:

```
user@linux$ ./run_example.sh
```
You can modify the script to run on your own data.

Meanwhile, please make sure that your RNA-seq data is processed using [TCGA's RNA-seq processing pipeline](https://webshare.bioinf.unc.edu/public/mRNAseq_TCGA/UNC_mRNAseq_summary.pdf) (i.e., reads were
aligned to the human genome using [MapSplice](https://academic.oup.com/nar/article/38/18/e178/1068935) and expression was quantified/normalized using [RSEM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323) against UCSC genes).


## Dependencies

The models are trained on the [TCGA Pancancer Atlas data](https://gdc.cancer.gov/about-data/publications/pancanatlas), using [glmnet](https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html) package in R.

## Contact
For technical issues please contact Kevin via email: qmu@connect.ust.hk
