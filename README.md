SAIGEgds: Scalable Implementation of Generalized mixed models using GDS files
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Features

Scalable and accurate implementation of generalized mixed mode with the support of Genomic Data Structure ([GDS](https://github.com/zhengxwen/SeqArray)) files and highly optimized C++ implementation. It is designed for single variant tests in large-scale phenome-wide association studies (PheWAS) with millions of variants and hundreds of thousands of samples, e.g., [UK Biobank genotype data](https://www.ukbiobank.ac.uk).

The implementation of SAIGEgds is based on the original [SAIGE](https://github.com/weizhouUMICH/SAIGE) R package (v0.29.4.4). It is implemented with optimized C++ codes taking advantage of sparse structure of genotypes. All of the calculation with single-precision floating-point numbers in [SAIGE](https://github.com/weizhouUMICH/SAIGE) are replaced by the double-precision calculation in SAIGEgds. SAIGEgds also implements some of the [SPAtest](https://cran.r-project.org/web/packages/SPAtest/index.html) functions in C to speed up the calculation of Saddlepoint Approximation.

Benchmarks using the UKBiobank White British genotype data (N=430K) with coronary heart disease and simulated cases, show that SAIGEgds is 5 to 6 times faster than the SAIGE R package in the steps of fitting null models and p-value calculations. When used in conjunction with high-performance computing (HPC) clusters and/or cloud resources, SAIGEgds provides an efficient analysis pipeline for biobank-scale PheWAS.


## Package Maintainer

Dr. [Xiuwen Zheng](xiuwen.zheng@abbvie.com)


## Installation

* Requires R (≥ v3.5.0), [SeqArray](http://www.bioconductor.org/packages/SeqArray) (≥ v1.24.1)

* Bioconductor repository (available soon)
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SAIGEgds")
```

* Development version from Github
```R
library("devtools")
install_github("AbbVie-ComputationalGenomics/SAIGEgds")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.

* Package vignette

If the package is installed from Bioconductor repository, users can start R and enter to view documentation:
```R
browseVignettes("SAIGEgds")
```


## Examples

```R
library(SeqArray)
library(SAIGEgds)

# open the GDS file for genetic relationship matrix
grm_fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
(grm_gds <- seqOpen(grm_fn))

# load phenotype
phenofn <- system.file("extdata/pheno.txt.gz", package="SAIGEgds")
pheno <- read.table(phenofn, head=TRUE, as.is=TRUE)
head(pheno)
##   y      x1 x2 sample.id
## 1 0  1.5118  1        s1
## 2 0  0.3898  1        s2
## 3 0 -0.6212  1        s3
## ...

# fit the null model
glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, grm_gds, num.thread=2)
## SAIGE association analysis:
## Filtering variants:
## [==================================================] 100%, completed (0s)
## Fit the null model: y ~ x1 + x2 + var(GRM)
##     # of samples: 1,000
##     # of variants: 9,976
##     using 2 threads
## ...

# close the file
seqClose(grm_gds)



################################

# open the GDS file for association testing
geno_fn <- system.file("extdata/assoc_100snp.gds", package="SAIGEgds")
(geno_gds <- seqOpen(geno_fn))

# p-value calculation
assoc <- seqAssocGLMM_SPA(geno_gds, glmm, mac=10, parallel=2)
## SAIGE association analysis:
##     # of samples: 1,000
##     # of variants: 100
##     p-value threshold for SPA adjustment: 0.05
##     variance ratio for approximation: 0.9410486
##     # of processes: 2
## [==================================================] 100%, completed, 0s
## # of variants after filtering MAF/MAC: 38
## Done.

head(assoc)
##   id chr pos rs.id ref alt AF.alt AC.alt  num       beta       SE      pval pval.noadj converged
## 1  4   1   4   rs4   A   C 0.0100     20 1000 -0.0748377 0.790873 0.9246113  0.9246113      TRUE
## 2 12   1  12  rs12   A   C 0.0150     30 1000 -0.0908154 0.656467 0.8899720  0.8899720      TRUE
## 3 14   1  14  rs14   A   C 0.0375     75 1000 -0.0753012 0.433707 0.8621624  0.8621624      TRUE
## ...

# close the file
seqClose(geno_gds)
```


## Citations

Zheng X, Davis J.Wade, SAIGEgds -- an efficient statistical tool for large-scale PheWAS with mixed models; (Abstract 1920276). Presented at the 69th Annual Meeting of The American Society of Human Genetics (ASHG), Oct 15-19, Houston, US.

Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN, LeFaive J, VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ, Denny JC, Lin M, Hveem K, Kang HM, Abecasis GR, Willer CJ, Lee S. Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. *Nat Genet* (2018). Sep;50(9):1335-1341. [DOI: 10.1038/s41588-018-0184-y](https://www.nature.com/articles/s41588-018-0184-y).

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D. SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics* (2017). [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).


## Also See

[SeqArray](https://www.bioconductor.org/packages/SeqArray): Data management of large-scale whole-genome sequence variant calls

[gds2bgen](https://github.com/zhengxwen/gds2bgen): Format conversion from bgen to gds

