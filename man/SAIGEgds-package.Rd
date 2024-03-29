\name{SAIGEgds-package}
\alias{SAIGEgds-package}
\alias{SAIGEgds}
\docType{package}
\title{Scalable Implementation of Generalized mixed models in Phenome-Wide
    Association Studies using GDS files}
\description{
    Scalable and accurate implementation of generalized mixed mode with the
support of Genomic Data Structure (GDS) files and highly optimized C++
implementation. It is designed for single variant tests in large-scale
phenome-wide association studies (PheWAS) with millions of variants and hundreds
of thousands of samples, e.g., UK Biobank genotype data, controlling for
case-control imbalance and sample structure in single variant association
studies.

    The implementation of SAIGEgds is based on the original SAIGE R package
(v0.29.4.4) [Zhou et al. 2018]
\url{https://github.com/weizhouUMICH/SAIGE/releases/tag/v0.29.4.4}. All of the
calculation with single-precision floating-point numbers in SAIGE are replaced
by the double-precision calculation in SAIGEgds. SAIGEgds also implements some
of the SPAtest functions in C to speed up the calculation of Saddlepoint
Approximation.
}
\details{
    \tabular{ll}{
        Package: \tab SAIGEgds \cr
        Type: \tab Package \cr
        License: \tab GPL version 3 \cr
    }
}
\references{
    Zheng X, Davis J.Wade. SAIGEgds -- an efficient statistical tool for
large-scale PheWAS with mixed models. *Bioinformatics* (2020).
DOI: 10.1093/bioinformatics/btaa731.

    Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN,
LeFaive J, VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ, Denny JC,
Lin M, Hveem K, Kang HM, Abecasis GR, Willer CJ, Lee S. Efficiently controlling
for case-control imbalance and sample relatedness in large-scale genetic
association studies. *Nat Genet* (2018). Sep;50(9):1335-1341.

    Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C,
Levine D. SeqArray -- A storage-efficient high-performance data format for WGS
variant calls. *Bioinformatics* (2017). DOI: 10.1093/bioinformatics/btx145.
}
\author{
    Xiuwen Zheng \email{xiuwen.zheng@abbvie.com},
    Wei Zhou (the original author of the SAIGE R package,
\url{https://github.com/weizhouUMICH/SAIGE})
}

\examples{
# open the GDS file
fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
gdsfile <- seqOpen(fn)

# load phenotype
phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)
head(pheno)

# fit the null model
glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile, trait.type="binary")

# p-value calculation
assoc <- seqAssocGLMM_SPA(gdsfile, glmm, mac=10)

head(assoc)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
