\name{seqAssocGLMM_SPA}
\alias{seqAssocGLMM_SPA}
\title{P-value calculation}
\description{
    P-value calculations using variance approximation and an adjustment of
Saddlepoint approximation.
}

\usage{
seqAssocGLMM_SPA(gdsfile, modobj, maf=NaN, mac=10, missing=0.1, dsnode="",
    spa.pval=0.05, var.ratio=NaN, res.savefn="", res.compress="LZMA",
    parallel=FALSE, verbose=TRUE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename, or a GDS object}
    \item{modobj}{an R object for SAIGE model parameters}
    \item{maf}{minor allele frequency threshold (checking >= maf), \code{NaN}
        for no filter}
    \item{mac}{minor allele count threshold (checking >= mac), \code{NaN}
        for no filter}
    \item{missing}{missing threshold for variants (checking <= missing),
        \code{NaN} for no filter}
    \item{dsnode}{"" for automatically searching the GDS nodes "genotype" and
        "annotation/format/DS", or use a user-defined GDS node in the file}
    \item{spa.pval}{the p-value threshold for SPA adjustment, 0.05 by default
        (since normal approximation performs well when the test statistic is
        close to the mean)}
    \item{var.ratio}{\code{NaN} for using the estimated variance ratio in the
        model fitting, or a user-defined variance ratio}
    \item{res.savefn}{an RData or GDS file name, "" for no saving}
    \item{res.compress}{the compression method for the output file, it should
        be one of LZMA, LZMA_RA, ZIP, ZIP_RA and none}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), a numeric value for the number of cores, or other value;
        \code{parallel} is passed to the argument \code{cl} in
        \code{\link{seqParallel}}, see \code{\link{seqParallel}} for more
        details}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    The original SAIGE R package uses 0.05 as a threshold for unadjusted
p-values (based on asymptotic normality) to further calculate adjusted
p-values (Saddlepoint approximation, SPA). If \code{var.ratio=NaN},
the average of variance ratios (\code{mean(modobj$var.ratio$ratio)}) is used
instead.
    For more details of SAIGE algorithm, please refer to the SAIGE paper
[Zhou et al. 2018] (see the reference section).
}

\value{
    Return a \code{data.frame} with the following components if not saving to
a file:
    \item{id}{variant ID in the GDS file;}
    \item{chr}{chromosome;}
    \item{pos}{position;}
    \item{rs.id}{the RS IDs if it is available in the GDS file;}
    \item{ref}{the reference allele;}
    \item{alt}{the alternative allele;}
    \item{AF.alt}{allele frequency for the alternative allele; the minor allele
        frequency is \code{pmin(AF.alt, 1-AF.alt)}}
    \item{mac}{minor allele count; the allele count for the alternative allele
        is \code{ifelse(AF.alt<=0.5, mac, 2*num-mac)}}
    \item{num}{the number of samples with non-missing genotypes;}
    \item{beta}{beta coefficient, odds ratio if binary outcomes (alternative
        allele vs. reference allele);}
    \item{SE}{standard error for beta coefficient;}
    \item{pval}{adjusted p-value with SPA;}
    \item{pval.noadj}{p-value without adjustment (it is used for model
        checking only);}
    \item{converged}{whether the SPA algorithm converges or not for adjusted
        p-values.}
}

\references{
    Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN,
LeFaive J, VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ,
Denny JC, Lin M, Hveem K, Kang HM, Abecasis GR, Willer CJ, Lee S.
Efficiently controlling for case-control imbalance and sample relatedness in
large-scale genetic association studies. Nat Genet (2018). Sep;50(9):1335-1341.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAssocGLMM_SPA}}, \code{\link{seqSAIGE_LoadPval}}
}

\examples{
# open a GDS file
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
