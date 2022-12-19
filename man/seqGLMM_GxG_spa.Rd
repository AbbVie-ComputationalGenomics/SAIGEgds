\name{seqGLMM_GxG_spa}
\alias{seqGLMM_GxG_spa}
\title{SNP Interaction Testing}
\description{
    SNP interaction testing with Saddlepoint approximation method in the mixed framework.
}

\usage{
seqGLMM_GxG_spa(formula, data, gds_grm, gds_assoc, snp_pair,
    trait.type=c("binary", "quantitative"), sample.col="sample.id", maf=0.005,
    missing.rate=0.01, max.num.snp=1000000L, variant.id=NULL, inv.norm=TRUE,
    X.transform=TRUE, tol=0.02, maxiter=20L, nrun=30L, tolPCG=1e-5,
    maxiterPCG=500L, tau.init=c(0,0), use_approx_tau=FALSE, glm_threshold=FALSE,
    traceCVcutoff=0.0025, ratioCVcutoff=0.001, geno.sparse=TRUE, num.thread=1L,
    model.savefn="", seed=200L, fork.loading=FALSE, verbose=TRUE,
    verbose.detail=TRUE)
}

\arguments{
    \item{formula}{an object of class \code{formula} (or one that can be
        coerced to that class), e.g., \code{y ~ x1 + x2}, see \code{\link{lm}}}
    \item{data}{a data frame for the formulas}
    \item{gds_grm}{a SeqArray GDS filename, or a GDS object}
    \item{gds_assoc}{a SeqArray GDS filename, a GDS object, or a 0/1/2/NA
        matrix with row names for sample IDs}
    \item{snp_pair}{a \code{data.frame} with the first two columns for the
        variant IDs (i.e., 'variant.id') in \code{gds_assoc}}
    \item{trait.type}{"binary" for binary outcomes, "quantitative" for
        continuous outcomes}
    \item{sample.col}{the column name of sample IDs corresponding to the GDS
        file}
    \item{maf}{minor allele frequency for imported genotypes (checking >= maf),
        if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{missing.rate}{threshold of missing rate (checking <= missing.rate),
        if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{max.num.snp}{the maximum number of SNPs used in the calculation of
        GRM, or -1 for no limit}
    \item{variant.id}{a list of variant IDs, used to construct GRM}
    \item{inv.norm}{if \code{TRUE}, perform inverse normal transformation on
        residuals for quantitative outcomes, see the reference [Sofer, 2019]}
    \item{X.transform}{if \code{TRUE}, perform QR decomposition on the design
        matrix}
    \item{tol}{overall tolerance for model fitting}
    \item{maxiter}{the maximum number of iterations for model fitting}
    \item{nrun}{the number of random vectors in the trace estimation}
    \item{tolPCG}{tolerance of PCG iterations}
    \item{maxiterPCG}{the maximum number of PCG iterations}
    \item{tau.init}{a 2-length numeric vector, the initial values for variance
        components, tau; for binary traits, the first element is always be set
        to 1. if \code{tau.init} is not specified, the second element will be
        0.5 for binary traits}
    \item{use_approx_tau}{if \code{TRUE}, fit the model defined in
        \code{formula} without any SNP markers for the interactions to provide
        the estimated tau value (variance component estimates)}
    \item{glm_threshold}{FALSE, TRUE or a numeric value for p-value threshold;
        if TRUE use 0.01 as a threshold}
    \item{traceCVcutoff}{the threshold for coefficient of variation (CV) for
        the trace estimator, and the number of runs for trace estimation will
        be increased until the CV is below the threshold}
    \item{ratioCVcutoff}{the threshold for coefficient of variation (CV) for
        estimating the variance ratio, and the number of randomly selected
        markers will be increased until the CV is below the threshold}
    \item{geno.sparse}{if \code{TRUE}, store the sparse structure for genotypes;
        otherwise, save genotypes in a 2-bit dense matrix; see details}
    \item{num.thread}{the number of threads}
    \item{model.savefn}{the filename of model output, R data file '.rda',
        '.RData', '.rds', '.txt' or '.csv'}
    \item{seed}{an integer as a seed for random numbers}
    \item{fork.loading}{load genotypes via forking or not; forking processes
        in Unix can reduce loading time of genotypes, but may double the memory
        usage; not applicable on Windows}
    \item{verbose}{if \code{TRUE}, show information}
    \item{verbose.detail}{if \code{TRUE}, show the details for model fitting}
}

\details{
    For more details of SAIGE algorithm, please refer to the SAIGE paper
[Zhou et al. 2018] (see the reference section).
}

\value{
    Return a \code{data.frame} with the following components:
    \item{id1}{variant ID for the first SNP in the GDS file;}
    \item{snp1}{includes chromosome, position, reference & alterative alleles for SNP1;}
    \item{maf1}{minor allele frequency for the first SNP;}
    \item{id2}{variant ID for the second SNP in the GDS file;}
    \item{snp2}{includes chromosome, position, reference & alterative alleles for SNP2;}
    \item{maf2}{minor allele frequency for the second SNP;}
    \item{beta}{beta coefficient, odds ratio if binary outcomes;}
    \item{SE}{standard error for beta coefficient;}
    \item{n_nonzero}{the number of non-zero values in the interaction term;}
    \item{pval}{adjusted p-value with the Saddlepoint approximation method;}
    \item{p.norm}{p-values based on asymptotic normality (could be 0 if it
        is too small, e.g., \code{pnorm(-50) = 0} in R; used for checking only}
    \item{converged}{whether the SPA algorithm converges or not for adjusted
        p-values.}
    \item{p.glm}{glm p-value with SPA calculation}
    \item{p.glm.norm}{glm p-value without SPA calculation}
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}, \code{\link{seqAssocGLMM_SPA}}
}

\examples{
# open the GDS file for genetic relationship matrix (GRM)
grm_fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
(grm_gds <- seqOpen(grm_fn))

# load phenotype
phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)
head(pheno)

# define the SNP pairs
snp_pair <- data.frame(s1=2:3, s2=6:7, note=c("F1", "F2"))

seqGLMM_GxG_spa(y ~ x1 + x2, pheno, grm_gds, grm_fn, snp_pair,
    trait.type="binary", verbose.detail=FALSE)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
\keyword{interaction}
