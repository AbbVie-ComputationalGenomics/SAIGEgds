\name{seqAssocGLMM_spaBurden}
\alias{seqAssocGLMM_spaBurden}
\title{Burden tests}
\description{
    Burden p-value calculations using mixed models and the
Saddlepoint approximation method for case-control imbalance.
}

\usage{
seqAssocGLMM_spaBurden(gdsfile, modobj, units, wbeta=AggrParamBeta,
    summac=3, dsnode="", spa.pval=0.05, var.ratio=NaN, res.savefn="",
    res.compress="LZMA", parallel=FALSE, verbose=TRUE, verbose.maf=TRUE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename, or a GDS object}
    \item{modobj}{an R object for SAIGE model parameters}
    \item{units}{a list of units of selected variants, with S3 class
        \code{"SeqUnitListClass"} defined in the SeqArray package}
    \item{wbeta}{weights for per-variant effect, using beta distribution
        \code{dbeta()} according to variant's MAF; a length-two vector, or a
        matrix with two rows for multiple beta parameters; by default, using
        beta(1,1) and beta(1,25) both}
    \item{summac}{a threshold for the weighted sum of minor allele counts
        (checking \code{>= summac})}
    \item{dsnode}{"" for automatically searching the GDS nodes "genotype" and
        "annotation/format/DS", or use a user-defined GDS node in the file}
    \item{spa.pval}{the p-value threshold for SPA adjustment, 0.05 by default}
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
    \item{verbose.maf}{if \code{TRUE}, show summary of MAFs in units}
}

\details{
    The original SAIGE R package uses 0.05 as a threshold for unadjusted
p-values to further calculate SPA-adjusted p-values. If \code{var.ratio=NaN},
the average of variance ratios (\code{mean(modobj$var.ratio$ratio)}) is used
instead.
    For more details of SAIGE algorithm, please refer to the SAIGE paper
[Zhou et al. 2018] (see the reference section).
}

\value{
    Return a \code{data.frame} with the following components if not saving to
a file:
    \code{chr}, chromosome;
    \code{start}, a starting position;
    \code{end}, an ending position;
    \code{numvar}, the number of variants in a window;
    \code{summac}, the weighted sum of minor allele counts;
    \code{beta}, beta coefficient, odds ratio if binary outcomes);
    \code{SE}, standard error for beta coefficient;
    \code{pval}, adjusted p-value with Saddlepoint approximation;
    \item{p.norm}{p-values based on asymptotic normality (could be 0 if it
        is too small, e.g., \code{pnorm(-50) = 0} in R; used for checking only}
    \code{cvg}, whether the SPA algorithm converges or not for adjusted p-value.
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
    \code{\link{seqAssocGLMM_spaACAT_V}}, \code{\link{seqAssocGLMM_spaACAT_O}}
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

# get a list of variant units for burden tests
units <- seqUnitSlidingWindows(gdsfile, win.size=500, win.shift=250)

assoc <- seqAssocGLMM_spaBurden(gdsfile, glmm, units)
head(assoc)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
