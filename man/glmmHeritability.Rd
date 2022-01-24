\name{glmmHeritability}
\alias{glmmHeritability}
\title{Heritability estimation}
\description{
    Get the heritability estimate from the SAIGE model.
}

\usage{
glmmHeritability(modobj, adjust=TRUE)
}

\arguments{
    \item{modobj}{an R object for SAIGE model parameters}
    \item{adjust}{if \code{TRUE} and binary outcomes, uses adjusted tau
        estimate for the heritability estimation}
}

\details{
    In SAIGE, penalized quasi-likelihood (PQL) is used to estimate the variance
component parameter tau. It is known to produce biased estimate of the variance
component tau using PQL. If \code{adjust=TRUE} for binary outcomes, tau is
adjusted based prevalence and observed tau using the data in Supplementary
Table 7 (Zhou et al. 2018) to reduce the bias of PQL estimate of variance
component.
}

\value{
    Return a liability scale heritability.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}
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

glmmHeritability(glmm)

seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
