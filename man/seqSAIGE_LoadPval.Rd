\name{seqSAIGE_LoadPval}
\alias{seqSAIGE_LoadPval}
\title{Load the association results}
\description{
    Load the association results from an RData or GDS file.
}

\usage{
seqSAIGE_LoadPval(fn, varnm=NULL, index=NULL, verbose=TRUE)
}

\arguments{
    \item{fn}{RData or GDS file names, merging datasets if multiple files}
    \item{varnm}{NULL, or a character vector to include the column names; e.g.,
        \code{c("chr", "position", "rs.id", "ref", "alt", "pval")}}
    \item{index}{NULL, or a logical/numeric vector for a set of rows}
    \item{verbose}{if \code{TRUE}, show information}
}

\value{
    Return a \code{data.frame} including p-values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}, \code{\link{seqAssocGLMM_SPA}}
}

\examples{
(fn <- system.file("unitTests", "saige_pval.rda", package="SAIGEgds"))
pval <- seqSAIGE_LoadPval(fn)

names(pval)
#  [1] "id"         "chr"        "pos"        "rs.id"      "ref"
#  [6] "alt"        "AF.alt"     "AC.alt"     "num"        "beta"
# [11] "SE"         "pval"       "pval.noadj" "converged"

head(pval)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
