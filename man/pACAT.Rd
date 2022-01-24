\name{pACAT}
\alias{pACAT}
\alias{pACAT2}
\title{Cauchy Combination Test}
\description{
    P-value calculation from Cauchy combination test.
}

\usage{
pACAT(p, w=NULL)
pACAT2(p, maf, wbeta=c(1,25))
}

\arguments{
    \item{p}{a numeric vector for p-values}
    \item{w}{weight for each p-value}
    \item{maf}{minor allele frequency for each p-value}
    \item{wbeta}{weights for per-variant effect, using beta distribution
        \code{dbeta()} according to variant's MAF}
}

\value{
    Return a single number for the combined p-value.
}

\references{
    Liu Y., Cheng S., Li Z., Morrison A.C., Boerwinkle E., Lin X.;
ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis
in Sequencing Studies. Am J Hum Genetics 104, 410-421 (2019).
}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}, \code{\link{seqAssocGLMM_SPA}}
}

\examples{
p1 <- 10^-4
p2 <- 10^-5
p3 <- 10^-(3:20)
sapply(p3, function(p) pACAT(c(p1, p2, p)))

pACAT2(c(10^-4, 10^-6), c(0.01, 0.005))
}

\keyword{association}
\keyword{Cauchy}
