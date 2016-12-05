\name{simes.test}
\alias{simes.test}

\title{
Simes's method for combining p-values
}
\description{
  Computes the combined p-value \eqn{\min(n p_{(i)}/i)} where
  \eqn{p_{(1)}, \ldots, p_{(n)}} are the ordered p-values.
}
\usage{
simes.test(x, returnstat = FALSE)
}

\arguments{
  \item{x}{
    a vector of p-values.
}
  \item{returnstat}{
    if TRUE then return the statistic as well as the p-value. (Very
    unusually, they are the same with this method.)
}
}
\value{
  if returnstat=TRUE then a pair (p,T) where p is the p-value and T is
  the test, otherwise just p.
}
\references{
Simes, R. J. (1986). An improved Bonferroni procedure for multiple tests of significance. Biometrika, 73(3):751-754.
}
\author{
Patrick Rubin-Delanchy <patrick.rubin-delanchy@bristol.ac.uk>
}

\seealso{
\code{\link{F.test}}, \code{\link{TMT.test}}
}
\examples{
## This should be uniformly distributed
simes.test(runif(10))
}
\concept{Combining p-values}
