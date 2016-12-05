\name{F.test}
\alias{F.test}

\title{
  Fisher's method for combining p-values
}
\description{
  Compute the p-value for \eqn{-2 \sum \log p_i}.
}
\usage{
F.test(x, returnstat = FALSE)
}

\arguments{
  \item{x}{
    a vector of p-values.
}
  \item{returnstat}{
    if TRUE then return the statistic as well as the p-value.
}
}
\value{
  if returnstat=TRUE then a pair (p,T) where p is the p-value and T is
  the test, otherwise just p.
}
\references{
  Mosteller, F. and Fisher, R. A. (1948). Questions and answers. The American Statistician, 2(5):pp. 30-31.
}
\author{
Patrick Rubin-Delanchy <patrick.rubin-delanchy@bristol.ac.uk>
}
\seealso{
\code{\link{TMT.test}}, \code{\link{simes.test}}
}
\examples{
## This should be uniformly distributed
F.test(runif(10))
}
\concept{Combining p-values}
