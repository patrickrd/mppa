\name{mkF}
\alias{mkF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of the cumulative intensity of a process using native \code{density}
}
\description{
  A function that estimates the density of the points
  and then creates a valid cumulative intensity function \code{F} for input to \code{\link[mppa]{corrtest}}.
}
\usage{
mkF(x, start = 0, end = 1, adjust = 1, disallow.zero = TRUE)
}

\arguments{
  \item{x}{
    a vector of event times.
}
  \item{start}{
    the start of the observation period: must be smaller than all elements
  of x.
}
  \item{end}{
    the end of the observation period: must be greater than all elements
  of x.
}
  \item{adjust}{
    parameter passed on to \code{\link{density}}.
}
  \item{disallow.zero}{
    if \code{TRUE} do not allow \code{F} to be constant. (Useful to
  avoid spurious results when using \code{\link[mppa]{corrtest}}.)
}
}
\value{
A non-decreasing function that can serve as input to \code{\link[mppa]{corrtest}}.
}
\references{
Patrick Rubin-Delanchy and Nicholas A Heard. ``A test for dependence between two point processes on the real line''.  arXiv:1408.3845.
}
\author{
Patrick Rubin-Delanchy <patrick.rubin-delanchy@bristol.ac.uk>
}
\seealso{
\code{\link[mppa]{corrtest}}
}
\examples{
start=0; end=1
A=rbeta(100, 1, 10)
B=rbeta(100, 1, 10)
##This will be extreme because A and B are not homogeneous over [0,1]
corrtest(A,B)
##If we use an estimate of F the p-value is less aggressive
corrtest(A,B,F=mkF(c(A,B), start=0, end=1))
##But we can still find evidence of A causing B:
Bc=c(B, sample(A, 10)+abs(rnorm(10, 0,.0001))); Bc = Bc[Bc>start&Bc<end]
corrtest(A,Bc,F=mkF(c(A,Bc), start=0, end=1))
}
\concept{Testing for dependence between point processes}
\concept{Cumulative intensity function}
