\name{TMT.test}
\alias{TMT.test}
\title{
A method for combining p-values based on a changepoint model
}
\description{
  Under the null hypothesis that  \eqn{p_1, \ldots, p_n} are independent
  and uniformly distributed on \eqn{[0,1]}, this function computes the
  p-value for the maximum likelihood of the changepoint model \eqn{d(p)
  \propto \lambda_1} if \eqn{p \le \tau} and \eqn{d(p) \propto \lambda_2} otherwise,  where \eqn{d} is a piecewise constant density with one changepoint, with unknown parameters \eqn{\lambda_1 \ge
    \lambda_2} and \eqn{\tau \in [0,1]}. The null hypothesis is
  \eqn{\lambda_1 = \lambda_2}.
}
\usage{
TMT.test(x, method = "AUTO", maxtau = 1, samples = 1000)
}

\arguments{
  \item{x}{
    a vector of p-values.
}
  \item{method}{
    if ``AUTO'' then the computation method will be chosen
  automatically. If ``NOE'', the p-value is computed exactly (effort is order
  \eqn{n^2}) by Noe's recursion. If "SIM", then the p-value is
  determined by simulation.
}
  \item{maxtau}{
    restrict the changepoint to a maximum maxtau. Useful if we are only
    interested in testing small p-values.  
}
  \item{samples}{
    number of samples used for simulation. Ignored if method is ``NOE''
  or method is ``AUTO'' with \eqn{n < 1000}.
}
}

\details{
  The p-value for the maximum likelihood can be restated as the solution
to
\deqn{P[u_1 \ge o_1, \ldots, u_n \ge o_n],}
where \eqn{o_1, \ldots, o_n} are a sequence determined from the maximum
  likelihood and \eqn{u_1, \ldots, u_n} are ordered uniform random
  variables. Computing this probability is harder than it looks, because
  simple analytical recursions fail due to catastrophic
  cancellation. We have implemented Noe's recursions, which are safe but
  expensive for large \eqn{n}, so we recommend simulation for \eqn{n
  \ge 1000}. This will be done automatically if method=``AUTO''.
}
\value{
  \item{lk}{the maximum likelihood.}
  \item{TMTi}{the rank of the largest p-value to the left of the
    estimated changepoint (the rank of the largest p-value in the group of
  `small' p-values).}
  \item{TMTu}{the largest p-value to the left of the changepoint (the
    largest p-value in the group of `small' p-values).}
  \item{p}{the p-value of the test. If maxtau was set, and no inputs were
  lower than maxtau, the p-value returned is 1.}
}
\references{
  Patrick Rubin-Delanchy and Nicholas A Heard. ``A test for dependence
  between two point processes on the real line''.  arXiv:1408.3845.
  
Noe, M. (1972). The calculation of distributions of two-sided
  Kolmogorov-Smirnov type statistics. The Annals of Mathematical
  Statistics, pages 58-64.
  
  Noe, M. and Vandewiele, G. (1968). The calculation of distributions of Kolmogorov-Smirnov type statistics including a table of significance points for a particular case. The Annals of Mathematical Statistics, 39(1):233-241.
}
\author{
Patrick Rubin-Delanchy <patrick.rubin-delanchy@bristol.ac.uk>
}

\seealso{
\code{\link{F.test}}, \code{\link{simes.test}}
}
\examples{
## This should be uniformly distributed
TMT.test(runif(10))$p
## Whenever no p-values fall inside [0,maxtau] the returned p-value is one
replicate(20,TMT.test(runif(10), maxtau=.1)$p)
##Use maxtau to gain extra detection power if only interested in a
##subset of very low p-values, e.g.
TMT.test(c(.04, .5))$p
##is larger than
TMT.test(c(.04, .5), maxtau=0.05)$p
}
\concept{Combining p-values}


