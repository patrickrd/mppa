\name{corrtest}
\alias{corrtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A test for dependence between two point processes on the real line.
}
\description{
  Given two point processes, \eqn{A} and \eqn{B},  simultaneously
  observed on the real line, this function computes a p-value for
  \eqn{B}'s dependence on \eqn{A}. The function can be parameterised to
  detect triggering (`\eqn{A} causes \eqn{B}'), correlation, inhibition
  or anti-correlation. Further details are available in the corresponding paper.
}
\usage{
corrtest(A, B, start = 0, end = 1, F = NULL, alternative = "causes",
      method = "timeout", transform = FALSE, usebeforeA1 = FALSE,
             maxtau = NA, careful = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
    a vector of event times. We condition on A and test B given A.
}
  \item{B}{
    a second vector of event times. B is tested conditional on A.
}
  \item{start}{
    the start of the observation period. Must be smaller than all of the
    event times in A and in B.
}
  \item{end}{
    the end of the observation period. Must be larger than all of the
    event times in A and in B.
}
  \item{F}{
    an increasing function, giving the cumulative intensity of B under
    the null hypothesis. If \code{F=NULL}, the
    cumulative intensity is assumed to be \eqn{F(t)=t} (a homogeneous Poisson process). See details.
}
  \item{alternative}{
    can be "causes" (testing for A causing B), "inhibits"
    (testing for A inhibiting B), "correlation" or "anticorrelation". 
}
  \item{method}{
    can be "timeout", "simes" or "fisher". These
    are different ways of testing the cumulative response times. See details.
}
  \item{transform}{
    if \code{TRUE}, \code{F} is first used to `homogenise' A and B by time
    change. See details.
}
  \item{usebeforeA1}{
    if \code{TRUE}, the event times before \code{A[1]} are assumed to follow the null model and
    included in the test. Ignored if alternative is "correlation" or "anticorrelation".
}
  \item{maxtau}{
    if method="timeout" then the user can set \code{maxtau} to limit the range of
    dependence being tested for. For example, if alternative is "causes"
    and \code{maxtau} is .1, then
    only events in B within .1 to the right of an event in A
    can be triggered.
}
  \item{careful}{
    if \code{TRUE} a number of error checks are made before starting the procedure.
}
}
\details{
Under the null hypothesis the events of \eqn{B} are assumed to follow a
non-homogeneous Poisson process. This has an intensity function
\eqn{\lambda(t)}. The intensity is provided by the user via the input
\code{F}. If \code{F} is \code{NULL} (the default), \eqn{B} follows a
homogeneous Poisson process under the null hypothesis. Otherwise
\code{F} is interpreted as \eqn{F(t) = \int_{\mbox{start}}^t
  \lambda(x) dx}. The tests are invariant to any affine transformation of the input function
\eqn{F^*(t) = \alpha+\beta F(t)}. Note that \code{F=NULL} is
equivalent to \code{F=function(t)t}, although the latter is not
recommended because the algorithm is faster if it knows the null
hypothesis is homogeneous. To choose \code{F} the users are invited to
estimate the non-homogeneous intensity of \eqn{B} via their favourite method, or
try \code{mkF} which uses R's in-built \code{density} estimation. Note
that the non-homogeneous Poisson assumption is a bit stronger than needed, see
corresponding paper.

The procedure computes so-called cumulative response times. These are a
sequence of values between 0 and 1, one for each event time in \eqn{B}  (or each B after
\code{A[1]} in the case of causes/inhibition with
\code{usebeforeA1=FALSE}), which should be ordered uniform variables
under the null hypothesis but `closer to zero' under the
alternative.

To test the cumulative response times, the user has the choice between three methods: a timeout
test (see \code{\link{TMT.test}}), Fisher's method (see
\code{\link{F.test}}) and Simes's test (see \code{\link{simes.test}}).

By default \eqn{A} and \eqn{B} are analysed as given. If
\code{transform=TRUE}, \eqn{F} is first used to transform \eqn{A} and
\eqn{B} so that \eqn{B} is homogeneous, and then the transformed processes are tested
against a homogeneous null hypothesis.

For further details see corresponding paper.
}
\value{
An S3 object of class ``htest'' with slots
\item{p}{the p-value.}
\item{T}{the test statistic.}
\item{data.name}{the input processes.}
\item{t}{the cumulative response times.}
\item{TMTval}{if method="timeout", the largest cumulative response time
  estimated to have been affected by A, \code{NA} otherwise.}
\item{TMTi}{if method="timeout", the index of the largest cumulative
  response time estimated to have been affected by \code{A}, \code{NA} otherwise.}
}
\references{
Patrick Rubin-Delanchy and Nicholas A Heard. ``A test for dependence between two point processes on the real line''.  arXiv:1408.3845.
}
\author{
Patrick Rubin-Delanchy <patrick.rubin-delanchy@bristol.ac.uk> and
Nicholas A Heard <n.heard@imperial.ac.uk>
}


\seealso{
\code{\link{mkF}}, \code{\link{TMT.test}}, \code{\link{F.test}}, \code{\link{simes.test}}
}
\examples{
A = runif(20)
Braw=runif(20)
##around ten B events are caused by A 
B1=c(Braw, sample(A, 10)+abs(rnorm(10, 0,.01))); B1 = B1[B1>0&B1<1]
##about ten B events are correlated to A
B2=c(Braw, sample(A, 10)+rnorm(10, 0,.01)); B2 = B2[B2>0&B2<1]
##the ten closest B events to A are deleted (anticorrelation)
d=sapply(Braw, function(b) min(abs(b-A))); B3=Braw[-order(d)[1:10]]
##The ten B events that are closest to but after an A event are deleted (inhibition).
##(Adding an A event at zero to be sure there are 10.)
A=c(0,A); d=sapply(Braw, function(b) b-max(A[A<b])); B4=Braw[-order(d)[1:10]]

alternatives=c("causes", "correlation", "anticorrelation", "inhibits")
ps=c()
for (B in list(B1, B2, B3, B4)){
    for (alternative in alternatives){
        p=corrtest(A,B, alternative=alternative)$p
        ps = c(ps, p)
    }
}
M=matrix(ps, nrow=4, ncol=4)
colnames(M) = c("causal_data", "correlated_data", "anticorrelated_data", "inhibited_data")
rownames(M) = c("causes_test", "correlation_test", "anticorrelation_test", "inhibition_test")
##should (hopefully!) see low p-values on the diagonal:
M
## and high p-values for opposite data versus test combinations,
## e.g. testing for A inhibiting B when in fact A causes B:
M["inhibition_test", "causal_data"]

##Now for a non-homogeneous example. A and B have a common daily pattern:
##their intensity is a sinusoidal curve lambda(t) = 1+sin(2*pi*t)
start=0; end=365 #A year
##the cumulative intensity is 
F=function(t){
    t\%/\%1+t\%\%1+(1-cos(2*pi*t\%\%1))/(2*pi)
}
##Dropping 365 A and B points according to F
A=sapply(runif(365), function(u){uniroot(function(x) F(x)/365-u, interval=c(0,365))$root})
##B is independent of A aside from common periodicity
B=sapply(runif(365), function(u){uniroot(function(x) F(x)/365-u, interval=c(0,365))$root})
##Bc also has some caused events
Bc=c(B, sample(A, 10)+abs(rnorm(10, 0,.001))); Bc = Bc[Bc>start&Bc<end]

##If we do not account for the common periodicity of A and B we get
##spuriously strong evidence for A causing B (using Fisher's method for speed):
corrtest(A,B, start=start, end=end, method="fisher")
##On the other hand with the correct F, the p-value is uniformly distributed:
corrtest(A,B, start=start, end=end, F=F, method="fisher")

##For reference, with the truly caused process Bc we get
corrtest(A,Bc, start=start, end=end, method="fisher")
corrtest(A,Bc, start=start, end=end, F=F, method="fisher")
}
\concept{Testing for dependence between point processes}
