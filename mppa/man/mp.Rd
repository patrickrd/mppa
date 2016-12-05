\name{mp}
\alias{mp}

\title{
  Create objects of class \code{mproc}
}
\description{
  A function to create objects of class \code{mproc}, allowing
  manipulation of multiple simultaneous point processes.
}
\usage{
mp(..., start = numeric(0), end = numeric(0))
}

\arguments{
  \item{\dots}{
    either a list of vectors of event times, or the same vectors given as a
    sequence of arguments.
}
  \item{start}{
    the start of observations, must be smaller than all event times. If not
     specified the start will be set to the smallest event time.
}
  \item{end}{
    the end of observations, must be greater than all event times. If not
    specified the end will be set to the greatest event time.
}
}
\value{
  An object of class \code{mproc} with slots:
  \item{events}{a list of vectors of event times.}
  \item{start}{the start of observations.}
  \item{end}{the end of observations.}
}
\author{
Patrick Rubin-Delanchy <patrick.rubin-delanchy@bristol.ac.uk>
}

\examples{
mp(runif(10))
mp(runif(10), runif(10), start=0, end=1)
mp(list(runif(10), runif(10)), start=0, end=1)
}
\concept{Multiple point processes}
