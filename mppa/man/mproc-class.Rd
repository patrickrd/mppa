\name{mproc-class}
\Rdversion{1.1}
\docType{class}
\alias{mproc-class}

\title{Class \code{mproc}}
\description{
  Manipulation of sets of one-dimensional processes - Internal.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("mproc", ...)},
  or use \code{\link{mp}} (recommended).
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{events}:}{Object of class \code{list}}.
    \item{\code{start}:}{Object of class \code{numeric}}.
    \item{\code{end}:}{Object of class \code{numeric}}.
  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "mproc", y = "missing")}: ... }
	 }
}
\author{
  Patrick Rubin-Delanchy <patrick.rubin-delanchy@bristol.ac.uk>
}

\seealso{
\code{\link{mp}}
}
\examples{
showClass("mproc")
}
\keyword{classes}
