\name{equivt}
\alias{equivt}
\alias{print.equivt}
\title{Equivalence Test}
\usage{
  equivt(x, y, equivint, alpha = 0.05, varequiv = FALSE,
    na.rm = TRUE, ...)

  \method{print}{equivt} (x)
}
\arguments{
  \item{x}{Numeric vector for group 1.}

  \item{y}{Numeric vector for group 2.}

  \item{equivint}{Equivalence interval.}

  \item{alpha}{Alpha value}

  \item{varequiv}{Logical. Are there equal variances?
  Default is \code{FALSE}.}

  \item{na.rm}{logical; remove missing data?}

  \item{...}{Additional arguments to be passed to the
  function.}
}
\description{
  Test for the equivalence of two groups with a given
  tolerable interval. If a significant effect is returned
  for both t-values then there is evidence of equivalence.
  Allows for equal or unequal variances.
}
\section{\code{equivint} must be specified in the metric of the variables (i.e. not as a standard
deviation, percent different, etc).}{

}
\examples{
\dontrun{
data(nonnorm_hetvar)
attach(nonnorm_hetvar)

#equivalence within 2 points
eq <- 2
equivt(depres[Sex=="female"],depres[Sex=="male"], eq)
}
}
\author{
  Rob Cribbie and Phil Chalmers
}

