\name{equivt}
\alias{equivt}

\title{Equivalence Test}

\description{
     Test for the equivalence of two groups with a given tolerable interval. If a significant effect is returned for both t-values then there is evidence of equivalence. Allows for equal or unequal variances.}

\usage{equivt(x, y, equivint, alpha = 0.05, varequiv = FALSE, na.rm = TRUE, ...)}

\arguments{
    \item{x}{Numeric vector for group 1.}
    \item{y}{Numeric vector for group 2.}
    \item{equivint}{Equavalence interval.}
    \item{alpha}{Alpha value}
    \item{varequiv}{Logical. Are there equal variances? Default is \code{FALSE}.}
    \item{na.rm}{Logical. Remove missing data?}
    \item{\dots}{Additional arguments to be passed to the function.}
}

\details{
    \code{equivint} must be specified in the metric of the variables (i.e. not as a standard deviation, percent different, etc).}
    
\author{Rob Cribbie}

\examples{
data(nonnorm_hetvar)
attach(nonnorm_hetvar)

#equivilence within 2 points
eq <- 2
equivt(depres[Sex=="female"],depres[Sex=="male"], eq)
}

