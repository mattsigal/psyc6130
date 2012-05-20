\name{t2way}
\alias{t2way}
\title{Trimmed means for 2-way ANOVA}
\usage{
  t2way(J, K, data, tr = 0.2, grp = c(1:p), p = J * K,
    MAT = F, lev.col = c(1:2), var.col = 3, pr = T)
}
\arguments{
  \item{J}{number of levels for factor 1}

  \item{K}{number of levels for factor 2}

  \item{data}{contains the raw data stored in list mode, or
  a matrix with columns corresponding to groups.  If stored
  in list mode, data[[1]] contains the data for the first
  level of both factors: level 1,1.  data[[2]] is assumed
  to contain the data for level 1 of the first factor and
  level 2 of the second factor: level 1,2, etc. It is
  assumed that data has length JK, the total number of
  groups being tested.}

  \item{tr}{proportion to be trimmed, default is .2}

  \item{grp}{\code{1:p}}

  \item{p}{is a numberic value J*K}

  \item{MAT}{logical; if \code{TRUE} assumes data are
  stored in matrix with 3 columns with two of the columns
  indicated by the argument \code{lev.col} specifying the
  columns of x containing the values of the levels of the
  two factors.}

  \item{lev.col}{something}

  \item{var.col}{something}

  \item{pr}{logical; print the results?}
}
\description{
  A heteroscedastic two-way ANOVA for trimmed means using a
  generalization of Welch's method.
}
\section{By default,}{

}
\examples{
\dontrun{
data(twowaythreeway)
attach(twowaythreeway)

}
}
\author{
  Rob Cribbie and Phil Chalmers
}

