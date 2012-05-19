\name{tukey}
\alias{tukey}
\title{Tukey and Games-Howell tests}
\usage{
  tukey(data, group, method = c("Tukey", "Games-Howell"))
}
\arguments{
  \item{data}{the data}

  \item{group}{logical; remove missing data?}

  \item{method}{either \code{'Tukey'} or
  \code{'Games-Howell'}}
}
\description{
  Tukey and Games-Howell post- hoc tests.
}
\examples{
\dontrun{
data(nonnorm_hetvar)
attach(nonnorm_hetvar)

tukey(depres, Sex, method = 'Games-Howell')
}
}
\author{
  Rob Cribbie and Phil Chalmers
}

