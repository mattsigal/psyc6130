\name{paired_resample}
\alias{paired_resample}
\title{Paired resampling}
\usage{
  paired_resample(level1, level2, ntails = 2,
    nresamp = 10000)
}
\arguments{
  \item{level1}{numeric vector for the first level}

  \item{level2}{numeric vector for the second level}

  \item{ntails}{option for one or two tailed, default is 2}

  \item{nresamp}{number of samples to draw, default is
  10000}
}
\description{
  Run a paired resampling scheme to obtain one- or
  two-tailed p-values.
}
\examples{
\dontrun{
x<-rnorm(10,mean=-.2)
y<-rnorm(10)
paired_resample(x,y,ntails=1,nresamp=10)
}
}
\author{
  Rob Cribbie and Phil Chalmers
}

