\name{data_files}
\alias{data_files}

\title{Data sets used in the course}

\description{A list of data files used throughout the course.}
\usage{
data(datasetName)}

\arguments{
    \item{contrasts}{Data-set for linear contrasts.}
    \item{nhl_data}{NHL data from 2010.}
    \item{nonnorm_hetvar}{Non-normal data with heterogenous variance.}
    \item{one_way}{One way ANOVA data-set.}
}

\value{
  All data-sets are \code{data.frame} objects.
}

\author{Rob Cribbie}

\examples{
data(one_way)
attach(one_way)
}


