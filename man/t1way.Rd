\name{t1way}
\alias{print.t1way}
\alias{t1way}
\title{Trimmed means ANOVA}
\usage{
  t1way(x, tr = 0.2, grp = NA)

  \method{print}{t1way} (x)
}
\arguments{
  \item{x}{A \code{list} consisting of the data, where
  \code{Length(x)} must correspond to the total number of
  groups.}

  \item{tr}{Proportion to be trimmed; default = .2.}

  \item{grp}{Vector indicating which subset of the groups
  to be compared. Default is \code{NA}}
}
\value{
  Returns a \code{list} with items \item{TEST}{Test
  statistic.} \item{nu1}{Degrees of freedom 1.}
  \item{nu2}{Degrees of freedom 2.}
  \item{p.value}{Probability value.}
}
\description{
  A heteroscedastic one-way ANOVA for trimmed means using a
  generalization of Welch's method.
}
\section{By default, the null hypothesis is that all groups have a common mean.
To compare a subset of the groups, use \code{grp} to indicate which
groups are to be compared. For example, if you type the
command \code{grp<-c(1,3,4)}, and then execute this function, groups
1, 3, and 4 will be compared with the remaining groups ignored.

Missing values are automatically removed.}{

}
\examples{
\dontrun{
data(nhl_data)
attach(nhl_data)

#make a list containing the variable positions
ptslist<-list(Points[Pos=="Center"],Points[Pos=="Left Wing"],
Points[Pos=="Right Wing"], Points[Pos=="Defense"])
t1way(ptslist)
}
}
\author{
  Rand Wilcox, Rob Cribbie and Phil Chalmers
}

