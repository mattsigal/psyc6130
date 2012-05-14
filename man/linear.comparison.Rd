\name{linear.comparison}
\alias{linear.comparison}

\title{Linear comparisons}

\description{
    Computes an F test for a linear comparison among group means. Uses procedures described in Maxwell & Delaney (2004; Designing Experiments and Analyzing Data), pages 1157-1162. This routine is suitable for one-way designs.
}

\usage{
linear.comparison(y, group, c.weights, alpha = 0.05,
    var.equal = TRUE, print.output = TRUE)
}

\arguments{
  \item{y}{Dependent variable}
  \item{group}{Grouping variable}
  \item{c.weights}{A list containing weights for multiple linear contrasts; or a numeric vector containing weights for 1 contrast}
  \item{alpha}{compute 100*(1-alpha)\% confidence interval; default is 0.05}
  \item{var.equal}{Logical. Are the variances equal? Default is \code{TRUE}.}  
}

\details{
    If variances are assumed to differ across groups, the routine calculates MS-within and adjusts df-within using procedures developed by Welch (1938) & Shatterwaite (1946), as described in Maxwell & Delaney (2004; Designing Experiments and Analyzing Data), pages 165-168.

}

\author{Rob Cribbie}

\seealso{
    \code{\link{lc.var.equal}}, \code{\link{lc.var.unequal}}}

\examples{
data(contrasts)
attach(contrasts)

linear.comparison(achiev,major,c.weights=list(c(1,0,-.5,-.5,0),
    c(1,-2,0,0,1),c(1,0,-1,-1,1),c(0,1,1,-1,-1)),var.equal=F)
}
