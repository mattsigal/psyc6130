\name{yuen}
\alias{yuen}

\title{Yuen's test}

\description{Perform Yuen's test for trimmed means on the data in x and y.}
\usage{yuen(x, y, tr = 0.2, alpha = 0.05)}

\arguments{
    \item{x}{Numeric data vector for group 1.  }
    \item{y}{Numeric data for vector group 2.  }
    \item{tr}{  Amount of data to truncate; default is .2.  }
    \item{alpha}{  Alpha value; default is .05.  }
}

\value{
  Function automatically removes missing data. Returns a \code{list} with
  \item{ci }{Lower and upper confidence intervals for the trimmed mean.}
  \item{p.value }{Yuen's p-value.}
  \item{dif }{Mean difference}
  \item{se }{Standard Error}
  \item{teststat }{The calculated test statistic}
  \item{crit}{The critical value}
  \item{df}{The degrees of freedom}
}

\author{Rand Wilcox and Rob Cribbie}

\examples{
hs_less<-c(8,3,6,6,9,7,11,7,4,8,17,6,21,6,4,10,6,9)
greater_hs<-c(7,12,7,6,7,18,15,11,8,7,8,13,12,11,16,15,30,11,10,7,8,25,3,9,5)

### Run a Welch trimmed t-test
yuen(hs_less,greater_hs)
}


