\name{mcpqest}
\alias{mcpqest}

\title{
Point and variance estimation for quantiles of independent groups of samples
}
\description{
Computes the pth quantile and variances for groups of given samples in one-way anova layout. It has option for right censored data.
}
\usage{
mcpqest(y, f, event = NULL, 
Right.Censored = FALSE, p = 0.5, ...)
}

\arguments{
  \item{y}{
a numeric vector, the response variable. If Right.Censored = True, y is non-negative follow up time for
right censored in survival data.
}
  \item{f}{
a factor variable of the same length as y, assigning the observations in y into k groups.
}
  \item{event}{
a binary variable indicating status for right censored data. Usually, 1 if event of interest has occurred
(death = 1) and 0 otherwise (alive = 0).
}
  \item{Right.Censored}{
a logical expression indicating right-censored data is being used for constructing simultaneous
confidence interval.
}
  \item{p}{
a single numeric value between 0 and 1 indicating the level of quantile for the contrasts. The default
is p = 0.5 (the median).
}
  \item{\dots}{
further arguments to be passed to the internal methods, in particular: \code{bw.selec} is a single character
string specifying the method of bandwidth selection when using right censored survival data; 
\code{bw.selec= "plug-in"}.
}
}
\details{
Mainly for internal use.
}
\value{
a list with elements:
\item{quantileEST }{a numeric vector, the point estimates of quantiles for each factor level.}
\item{varEST }{a numeric vector, the variance estimates for each factor level.}
\item{n}{a numeric vector, the sample size of each factor level}
}

\author{
Lawrence S. Segbehoe, Gemechis Dilba Djira,  Frank Schaarschmidt (package inclusion)
}

\keyword{ htest }

