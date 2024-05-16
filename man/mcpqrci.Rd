\name{mcpqrci}
\alias{mcpqrci}
\alias{mcpqdci}

\title{
Simultaneous confidence intervals for contrasts of quantiles
}
\description{

The following functions construct simultaneous confidence intervals for multiple constrasts of quantiles (for ratios and differences) in a one-way layout. 
The "mcpqrci" is for ratios and "mcpqdci" is for differences of quantiles. Both functions have also options for right censored data.
}

\usage{
mcpqrci(y, f, event = NULL, Right.Censored = FALSE,
p = 0.5, conf.level = 0.95, type = "Dunnett",
base = 1, Num.cmat = NULL, Den.cmat = NULL,
method = c("Wald", "Fieller"), ...)

mcpqdci(y, f, event = NULL, Right.Censored = FALSE,
p = 0.5, conf.level = 0.95, type = "Dunnett", 
base = 1 , cmat = NULL,...)
}

\arguments{
  \item{y}{
a numeric vector, the response variable. If \code{Right.Censored = TRUE}, y is non-negative follow up time for right censored in survival data.
}
  \item{f}{
a factor variable of the same length as y, assigning the observations in y into k groups.
}
  \item{event}{
a binary variable indicating status for right censored data. Usually, 1 if event of interest has occurred (death = 1) and 0 otherwise (alive = 0); (optional: only if y is survival data).
}
  \item{Right.Censored}{
a logical expression indicating right-censored data is being used for constructing simultaneous confidence intervals, (optional: only if y is survival data).
}
  \item{p}{
a single numeric value between 0 and 1 indicating the level of quantile for the contrasts. The default is p = 0.5 (the median).
}
  \item{conf.level}{
a single numeric value between 0 and 1 indicating the level of confidence interval.
}
  \item{type}{
a single character string, naming a contrast type, see contrMat and contrMatRatio, for the options; this argument is ignored if a contrast matrix is specified in cmat or Num.cmat and Den.cmat.
}
  \item{base}{
a positive integer specifying the control group for the Dunnett contrasts, ignored otherwise. When base is not given the first group in terms of an alphanumeric order is taken as the control group.
}
  \item{cmat}{
(optional) a matrix with numeric entries, containing contrast coefficients defining differences of quantiles in function \code{mcpqdci}; if there are k levels in f, the matrix should have k columns. type is ignored if cmat is specified
}
  \item{Num.cmat}{
(optional) Numerator contrast matrix for ratios of quantiles in function \code{mcpqrci}, where the columns correspond to k groups and rows correspond to m contrasts.
}
  \item{Den.cmat}{
(optional) Denominator contrast matrix for ratios of quantiles in function \code{mcpqrci}, where the columns correspond to k groups and rows correspond to m contrasts. type is ignored if Den.cmat and Den.cmat is specified.
}
  \item{method}{
a single character string, naming the method by which to compute the confidence intervals for ratios of quantiles. Default is "Wald". Note if the calculated lower confidence limit is negative and the ratio cannot be negative, set the lower confidence limit to zero.
}

  \item{\dots}{
further arguments to be passed to the internal methods, in particular: dist must be a single character string invoking the use of multivariate normal quantiles; \code{dist="MVN"} or multivariate normal quantiles; \code{dist="MVT"}. \code{bw.selec} is a single character string specifying the method of bandwidth selection when using right censored survival data; \code{bw.selec = "plug-in"}.
}
}
\details{

The interest is to construct simultaneous confidence intervals for several contrast of quantiles in a one-way layout. An asymptotic approach is used in estimating the variance of estimated quantiles. The \code{mcpqrci} handles ratios of multiple contrasts of quantiles and 
\code{mcpqdci} handles differences of multiple contrast of quantiles.


If \code{event} argument is provided and \code{Right.Censored = TRUE}, the functions computes simultaneous confidence intervals for right censored data in y. 
The type argument defines the type of contrast matrix to use. Users can also define a preferred contrast matrix, cmat.

}
\value{
a list with elements

\item{cmat}{Matrix of contrast used for contrast differences.}
\item{Num.Contrast}{Matrix of contrast used for the numerator of ratios.}
\item{Den.Contrast}{Matrix of contrast used for the denominator of ratios.}
\item{conf.level}{A numeric value, as input.}
\item{estimate}{a column vector, containing the point estimates of the contrasts.}
\item{std.err}{a column vector, containing the standard error of the contrast estimates.}
\item{conf.int}{a Mx2 matrix of confidence bounds, if M comparisons among the K samples are invoked.}

}

\author{
Lawrence S. Segbehoe, Gemechis Dilba Djira, and Frank schaarschmidt (inclusion in the package)
}


\seealso{
\code{sciratio} for simultaneous confidence intervals for ratios of linear combinations of means
}
\examples{

data("DiabeticMice")
response <- DiabeticMice$response
group <- DiabeticMice$group
## Example 1
Num.cmat <- matrix(c(1,1,0,0,0,1,0,0,0),3)
Den.cmat <- matrix(c(0,0,0,1,0,0,0,1,1),3)
mcpqdci(y = response, f = group, cmat = (Num.cmat + -1*Den.cmat))
mcpqdci(y = response, f = group, cmat = (Num.cmat + -1*Den.cmat)[-1,])
mcpqrci(y = response, f = group, Num.cmat = Num.cmat, Den.cmat = Den.cmat )
mcpqrci(y = response, f = group, Num.cmat = Num.cmat[-1,], Den.cmat = Den.cmat[-1,] )


## Example 2 

data("bnct")
mcpqrci(y = bnct$time, f = bnct$trt, event = bnct$death, Right.Censored=TRUE)

## Sampled data:

y <- c(rnorm(20),rnorm(16,3),rnorm(24,7,2))
f <- rep(paste0("group", 1:3), c(20,16, 24))
event <- rbinom(60,1,0.8)

mcpqdci(y=y, f=f, method = "Fieller", base = 3)

mcpqrci(y=abs(y), f=f, event=event, Right.Censored=TRUE,
Num.cmat = cbind(c(1,1), 0*diag(2)),
Den.cmat = cbind(c(0,0), diag(2)))

cmat <- cbind(-c(1,1),diag(2))
mcpqdci(y=y, f=f, method = "Fieller", cmat = cmat)
}

\keyword{ htest }
\keyword{ survival }