\name{bnct}
\alias{bnct}
\docType{data}
\title{
Boron neutron capture therapy (BNCT)
}
\description{
Death times (in days) from a study to determine the efficacy of BNCT in treating therapeutically refractory F98 glioma.
}
\usage{data("bnct")}
\format{
  A data frame with 30 observations on the following 3 variables.
  \describe{
    \item{\code{trt}}{a numeric vector: Treatment (1=untreated, 2=radiated, 3=radiated + BPA)}
    \item{\code{time}}{a numeric vector: Death time or on-study time, days}
    \item{\code{death}}{a numeric vector: Death indicator (1=dead, 0=alive)}
  }
}
\details{
A right censored data from a study performed to determine the efficacy of boron neutron capture therapy (BNCT) in treating the therapeutically refractory F98 glioma, using boronophenylalanine (BPA) as the capture agent. F98 glioma cells were implanted into the brains of rats. Three groups of rats each with 10 rats were studied. One group went untreated, another was treated only with radiation, and the third group received radiation plus an appropriate concentration of BPA.
}
\source{
Klein and Moeschberger (2006). \emph{Survival Analysis: Techniques for Censored and truncated data,} 2nd edition. Springer, New York.
}

\examples{
data(bnct)
 str(bnct)
 
with(bnct, mcpqdci(y = time, f = trt, event = death, TRUE))
with(bnct, mcpqrci(y = time, f = trt, event = death, TRUE))
}
\keyword{datasets}
