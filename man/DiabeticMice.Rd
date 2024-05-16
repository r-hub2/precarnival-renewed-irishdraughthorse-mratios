\name{DiabeticMice}
\alias{DiabeticMice}
\docType{data}
\title{
Serum albumin of diabetic mice
}
\description{
The amounts of nitrogen-bound bovine serum albumen produced by three groups of diabetic mice
}
\usage{data("DiabeticMice")}
\format{
  A data frame with 57 observations on the following 2 variables.
  \describe{
    \item{\code{group}}{a factor with levels \code{alloxan} \code{insulin} \code{normal}}
    \item{\code{response}}{Amounts of nitrogen-bound bovine serum albumen produced by the mice}
  }
}
\details{
The 57 observations of the amounts of nitrogen-bound bovine serum albumen produced by three groups of diabetic mice, these being normal, alloxan diabetic and alloxan diabetic treated with insulin.
}
\source{
Hand, D. J., Daly, F., McConway, K., Lunn, D. and Ostrowski, E. (1994). \emph{A Handbook of Small Data Sets.} Chapman & Hall/CRC, London.
}

\examples{
data(DiabeticMice)
str(DiabeticMice)
boxplot(response~group, data = DiabeticMice)

y <- DiabeticMice$response
f <- DiabeticMice$group
mcpqdci(y, f)
mcpqrci(y, f)

}
\keyword{datasets}
