\name{estimateRsqr.pRRopheticCv}
\alias{estimateRsqr.pRRopheticCv}
\title{R^2 from "pRRopheticCv" object.}
\usage{
estimateRsqr.pRRopheticCv(cvOutData, powerTranform = TRUE)
}
\arguments{
  \item{cvOutData}{an object of class "pRRopheticCv", i.e.
  the outpout of the "predictionAccuracyByCv()" function}
}
\value{
a numeric vector containing the R squared value from the
cross validation.
}
\description{
Given an object of class "pRRopheticCv", i.e. the output of
cross validation, calculate the R^2 value for the
prediction (an estimate of prediction accuracy).
}
\keyword{r-squared}

