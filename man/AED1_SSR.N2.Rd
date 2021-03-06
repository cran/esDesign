% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esDesign.R
\name{AED1_SSR.N2}
\alias{AED1_SSR.N2}
\title{Calculate the sample size required at the second stage of the adaptive
  enrichment design (Strategy1) with Sample Size Re-estimation Procedure}
\usage{
AED1_SSR.N2(c, z1, N1, beta)
}
\arguments{
\item{c}{The critical value used at the final analysis}

\item{z1}{The test statistic obtained at the interim analysis}

\item{N1}{The sample size used at the first stage}

\item{beta}{The (1 - power)}
}
\value{
The Value of the re-estimated sample size
}
\description{
The \code{AED1_SSR.N2()} is used to calculated the sample size
  required at the second stage of the Adaptive Enrichment Design (Strategy 1)
  with Sample Size Re-estimation Procedure.
}
\examples{
c <- 2.258
z1 <- 1.974
N1 <- 248
beta <- 0.2
AED1_SSR.N2(c = c, z1 = z1, N1 = N1, beta = beta)

}
\references{
\itemize{
  \item Zhao Yang, Ruitao Lin, Guosheng Yin and Ying Yuan. (2018) Sample Size
    Re-estimation in Adaptive Enrichment Trials. (Under Preparation)
}
}
