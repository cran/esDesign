% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esDesign.R
\name{AED3_SSR.CP}
\alias{AED3_SSR.CP}
\title{Calculate the \eqn{N2} and the critical value \eqn{C} in the Adaptive
  Enrichment Design (Strategy 3) with Sample Size Re-estimation Procedure}
\usage{
AED3_SSR.CP(
  Z1 = NULL,
  delta = NULL,
  N1 = NULL,
  pstar,
  rho,
  alpha,
  beta,
  N2 = NULL
)
}
\arguments{
\item{Z1}{The test statistic obtained at the interim analysis}

\item{delta}{The standardized size of treatment effect, which can be estimated
by using \eqn{(\mu_{X} - \mu_{Y})/\sqrt{\sigma^2}}.}

\item{N1}{The sample size used at the first stage}

\item{pstar}{The \code{(1 - power)} of accepting the null hypothesis at the
interim analysis.}

\item{rho}{The proportion of subgroup 1}

\item{alpha}{The overall Type I error rate}

\item{beta}{The \code{(1 - Power)}}

\item{N2}{The pre-specified sample size used at the second stage, which is used
to conduct the conditional power analysis}
}
\value{
A list contains
\itemize{
  \item N2 The pre-specified sample size used at the second stage, which is
          used to implement the conditional power analysis
  \item Conditional.Power The value of conditional power given the value of \code{N2} in the
          conditional power analysis
  \item P.Value The corresponding P-Value used at the final analysis in the conditional
          power analysis
  \item N2.CP The re-estimated sample size of \code{N2} to ensure an adequate
          conditional power
  \item c.CP The estimated the critical value used at the final analysis based
          the conditional power
}
}
\description{
The \code{AED3_SSR.CP()} is used to calculate the sample size required
   at the second stage and the critical value used at the final analysis in the
   Adaptive Enrichment Design with Sample Size Re-estimation Procedure. In
   addition, this function can also used to conduct the conditional power
   analysis in terms of \eqn{N2}
}
\examples{
Z1 <- 1.974
delta <- 0.355
N1 <- 248
pstar <- 0.15
alpha <- 0.05
rho <- 0.5
beta <- 0.20
N2 <- 108
AED3_SSR.CP(Z1 = Z1, delta = delta, N1 = N1, pstar = pstar,
           alpha = alpha, rho = rho, beta = beta, N2 = N2)

}
