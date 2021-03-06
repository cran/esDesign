% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esDesign.R
\name{sSize.norm}
\alias{sSize.norm}
\title{Sample size calculation for the standard design with continuous endpoint}
\usage{
sSize.norm(alpha, beta, theta, side, r, sigma2)
}
\arguments{
\item{alpha}{The Type I error rate or the significant level}

\item{beta}{beta The \code{(1 -Power)}}

\item{theta}{The size of treatment effect}

\item{side}{One-sided or two-sided Test}

\item{r}{The ratio of sample size between the experimental and control arms}

\item{sigma2}{The variance of the treatment effect}
}
\value{
A list contains the total sample size, and the sample sizes required
          for the experimental and control arms.
}
\description{
The \code{sSize.norm()} is used to calculate the sample size
    used in the standard design with continuous endpoint.
}
\examples{
alpha <- 0.05
beta <- 0.2
theta <- 0.2
side <- 1
r <- 1
sigma2 <- 0.8
sSize.norm(alpha = alpha, beta = beta, theta = theta,
side = side, r = r, sigma2 = sigma2)

}
