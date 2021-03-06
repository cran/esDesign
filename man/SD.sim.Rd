% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esDesign.R
\name{SD.sim}
\alias{SD.sim}
\title{Conduct the simulation studies of the standard design}
\usage{
SD.sim(N, rho, alpha, beta, theta, theta0, sigma0, nSim, Seed)
}
\arguments{
\item{N}{The total sample size required}

\item{rho}{The proportion of subgroup 1}

\item{alpha}{The overall Type I error rate}

\item{beta}{The \code{(1 -Power)}}

\item{theta}{The sizes of treatment effects for subgroups 1 and 2 in experimental arm}

\item{theta0}{The size of treatment effect for the control arm}

\item{sigma0}{The variance of the treatment effect}

\item{nSim}{The number of simulated studies}

\item{Seed}{The random seed}
}
\value{
A list contains,
\itemize{
  \item nTotal the total sample used
  \item The power of the specified trial. Here, the power is defined as the
          probability of rejecting the null hypothesis.
}
}
\description{
The \code{SD.sim()} is used to implement the simulation studies
   of the standard design.
}
\examples{
N <- 620
rho <- 0.5
alpha <- 0.05
beta <- 0.2
theta <- c(0.2,0.0)
theta0 <- 0
sigma0 <- 1
nSim <- 1000
Seed <- 6
SD.sim(N = N, rho = rho,
       alpha = alpha, beta = beta, theta = theta, theta0 = theta0,
       sigma0 = sigma0, nSim = nSim, Seed = Seed)

}
