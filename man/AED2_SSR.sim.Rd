% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esDesign.R
\name{AED2_SSR.sim}
\alias{AED2_SSR.sim}
\title{Conduct the simulation studies of the Adaptive Enrichment Design (Strategy 2)
  with Sample Size Re-estimation Procedure}
\usage{
AED2_SSR.sim(
  N1,
  rho,
  alpha,
  beta,
  pstar,
  theta,
  theta0,
  sigma0,
  epsilon,
  nSim,
  Seed
)
}
\arguments{
\item{N1}{The sample size used in the first stage}

\item{rho}{The proportion of subgroup 1}

\item{alpha}{The overall Type I error rate}

\item{beta}{The (1 - power)}

\item{pstar}{The \code{(1 - power)} of accepting the null hypothesis at the
interim analysis.}

\item{theta}{The sizes of treatment effect in subgroups 1 and 2 with the
experimental treatment}

\item{theta0}{The size of treatment effect with the standard treatment}

\item{sigma0}{The variance of the treatment effect}

\item{epsilon}{The threshold of the difference between subgroup-specific test
statistics}

\item{nSim}{The number of simulated studies}

\item{Seed}{The random seed}
}
\value{
A list contains
\itemize{
  \item nTotal The average expected sample size
  \item H00 The probability of rejecting the null hypothesis of \eqn{H_{00}}
  \item H01 The probability of rejecting the null hypothesis of \eqn{H_{01}}
  \item H02 The probability of rejecting the null hypothesis of \eqn{H_{02}}
  \item H0  The probabilities of rejecting at least one of the null hypothesis
  \item ESF The probability of early stopping for futility
  \item ESE The probability of early stopping for efficacy
  \item Enrich01 The prevalence of adaptive enrichment of subgroup 1
  \item Enrich02 The prevalence of adaptive enrichment of subgroup 2
  \item Trigger03 The prevalence of no enrichment
}
}
\description{
The \code{AED2_SSR.sim()} is used to conduct the simulation studies
  of the Adaptive Enrichment Design (Strategy) with sample size re-estimation
  procedure. The AED2-SSR is different from the AED3-SSR, in which an
  \eqn{\epsilon}-rule is introduced to select the subgroup with larger
  subgroup-specific test statistic.
}
\examples{
N <- 310
rho <- 0.5
alpha <- 0.05
beta <- 0.2
theta <- c(0,0)
theta0 <- 0
sigma0 <- 1
epsilon <- 0.5
pstar <- 0.20
nSim <- 1000
Seed <- 6
res <- AED2_SSR.sim(N1 = N, rho = rho, alpha = alpha,
             beta = beta, theta = theta, theta0 = theta0,
             sigma0 = sigma0, pstar = pstar, epsilon = epsilon,
             nSim = nSim, Seed = Seed)
}
