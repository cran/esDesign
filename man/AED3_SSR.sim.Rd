% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esDesign.R
\name{AED3_SSR.sim}
\alias{AED3_SSR.sim}
\title{Conduct the simulation studies of the Adaptive Enrichment Design
       (Strategy 3) with Sample Size Re-estimation Procedure based on
       Futility and Efficacy Stopping Boundaries for the continuous
       endpoint}
\usage{
AED3_SSR.sim(N1, rho, alpha, beta, theta, theta0, sigma0, pstar, nSim, Seed)
}
\arguments{
\item{N1}{The sample size used at the first stage}

\item{rho}{The proportion of subgroup 1 among the overall patients}

\item{alpha}{The overall Type I error rate}

\item{beta}{The \code{(1 - Power)}}

\item{theta}{The sizes of treatment effect in subgroups 1 and 2 with experimental treatment}

\item{theta0}{The size of treatment effect in standard treatment}

\item{sigma0}{The known variance of the treatment effect}

\item{pstar}{The \code{(1 - power)} of accepting the null hypothesis at the
interim analysis.}

\item{nSim}{The number of simulated studies.}

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
  \item Enrich01 The prevalence of adaptive enrichment of subgroup 1
  \item Enrich02 The prevalence of adaptive enrichment of subgroup 2
  \item Trigger03 The prevalence of early stopping for the situation, in which
     the treatment effect in subgroup 1 is superiority, while the treatment
     effect in subgroup 2 is inconclusive
  \item Trigger04 The prevalence of early stopping for the situation, in which
     the treatment effect in subgroup 2 is superiority, while the treatment
     effect in subgroup 2 is inconclusive
  \item ESF The probability of early stopping for futility
  \item ESE The probability of early stopping for efficacy
}
}
\description{
The \code{AED3_SSR.sim()} is used to conduct the adaptive enrichment
  design with Sample Size Re-estimation, in which futility and efficacy stopping
  boundaries are used to guide the adaptive enrichment process. For the
  adaptively enriched subgroup, we re-estimate the sample size to maintain an
  adequate conditional power meanwhile protect the overall Type I error rate.
}
\examples{
N <- 310
rho <- 0.5
alpha <- 0.05
beta <- 0.2
theta <- c(0,0)
theta0 <- 0
sigma0 <- 1
pstar <- 0.20
nSim <- 100
Seed <- 6
res <- AED3_SSR.sim(N1 = N, rho = rho, alpha = alpha,
             beta = beta, theta = theta, theta0 = theta0,
             sigma0 = sigma0, pstar = pstar, nSim = nSim,
             Seed = Seed)
}
