#' @title Commonly used \eqn{\alpha}-spending functions
#'
#' @description The \code{SigP()} is used to calculate the reduced significant
#'    level based on several widely used \eqn{\alpha}-spending functions, such as
#'    the "Pocock", "Lan-DeMets", "O`Brein-Fleming" and "Power" functions.
#'
#' @param alpha The overall Type I error rate
#' @param Info The fraction of the observed information
#' @param esFunction The specific \eqn{\alpha}-spending function. For example,
#'          \code{esFunction = "Pocock"} for the Pocock method,
#'          \code{esFunction = "LD"} for the Lan-Demets method,
#'          \code{esFunction = "OF"} for the O'Brein-Fleming method, and
#'          \code{esFunction = "Power"} for the Power method.
#' @param gamma The parameter used in the Power method. The default value is
#'          \code{gamma = 1}.
#'
#' @import stats
#'
#' @return The reduced significant level
#'
#' @examples
#' alpha <- 0.05
#' Info <- 0.5
#' esFunction = "OF"
#' SigP(alpha = alpha, Info = Info, esFunction = esFunction)
#' @export

SigP <- function(alpha, Info, esFunction = "Pocock", gamma = 1) {
  x <- NA
  if (esFunction == "LD") {
    x <- 4*(1 -  pnorm(abs(qnorm(alpha/4))/(Info)))
  } else if (esFunction == "OF") {
    x <- 2*(1 -  pnorm(abs(qnorm(alpha/2))/(Info)))
  } else if (esFunction == "Pocock") {
    x <- alpha * log(1 + (exp(1) - 1)*(Info))
  } else if (esFunction == "Power"){
    x <- alpha*Info**gamma
  } else {
    stop("Warnning: Please choose a proper method!")
  }
  return(x)
}


#' @title Sample size calculation for the standard design with continuous endpoint
#'
#' @description The \code{sSize.norm()} is used to calculate the sample size
#'     used in the standard design with continuous endpoint.
#'
#' @param alpha The Type I error rate or the significant level
#' @param beta beta The \code{(1 -Power)}
#' @param theta The size of treatment effect
#' @param side One-sided or two-sided Test
#' @param r The ratio of sample size between the experimental and control arms
#' @param sigma2 The variance of the treatment effect
#'
#' @import stats
#'
#' @return A list contains the total sample size, and the sample sizes required
#'           for the experimental and control arms.
#'
#' @examples
#' alpha <- 0.05
#' beta <- 0.2
#' theta <- 0.2
#' side <- 1
#' r <- 1
#' sigma2 <- 0.8
#' sSize.norm(alpha = alpha, beta = beta, theta = theta,
#' side = side, r = r, sigma2 = sigma2)
#'
#' @export

sSize.norm <- function(alpha, beta, theta, side, r, sigma2) {
  n2 <- (1 + 1/r)*sigma2*(qnorm(alpha/side) + qnorm(beta))^2/(theta^2)
  n1 <- r*n2
  return(list(n = ceiling(n1) + ceiling(n2),
              n2 = ceiling(n2),
              n1 = ceiling(n1)))
}
sSize.norm(alpha = 0.04, beta = 0.2, theta = 0.2, side = 1, r = 1, sigma2 = 0.8)


#' @title Conduct the simulation studies of the standard design
#'
#' @description The \code{SD.sim()} is used to implement the simulation studies
#'    of the standard design.
#'
#' @param N The total sample size required
#' @param rho The proportion of subgroup 1
#' @param alpha The overall Type I error rate
#' @param beta The \code{(1 -Power)}
#' @param theta The sizes of treatment effects for subgroups 1 and 2 in experimental arm
#' @param theta0 The size of treatment effect for the control arm
#' @param sigma0 The variance of the treatment effect
#' @param nSim The number of simulated studies
#' @param Seed The random seed
#'
#' @import stats
#'
#' @return A list contains,
#' \itemize{
#'   \item nTotal the total sample used
#'   \item The power of the specified trial. Here, the power is defined as the
#'           probability of rejecting the null hypothesis.
#' }
#'
#'
#' @examples
#' N <- 620
#' rho <- 0.5
#' alpha <- 0.05
#' beta <- 0.2
#' theta <- c(0.2,0.0)
#' theta0 <- 0
#' sigma0 <- 1
#' nSim <- 1000
#' Seed <- 6
#' SD.sim(N = N, rho = rho,
#'        alpha = alpha, beta = beta, theta = theta, theta0 = theta0,
#'        sigma0 = sigma0, nSim = nSim, Seed = Seed)
#'
#' @export
SD.sim <- function(N, rho, alpha, beta, theta, theta0, sigma0, nSim, Seed) {
  set.seed(Seed)
  r <- 1
  n1 <- r/(r + 1)*N
  c <- qnorm(1 - alpha)
  nTotal <- count1 <- 0
  theta.new <- rho*theta[1] + (1 - rho)*theta[2]
  for (i in 1:nSim) {
    x <- rnorm(n1, theta.new, sigma0)
    y <- rnorm(n1, theta0, sigma0)
    sigma <- var(c(x, y))
    z1 <- (mean(x) - mean(y))/sqrt(2 * sigma/n1)
    if (z1 > c) { count1 <- count1 + 1}
    nTotal <- nTotal + 2 * n1/nSim
  }
  return(list(nTotal = ceiling(nTotal),
              H0 = round(count1/nSim*100,1)))
}


#' @title Calculate the futility and efficacy stopping boundaries for Sample Size
#'   Re-estimation Procedure based on the conditional error function
#'
#' @description The \code{SSD.boundary()} is used to calculate the futility
#'    and efficacy stopping boundaries, meanwhile protect the overall Type I
#'    error rate at the pre-specified level.
#'
#' @param alpha The overall Type I error rate
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#'
#' @return A list contain
#' \itemize{
#'   \item upper.boundary The efficacy stopping boundary at the interim analysis
#'   \item lower.boundary The futility stopping boundary at the interim analysis
#' }
#'
#' @import stats
#'
#' @references
#' \itemize{
#'  \item Proschan MA, Hunsberger SA. Designed extension of studies based on
#'          conditional power. Biometrics 1995:1315-24. <doi:10.2307/2533262>
#' }
#'
#' @examples
#' alpha <- 0.05
#' pstar <- 0.2
#' res <- SSR.boundary(alpha = alpha, pstar = pstar)
#'
#' @export
SSR.boundary <- function(alpha, pstar) {
  zp <- qnorm(1 - pstar)
  f2 <- function(k) {
    f1 <- function(z){
      (1 - pnorm(sqrt(k^2 - z^2)))*dnorm(z)
    }
    integrate(f1, lower = zp, upper = k)$value + 1 - pnorm(k) - alpha
  }

  upper.boundary <- uniroot(f2, c(zp, 4))$root
  cat("\n")
  cat("Futility/Efficacy stopping boundaries for SSR:\n")
  cat("\n")
  cat("The futility stopping boundary = ", round(zp,3), ".\n")
  cat("The efficacy stopping boundary = ", round(upper.boundary,3), ".\n")
  return(list(upper.boundary = upper.boundary,
              lower.boundary = zp))
}

#' @title Conduct the simulation studies using SSR
#'
#' @description The \code{SSR.sim()} is used to implement the simulation studies
#'    based on the Sample Size Re-estimation Procedure.
#'
#' @param N The sample size used at the first stage. Note that this \code{N} is
#'    not the initial total sample size calculated using the standard design
#' @param rho The proportion of subgroup 1
#' @param alpha The overall Type I error rate
#' @param beta The \code{(1 - Power)}
#' @param theta The sizes of treatment effects for subgroups 1 and 2 in the
#'   experimental arm
#' @param theta0 The size of treatment effect in the control arm
#' @param sigma0 The variance of the treatment effect
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param nSim The number of simulated studies
#' @param Seed The random seed
#'
#' @return A list contains
#' \itemize{
#'   \item nTotal The average total sample size used in SSR
#'   \item H0 The power of SSR under the specific trial design. Here, the power
#'      is defined as the probability of rejecting the null hypothesis
#'   \item ESF The percentage of early stopping for futility
#'   \item ESE The percentage of early stopping for efficacy
#' }
#'
#' @import stats
#'
#' @references
#' \itemize{
#'  \item Proschan MA, Hunsberger SA. Designed extension of studies based on
#'     conditional power. Biometrics 1995:1315-1324. <doi:10.2307/2533262>
#' }
#'
#' @examples
#' N <- 310
#' rho <- 0.5
#' alpha <- 0.05
#' beta <- 0.2
#' pstar <- 0.2
#' theta <- c(0.2,0)
#' theta0 <- 0
#' sigma0 <- 1.0
#' nSim <- 1000
#' Seed <- 6
#' res <- SSR.sim(N = N, rho = rho, alpha = alpha, beta = beta, theta = theta,
#'         theta0 = theta0, sigma0 = sigma0, pstar = pstar,
#'         nSim = nSim, Seed = Seed)
#' @export

SSR.sim <- function(N, rho, alpha, beta, theta, theta0, sigma0, pstar, nSim, Seed) {
  set.seed(6)
  ## Sample size
  r <- 1
  n1 <- r/(r + 1) * N
  #n11 <- N * rho
  #n12 <- N * (1 - rho)
  ## Boundaries
  zp <- qnorm(1 - pstar) ## lower boundary of the intervals
  cat("\n")
  cat("## Summary of SSR: \n")
  k <- SSR.boundary(alpha = alpha, pstar = pstar)$upper.boundary
  zb <- qnorm(1 - beta)
  ##k <- 2.223702
  ## Results
  nTotal <- count1 <- count2 <- count3 <- 0
  #res <- NULL
  theta.new <- rho*theta[1] + (1 - rho)*theta[2]
  for (i in 1:nSim) {
    x1 <- rnorm(n1, theta.new, sigma0)  ## outcome for subgroup 1 in stage 1
    y1 <- rnorm(n1, theta0, sigma0)  ## outcome for control/2  in stage 1
    sigma <- var(c(x1, y1)) ## variance of the outcome

    ## Interim analysis for stage 1
    z1  <- (mean(x1) - mean(y1))/sqrt(2 * sigma/n1)    ## Test statistics for total population

    ## Sample Size Calculation
    nTotal <- nTotal + 2 * n1/nSim
    if (z1 < zp) {
      count1 <- count1 + 1    ## early stop for fulitity in subgroup 1
    } else if (z1 >= k) {
      count2 <- count2 + 1              ## (theta2 = 0) early stop for efficacy to subgroup 2
    } else {                              ## Sample Size Re-estimation for subgroup 2, in stage 2
      #mu <- mean(x1) - mean(y1)        ## Data obtained from stage 1
      A <- 1 - pnorm(sqrt(k^2 - z1^2))  ##
      za <- qnorm(1 - A)
      n2 <- n1 * (za + zb)^2/z1^2       ## Changed: original one n2 = n1/2 * (za + zb)^2/z1^2
      nTotal <- nTotal + 2 * n2/nSim
      c <- (z1^2 + za * (za + zb))/sqrt(z1^2 + (za + zb)^2)
      x2 <- rnorm(n2, theta.new, sigma0)
      y2 <- rnorm(n2, theta0, sigma0)
      x <- c(x1, x2)
      y <- c(y1, y2)
      sigma <- var(c(x, y))
      z <- (mean(x) - mean(y))/sqrt(2 * sigma/(n1 + n2))
      if (z >= c) { count3 <- count3 + 1 }
    }

  } # nSim
  cat("\n")
  cat("The expected sample size = ", ceiling(nTotal), "\n")
  cat("The overall power = ", round((count2 + count3)/nSim*100, 1), "%. \n")
  cat("\n")
  cat("Pr(Early stopping for futility) = ", round(count2/nSim*100,1), "%. \n")
  cat("Pr(Early stopping for efficacy) = ", round(count1/nSim*100,1), "%. \n")
  cat("\n")
  return(list(nTotal = ceiling(nTotal),
              H0 = round((count2 + count3)/nSim*100, 1),
              ESF = round(count1/nSim*100,1),
              ESE = round(count2/nSim*100,1)))

}

#' @title Calculate the \eqn{N2} and the critical value \eqn{C} in Sample Size
#'   Re-estimation Procedure
#'
#' @description The \code{SSR.CP()} is used to calculate the sample size required
#'    at the second stage and the critical value used at the final analysis. In
#'    addition, this function can also used to conduct the conditional power
#'    analysis in terms of \eqn{N2}
#'
#' @param Z1 The test statistic obtained at the interim analysis
#' @param delta The standardized size of treatment effect, which can be estimated
#'   by using \eqn{(\mu_{X} - \mu_{Y})/\sqrt{\sigma^2}}.
#' @param N1 The sample size used at the first stage
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param alpha The overall Type I error rate
#' @param beta The \code{(1 - Power)}
#' @param N2 The pre-specified sample size used at the second stage, which is used
#'   to conduct the conditional power analysis
#'
#' @return A list contains
#' \itemize{
#'   \item N2 The pre-specified sample size used at the second stage, which is
#'           used to implement the conditional power analysis
#'   \item Conditional.Power The value of conditional power given the value of \code{N2} in the
#'           conditional power analysis
#'   \item P.Value The corresponding P-Value used at the final analysis in the conditional
#'           power analysis
#'   \item N2.CP The re-estimated sample size of \code{N2} to ensure an adequate
#'           conditional power
#'   \item c.CP The estimated the critical value used at the final analysis based
#'           the conditional power
#' }
#'
#' @import stats
#'
#' @references
#' \itemize{
#'  \item Proschan MA, Hunsberger SA. Designed extension of studies based on
#'     conditional power. Biometrics 1995:1315-1324. <doi:10.2307/2533262>
#' }
#'
#' @examples
#' Z1 <- 1.527
#' delta <- 0.137
#' N1 <- 248
#' pstar <- 0.15
#' alpha <- 0.05
#' beta <- 0.2
#' res <- SSR.CP(Z1 = Z1, delta = delta, N1 = N1,
#'        pstar = pstar, alpha = alpha, beta = beta)
#'
#' @export
SSR.CP <- function(Z1 = NULL, delta = NULL, N1 = NULL,
                   pstar, alpha, beta, N2 = NULL){
  zp <- qnorm(1 - pstar)
  k <- SSR.boundary(alpha = alpha, pstar = pstar)$upper.boundary
  z1 <- Z1
  n1 <- N1*0.5
  A <- 0
  if (z1 >= zp & z1 < k) {
    A <- 1 - pnorm(sqrt(k^2 - z1^2))
  } else if (z1 >= k) {
    A <- 1
  }
  za <- qnorm(1 - A)
  cp <- p <- NULL
  if (!is.null(N2)) {
    cp <- 1 - pnorm(za - sqrt(N2/2)*delta)
    c <- (sqrt(N1)*z1 + sqrt(N2)*za)/(sqrt(N1 + N2))
    p <- 1 - pnorm(c)
    cat("\n")
    cat("Conditional power analysis: \n")
    cat("\n")
    cat("The conditional power = ", round(cp*100,1), "%, given N2 = ", N2, ".\n")
    cat("The exact critical value used at the final analysis = ", round(c, 3),
        ", with the corresponding P-Value = ",round(p,4), ".\n")
    cat("\n")
  }
  zb <- qnorm(1 - beta)
  N2.CP <- n1 * (za + zb)^2/z1^2
  c.CP <- (z1^2 + za*(za + zb))/(sqrt(z1^2 + (za + zb)^2))
  cat("\n")
  cat("The N2 required at the second stage = ", ceiling(N2.CP)*2, ".\n")
  cat("The estimated critical value used at final analysis = ", round(c.CP, 3), ".\n")
  cat("\n")

  return(list(N2 = N2, Conditional.Power = cp,
              P.Value = p, N2.CP = 2*ceiling(N2.CP), c.CP = c.CP))
}


#' @title Calculate the futility and efficacy stopping boundaries in Adaptive
#'   enrichment design (Strategy 3) with Sample Size Re-estimation Procedure
#'   for the continuous endpoint
#'
#' @description The \code{AED3_SSR.boundary()} is used to calculate the futility
#'   and efficacy stopping boundaries in the Adaptive Enrichment Design with
#'   Sample Size Re-estimation Procedure.
#'
#' @param rho The proportion of subgroup 1
#' @param alpha The overall Type I error rate
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#'
#' @return A list contains
#' \itemize{
#'   \item upper.boundary The upper or the efficacy stopping boundary
#'   \item lower.boundary The lower or the futility stopping boundary
#' }
#'
#' @examples
#' rho <- 0.5
#' alpha <- 0.05
#' pstar <- 0.15
#' res <- AED3_SSR.boundary(rho = rho, alpha = alpha, pstar = pstar)
#' @export
AED3_SSR.boundary <- function(rho, alpha, pstar) {
  N <- 500
  n11 <- N * rho  ## total sample size used in stage 1
  n12 <- N * (1 - rho) ## total sample size used to subgroup 2 in stage 2
  zp <- qnorm(1 - pstar)    ## quantile of power (conditional power)
  a <- sqrt(n11/(n11 + n12))  ## weight for subgroup 1
  b <- sqrt(n12/(n11 + n12))  ## weight for subgroup 2

  f <- function(k) {
    InnerFunc = function(y1, y2, a, b) {
      f <- 0
      if (y1 < k & y1 >= zp) {
        f <- 1/a * 1/b * dnorm(y2/a) * dnorm((y1 - y2)/b) * (1 - pnorm(sqrt(k^2 - y1^2)))
      }
      if (y1 >= k) {
        f <- 1/a * 1/b * dnorm(y2/a) * dnorm((y1 - y2)/b)
      }
      return(f)
    }

    InnerIntegral = function(y1) {
      sapply(y1, function(y1) {
        integrate(InnerFunc, zp * a, y1 - b * zp,
                  y1 = y1, a = a, b = b)$value
      })
    }

    InnerFunc1 <- function(z) {
      return((1 - pnorm(sqrt(k^2 - z^2))) * dnorm(z))
    }

    integrate(InnerIntegral, a * zp + b * zp, Inf)$value +
      (integrate(InnerFunc1, lower = zp, upper = k)$value +
         1 - pnorm(k)) * (1 - pstar) * 2 - alpha
  }

  k = uniroot(f, c(zp, 4))$root

  cat("\n")
  cat("Futility/Efficacy stopping boundaries for AED3-SSR:\n")
  cat("\n")
  cat("The futility stopping boundary = ", round(zp,3), ".\n")
  cat("The efficacy stopping boundary = ", round(k,3), ".\n")
  cat("\n")
  return(list(upper.boundary = k, lower.boundary = zp))
}


#' @title Conduct the simulation studies of the Adaptive Enrichment Design
#'        (Strategy 3) with Sample Size Re-estimation Procedure based on
#'        Futility and Efficacy Stopping Boundaries for the continuous
#'        endpoint
#'
#' @description The \code{AED3_SSR.sim()} is used to conduct the adaptive enrichment
#'   design with Sample Size Re-estimation, in which futility and efficacy stopping
#'   boundaries are used to guide the adaptive enrichment process. For the
#'   adaptively enriched subgroup, we re-estimate the sample size to maintain an
#'   adequate conditional power meanwhile protect the overall Type I error rate.
#'
#' @param N1 The sample size used at the first stage
#' @param rho The proportion of subgroup 1 among the overall patients
#' @param alpha The overall Type I error rate
#' @param beta The \code{(1 - Power)}
#' @param theta The sizes of treatment effect in subgroups 1 and 2 with experimental treatment
#' @param theta0 The size of treatment effect in standard treatment
#' @param sigma0 The known variance of the treatment effect
#' @param pstar  The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param nSim The number of simulated studies.
#' @param Seed The random seed
#'
#' @return A list contains
#' \itemize{
#'   \item nTotal The average expected sample size
#'   \item H00 The probability of rejecting the null hypothesis of \eqn{H_{00}}
#'   \item H01 The probability of rejecting the null hypothesis of \eqn{H_{01}}
#'   \item H02 The probability of rejecting the null hypothesis of \eqn{H_{02}}
#'   \item H0  The probabilities of rejecting at least one of the null hypothesis
#'   \item Enrich01 The prevalence of adaptive enrichment of subgroup 1
#'   \item Enrich02 The prevalence of adaptive enrichment of subgroup 2
#'   \item Trigger03 The prevalence of early stopping for the situation, in which
#'      the treatment effect in subgroup 1 is superiority, while the treatment
#'      effect in subgroup 2 is inconclusive
#'   \item Trigger04 The prevalence of early stopping for the situation, in which
#'      the treatment effect in subgroup 2 is superiority, while the treatment
#'      effect in subgroup 2 is inconclusive
#'   \item ESF The probability of early stopping for futility
#'   \item ESE The probability of early stopping for efficacy
#' }
#'
#' @import stats
#'
#' @examples
#' N <- 310
#' rho <- 0.5
#' alpha <- 0.05
#' beta <- 0.2
#' theta <- c(0,0)
#' theta0 <- 0
#' sigma0 <- 1
#' pstar <- 0.20
#' nSim <- 100
#' Seed <- 6
#' res <- AED3_SSR.sim(N1 = N, rho = rho, alpha = alpha,
#'              beta = beta, theta = theta, theta0 = theta0,
#'              sigma0 = sigma0, pstar = pstar, nSim = nSim,
#'              Seed = Seed)
#' @export
#'
AED3_SSR.sim <- function(N1, rho, alpha, beta, theta, theta0, sigma0, pstar, nSim, Seed) {
  set.seed(Seed)
  r <- 1
  n1 <- r/(r + 1) * N1
  n11 <- N1 * rho
  n12 <- N1 * (1 - rho)
  zp <- qnorm(1 - pstar)
  cat("\n")
  cat("## Summary of AED3-SSR:    \n")
  cat("\n")
  k <- AED3_SSR.boundary(rho = rho, alpha = alpha, pstar = pstar)$upper.boundary
  zb <- qnorm(1 - beta)
  nTotal <- count1 <- count2 <- count3 <- count4 <- count5 <- count6 <- 0
  count7 <- count8 <- count9 <- 0
  trigger02 <- trigger01 <- trigger03 <- trigger04 <- 0
  res <- NULL
  for (i in 1:nSim) {
    x1 <- rnorm(n11/2, theta[1], sigma0)
    y1 <- rnorm(n11/2, theta0, sigma0)
    x2 <- rnorm(n12/2, theta[2], sigma0)
    y2 <- rnorm(n12/2, theta0, sigma0)
    x <- c(x1, x2)
    y <- c(y1, y2)
    sigma <- var(c(x1, x2, y1, y2))
    z1  <- (mean(x) - mean(y))/sqrt(2 * sigma/n1)
    z11 <- (mean(x1) - mean(y1))/sqrt(4 * sigma/n11)
    z12 <- (mean(x2) - mean(y2))/sqrt(4 * sigma/n12)

    nTotal <- nTotal + 2 * n1/nSim
    if (z11 < zp) {
      if (z12 >= k) {
        count3 <- count3 + 1
        count6 <- count6 + 1
      } else if (z12 < k & z12 >= zp) {
        mu <- mean(x2) - mean(y2)
        A <- 1 - pnorm(sqrt(k^2 - z12^2))
        za <- qnorm(1 - A)
        n2 <- n12/2 * (za + zb)^2/z12^2
        nTotal <- nTotal + 2 * n2/nSim
        c <- (z12^2 + za * (za + zb))/sqrt(z12^2 + (za + zb)^2)
        x2 <- c(x2, rnorm(n2, theta[2], sigma0))
        y2 <- c(y2, rnorm(n2, theta0, sigma0))
        x <- c(x1, x2)
        y <- c(y1, y2)
        sigma <- var(c(x, y))
        z <- (mean(x2) - mean(y2))/sqrt(2 * sigma/(n12/2 + n2))
        if (z >= c) { count3 <- count3 + 1 }
        trigger02 <- trigger02 + 1
      } else {
        count4 <- count4 + 1
      }
    } else if (z11 > zp & z12 < zp) {
      if (z11 >= k) {
        count2 <- count2 + 1
        count7 <- count7 + 1
      } else if (z11 < k & z11 >= zp) {
        mu <- mean(x1) - mean(y1)
        A <- 1 - pnorm(sqrt(k^2 - z11^2))
        za <- qnorm(1 - A)
        n2 <- n11/2 * (za + zb)^2/z11^2
        nTotal <- nTotal + 2 * n2/nSim
        c <- (z11^2 + za * (za + zb))/sqrt(z11^2 + (za + zb)^2)
        x1 <- c(x1, rnorm(n2, theta[1], sigma0))
        y1 <- c(y1, rnorm(n2, theta0, sigma0))
        x <- c(x1, x2)
        y <- c(y1, y2)
        sigma <- var(c(x, y))
        z <- (mean(x1) - mean(y1))/sqrt(2 * sigma/(n11/2 + n2))
        if (z >= c) {
          count2 <- count2 + 1
        }
        trigger01 <- trigger01 + 1
      }
    } else if (z11 > zp & z12 > zp) {
      if (z1 >= k) {
        count1 <- count1 + 1
        count5 <- count5 + 1
      } else {
        if (z12 < k & z11 < k) {
          res <- c(res, 1 - pnorm(sqrt(k^2 - z1^2)))
          mu <- mean(x) - mean(y)
          A <- 1 - pnorm(sqrt(k^2 - z1^2))
          za <- qnorm(1 - A)
          n2 <- n1 * (za + zb)^2/z1^2
          n21 <- 2 * round(n2 * rho) + 2
          n22 <- 2 * round(n2 * (1 - rho)) + 2
          nTotal <- nTotal + (n21 + n22)/nSim
          c <- (z1^2 + za * (za + zb))/sqrt(z1^2 + (za + zb)^2)
          x1 <- c(x1, rnorm(n21/2, theta[1], sigma0))
          y1 <- c(y1, rnorm(n21/2, theta0, sigma0))
          x2 <- c(x2, rnorm(n22/2, theta[2], sigma0))
          y2 <- c(y2, rnorm(n22/2, theta0, sigma0))
          x <- c(x1, x2)
          y <- c(y1, y2)
          sigma <- var(c(x, y))
          z <- (mean(x) - mean(y))/sqrt(2 * sigma/(n1 + n2))
          if (z >= c) {
            count1 <- count1 + 1
          }
        } else if (z11 >= k & z12 >= k) {
          count1 <- count1 + 1
          count5 <- count5 + 1
        } else if (z11 < k & z12 >= k) {
          count3 <- count3 + 1
          count9 <- count9 + 1
          trigger04 <- trigger04 + 1
        } else if (z11 >= k & z12 < k) {
          count2 <- count2 + 1
          count8 <- count8 + 1
          trigger03 <- trigger03 + 1
        }
      }
    }
  } # nSim
  res <- matrix(NA, ncol = 5, nrow = 1)
  res[,1] <- ceiling(nTotal)
  res[,2] <- round((count1 + count2 + count3)/nSim*100,1)
  res[,3] <- round(count1/nSim*100,1)
  res[,4] <- round(count2/nSim*100,1)
  res[,5] <- round(count3/nSim*100,1)
  colnames(res) <- c("ESS", "H0", "H00", "H01", "H02")
  row.names(res) <- "%"
  cat("\n")
  cat("The expected sample size and overall power: \n")
  print(res)
  cat("\n")
  cat("Pr(Early stopping for futility) = ", round(count4/nSim*100,1), "%.\n")
  cat("Pr(Early stopping for efficacy) = ",
      round((count5 + count6 + count7 + count8 + count9)/nSim*100,1), "%.\n")
  cat("Pr(Enrich subgroup 1) = ", round(trigger01/nSim*100,1), "%.\n")
  cat("Pr(Enrich subgroup 2) = ", round(trigger02/nSim*100,1), "%.\n")

  return(list(nTotal = ceiling(nTotal), H00 = round(count1/nSim*100,1),
              H01 = round(count2/nSim*100,1), H02 = round(count3/nSim*100,1),
              H0 = round((count1 + count2 + count3)/nSim*100,1),
              Enrich01 = round(trigger01*100/nSim,1),
              Enrich02 = round(trigger02*100/nSim,1),
              Trigger03 = round(trigger03*100/nSim,1),
              Trigger04 = round(trigger04*100/nSim,1),
              ESF = round(count4/nSim*100,1),
              ESE = round((count5 + count6 + count7 + count8 + count9)/nSim*100,1)
  ))
}



#' @title Calculate the \eqn{N2} and the critical value \eqn{C} in the Adaptive
#'   Enrichment Design (Strategy 3) with Sample Size Re-estimation Procedure
#'
#' @description The \code{AED3_SSR.CP()} is used to calculate the sample size required
#'    at the second stage and the critical value used at the final analysis in the
#'    Adaptive Enrichment Design with Sample Size Re-estimation Procedure. In
#'    addition, this function can also used to conduct the conditional power
#'    analysis in terms of \eqn{N2}
#'
#' @param Z1 The test statistic obtained at the interim analysis
#' @param delta The standardized size of treatment effect, which can be estimated
#'   by using \eqn{(\mu_{X} - \mu_{Y})/\sqrt{\sigma^2}}.
#' @param N1 The sample size used at the first stage
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param rho The proportion of subgroup 1
#' @param alpha The overall Type I error rate
#' @param beta The \code{(1 - Power)}
#' @param N2 The pre-specified sample size used at the second stage, which is used
#'   to conduct the conditional power analysis
#'
#' @return A list contains
#' \itemize{
#'   \item N2 The pre-specified sample size used at the second stage, which is
#'           used to implement the conditional power analysis
#'   \item Conditional.Power The value of conditional power given the value of \code{N2} in the
#'           conditional power analysis
#'   \item P.Value The corresponding P-Value used at the final analysis in the conditional
#'           power analysis
#'   \item N2.CP The re-estimated sample size of \code{N2} to ensure an adequate
#'           conditional power
#'   \item c.CP The estimated the critical value used at the final analysis based
#'           the conditional power
#' }
#'
#' @import stats
#'
#' @examples
#' Z1 <- 1.974
#' delta <- 0.355
#' N1 <- 248
#' pstar <- 0.15
#' alpha <- 0.05
#' rho <- 0.5
#' beta <- 0.20
#' N2 <- 108
#' AED3_SSR.CP(Z1 = Z1, delta = delta, N1 = N1, pstar = pstar,
#'            alpha = alpha, rho = rho, beta = beta, N2 = N2)
#'
#' @export
AED3_SSR.CP <- function(Z1 = NULL, delta = NULL, N1 = NULL,
                        pstar, rho, alpha, beta, N2 = NULL){
  zp <- qnorm(1 - pstar)
  k <- AED3_SSR.boundary(rho = rho, alpha = alpha, pstar = pstar)$upper.boundary
  z1 <- Z1
  n1 <- N1*0.5
  A <- 0
  if (z1 >= zp & z1 < k) {
    A <- 1 - pnorm(sqrt(k^2 - z1^2))
  } else if (z1 >= k) {
    A <- 1
  }
  za <- qnorm(1 - A)
  if (!is.null(N2)) {
    cp <- 1 - pnorm(za - sqrt(N2/4)*delta)
    c <- (sqrt(N1)*z1 + sqrt(N2)*za)/(sqrt(N1 + N2))
    p <- 1 - pnorm(c)
    cat("\n")
    cat("Conditional power analysis: \n")
    cat("\n")
    cat("The conditional power = ", round(cp*100,1), "%, given N2 = ", N2, ".\n")
    cat("The exact critical value used at the final analysis = ", round(c, 3),
        ", with the corresponding P-Value = ", round(p, 4), ".\n")
    cat("\n")
  }

  zb <- qnorm(1 - beta)
  N2.CP <- n1/2 * (za + zb)^2/z1^2
  c.CP <- (z1^2 + za*(za + zb))/(sqrt(z1^2 + (za + zb)^2))
  cat("AED3-SSR: \n")
  cat("The N2 required at the second stage = ", ceiling(N2.CP)*2, ".\n")
  cat("The estimated critical value used at the final analysis = ", round(c.CP,3), ".\n")
  cat("\n")
  return(list(N2 = N2, Conditional.Power = cp, P.Value = p,
              N2.CP = 2*ceiling(N2.CP),
              c.CP = c.CP))
}


#' @title Calculate the futility and efficacy stopping boundaries of the Adaptive
#'   Enrichment Design (Strategy 2) with Sample Size Re-estimation Procedure
#'
#' @description The \code{AED2_SSR.boundary()} is used to calculate the futility
#'   and efficacy stopping boundaries of the Adaptive Enrichment Design (strategy 2)
#'   with Sample Size Re-estimation Procedure. In the AED2-SSR design, an
#'   \eqn{\epsilon}-rule is introduced to select the subgroup with larger test
#'   statistic. In practice, the value of \eqn{\epsilon} should be calibrated to
#'   fit the requirement of the trial.
#'
#' @param rho The proportion of subgroup 1
#' @param alpha The overall Type I error rate
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param epsilon The threshold of difference between the subgroup-specific test
#'    statistics
#'
#' @import stats
#'
#' @return A list contains
#' \itemize{
#'   \item upper.boundary The upper and efficacy stopping boundary
#'   \item lower.boundary The lower and futility stopping boundary
#' }
#'
#' @examples
#' rho <- 0.5
#' alpha <- 0.05
#' pstar <- 0.15
#' epsilon <- 0.5
#' AED2_SSR.boundary(rho = rho, alpha = alpha, pstar = pstar, epsilon = epsilon)
#' @export
AED2_SSR.boundary <- function(rho, alpha, pstar, epsilon) {
  N <- 500
  options(warn=-1)
  #alpha <- 0.05
  #p <- 0.2
  #epi <- 0.5
  zp <- qnorm(1 - pstar)
  n11 <- rho * N
  n12 <- (1 - rho) * N
  a <- sqrt(n11)/sqrt(n11 + n12)
  b <- sqrt(n12)/sqrt(n11 + n12)

  f <- function(k) {
    InnerFunc1 = function(y1, y2) {
      f <- dnorm(y1) * dnorm(y2)
      return(f)
    }
    InnerFunc2 = function(y1, y2, k) {
      f <- (1 - pnorm(sqrt(k^2 - y1^2))) * dnorm(y1) * dnorm(y2)
      return(f)
    }
    InnerFunc3 = function(y1, y2, a, b, k) {
      y = a * y1 + b * y2
      return(ifelse(y >= k,
                    dnorm(y1) * dnorm(y2),
                    (1 - pnorm(sqrt(k^2 - y^2))) * dnorm(y1) * dnorm(y2)))
    }


    InnerIntegral1 = function(y1) {
      sapply(y1, function(y1) {
        integrate(InnerFunc1,
                  -Inf, y1 - epsilon,
                  y1 = y1)$value
      })
    }
    InnerIntegral2 = function(y1) {
      sapply(y1, function(y1) {
        integrate(InnerFunc2,
                  -Inf, y1 - epsilon,
                  y1 = y1, k = k)$value
      })
    }
    InnerIntegral3 = function(y1) {
      sapply(y1, function(y1) {
        integrate(InnerFunc3,
                  y1 - epsilon, y1, y1 = y1,
                  a = a, b = b, k = k)$value
      })
    }
    InnerIntegral4 = function(y2) {
      sapply(y2, function(y2) {
        integrate(InnerFunc3,
                  y2 - epsilon, y2, y2 = y2,
                  a = a, b = b, k = k)$value
      })
    }

    integrate(InnerIntegral1, k, Inf)$value * 2 +
      integrate(InnerIntegral2, zp, k)$value * 2 +
      integrate(InnerIntegral3,zp, Inf)$value +
      integrate(InnerIntegral4,zp, Inf)$value - alpha

  }

  k = uniroot(f, c(zp, 4))$root

  cat("\n")
  cat("Futility/Efficacy stopping boundaries for AED2-SSR:\n")
  cat("\n")
  cat("The futility stopping boundary =", round(zp,3), ".\n")
  cat("The efficacy stopping boundary =", round(k,3), ".\n")

  return(list(upper.boundary = k,
              lower.boundary = zp))
}

#' @title Conduct the simulation studies of the Adaptive Enrichment Design (Strategy 2)
#'   with Sample Size Re-estimation Procedure
#'
#' @description The \code{AED2_SSR.sim()} is used to conduct the simulation studies
#'   of the Adaptive Enrichment Design (Strategy) with sample size re-estimation
#'   procedure. The AED2-SSR is different from the AED3-SSR, in which an
#'   \eqn{\epsilon}-rule is introduced to select the subgroup with larger
#'   subgroup-specific test statistic.
#'
#' @param N1 The sample size used in the first stage
#' @param rho The proportion of subgroup 1
#' @param alpha The overall Type I error rate
#' @param beta The (1 - power)
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param theta The sizes of treatment effect in subgroups 1 and 2 with the
#'    experimental treatment
#' @param theta0 The size of treatment effect with the standard treatment
#' @param sigma0 The variance of the treatment effect
#' @param epsilon The threshold of the difference between subgroup-specific test
#'   statistics
#' @param nSim The number of simulated studies
#' @param Seed The random seed
#'
#' @return A list contains
#' \itemize{
#'   \item nTotal The average expected sample size
#'   \item H00 The probability of rejecting the null hypothesis of \eqn{H_{00}}
#'   \item H01 The probability of rejecting the null hypothesis of \eqn{H_{01}}
#'   \item H02 The probability of rejecting the null hypothesis of \eqn{H_{02}}
#'   \item H0  The probabilities of rejecting at least one of the null hypothesis
#'   \item ESF The probability of early stopping for futility
#'   \item ESE The probability of early stopping for efficacy
#'   \item Enrich01 The prevalence of adaptive enrichment of subgroup 1
#'   \item Enrich02 The prevalence of adaptive enrichment of subgroup 2
#'   \item Trigger03 The prevalence of no enrichment
#' }
#'
#' @import stats
#'
#' @examples
#' N <- 310
#' rho <- 0.5
#' alpha <- 0.05
#' beta <- 0.2
#' theta <- c(0,0)
#' theta0 <- 0
#' sigma0 <- 1
#' epsilon <- 0.5
#' pstar <- 0.20
#' nSim <- 1000
#' Seed <- 6
#' res <- AED2_SSR.sim(N1 = N, rho = rho, alpha = alpha,
#'              beta = beta, theta = theta, theta0 = theta0,
#'              sigma0 = sigma0, pstar = pstar, epsilon = epsilon,
#'              nSim = nSim, Seed = Seed)
#' @export
#'
AED2_SSR.sim <- function(N1, rho, alpha, beta, pstar, theta, theta0,
                         sigma0, epsilon, nSim, Seed) {
  set.seed(Seed)
  r <- 1
  n1 <- r/(r + 1) * N1
  n11 <- N1 * rho
  n12 <- N1 * (1 - rho)
  zp <- qnorm(1 - pstar)
  cat("\n")
  cat("## Summary of AED2-SSR:    \n")
  cat("\n")
  k <- AED2_SSR.boundary(rho = rho, alpha = alpha, pstar = pstar,
                         epsilon = epsilon)$upper.boundary
  zb <- qnorm(1 - beta)
  res  <- NULL
  count1 <- count2 <- count3 <- count4 <- count5 <- count6 <- 0
  trigger01 <- trigger02 <- trigger03 <-  ESF <- ESE <- 0
  nTotal <- 0
  for (i in 1:nSim) {
    x1 <- rnorm(n11/2, theta[1], sigma0)
    y1 <- rnorm(n11/2, theta0, sigma0)
    x2 <- rnorm(n12/2, theta[2], sigma0)
    y2 <- rnorm(n12/2, theta0, sigma0)
    x <- c(x1, x2)
    y <- c(y1, y2)
    sigma <- var(c(x1, x2, y1, y2))
    z1 <- (mean(x) - mean(y))/sqrt(2 * sigma/n1)
    z11 <- (mean(x1) - mean(y1))/sqrt(4 * sigma/n11)
    z12 <- (mean(x2) - mean(y2))/sqrt(4 * sigma/n12)
    nTotal <- nTotal + 2 * n1/nSim

    opt <- which.max(c(z11, z12))
    if (opt == 1) {
      count4 <- count4 + 1
      if (z11 < zp) {
        ESF <- ESF + 1
        next
      } else if (z11 > k & z12 < (z11 - epsilon)) {
        ESE <- ESE + 1
        count2 <- count2 + 1
        next
      } else if (z11 > zp & z12 < (z11 - epsilon)) {
        mu <- mean(x1) - mean(y1)
        A <- 1 - pnorm(sqrt(k^2 - z11^2))
        za <- qnorm(1 - A)
        n2 <- n11/2 * (za + zb)^2/z11^2
        nTotal <- nTotal + 2 * n2/nSim
        c <- (z11^2 + za * (za + zb))/sqrt(z11^2 + (za + zb)^2)
        x1 <- c(x1, rnorm(n2, theta[1], sigma0))
        y1 <- c(y1, rnorm(n2, theta0, sigma0))
        x <- c(x1, x2)
        y <- c(y1, y2)
        sigma <- var(c(x, y))
        z <- (mean(x1) - mean(y1))/sqrt(2 * sigma/(n11/2 + n2))
        if (z >= c) {
          count2 <- count2 + 1
        }
        trigger01 <- trigger01 + 1
      } else {
        if (z1 >= k) {
          ESE <- ESE + 1
          count1 <- count1 + 1
          next
        }
        res <- c(res, 1 - pnorm(sqrt(k^2 - z1^2)))
        mu <- mean(x) - mean(y)
        A <- 1 - pnorm(sqrt(k^2 - z1^2))
        za <- qnorm(1 - A)
        n2 <- n1 * (za + zb)^2/z1^2
        n21 <- 2 * round(n2 * rho)
        n22 <- 2 * round(n2 * (1 - rho))
        nTotal <- nTotal + (n21 + n22)/nSim
        c <- (z1^2 + za * (za + zb))/sqrt(z1^2 + (za + zb)^2)
        x1 <- c(x1, rnorm(n21/2, theta[1], sigma0))
        y1 <- c(y1, rnorm(n21/2, theta0, sigma0))
        x2 <- c(x2, rnorm(n22/2, theta[2], sigma0))
        y2 <- c(y2, rnorm(n22/2, theta0, sigma0))
        x <- c(x1, x2)
        y <- c(y1, y2)
        sigma <- var(c(x, y))
        z <- (mean(x) - mean(y))/sqrt(2 * sigma/(n1 + n2))
        if (z >= c) {
          count1 <- count1 + 1
        }
        trigger03 <- trigger03 + 1
      }

    }

    if (opt == 2) {
      count5 <- count5 + 1
      if (z12 < zp) {
        ESF <- ESF + 1
        next
      } else if (z12 > k & z11 < (z12 - epsilon)) {
        ESE <- ESE + 1
        count3 <- count3 + 1
        next
      } else if (z12 > zp & z11 < (z12 - epsilon)) {
        mu <- mean(x2) - mean(y2)
        A <- 1 - pnorm(sqrt(k^2 - z12^2))
        za <- qnorm(1 - A)
        n2 <- n12/2 * (za + zb)^2/z12^2
        nTotal <- nTotal + 2 * n2/nSim
        c <- (z12^2 + za * (za + zb))/sqrt(z12^2 + (za + zb)^2)
        x2 <- c(x2, rnorm(n2, theta[2], sigma0))
        y2 <- c(y2, rnorm(n2, theta0, sigma0))
        x <- c(x1, x2)
        y <- c(y1, y2)
        sigma <- var(c(x, y))
        z <- (mean(x2) - mean(y2))/sqrt(2 * sigma/(n12/2 + n2))
        if (z >= c) {
          count3 <- count3 + 1
        }
        trigger02 <- trigger02 + 1
      } else {
        if (z1 >= k) {
          ESE <- ESE + 1
          count6 <- count6 + 1
          next
        }
        res <- c(res, 1 - pnorm(sqrt(k^2 - z1^2)))
        mu <- mean(x) - mean(y)
        A <- 1 - pnorm(sqrt(k^2 - z1^2))
        za <- qnorm(1 - A)
        n2 <- n1 * (za + zb)^2/z1^2
        n21 <- 2 * round(n2 * rho)
        n22 <- 2 * round(n2 * (1 - rho))
        nTotal <- nTotal + (n21 + n22)/nSim
        c <- (z1^2 + za * (za + zb))/sqrt(z1^2 + (za + zb)^2)
        x1 <- c(x1, rnorm(n21/2, theta[1], sigma0))
        y1 <- c(y1, rnorm(n21/2, theta0, sigma0))
        x2 <- c(x2, rnorm(n22/2, theta[2], sigma0))
        y2 <- c(y2, rnorm(n22/2, theta0, sigma0))
        x <- c(x1, x2)
        y <- c(y1, y2)
        sigma <- var(c(x, y))
        z <- (mean(x) - mean(y))/sqrt(2 * sigma/(n1 + n2))
        if (z >= c) {
          count6 <- count6 + 1
        }
        trigger03 <- trigger03 + 1
      }

    }

  }
  res <- matrix(NA, ncol = 5, nrow = 1)
  res[,1] <- ceiling(nTotal)
  res[,2] <- round((count1 + count2 + count3 + count6)/nSim*100,1)
  res[,3] <- round((count1 + count6)/nSim*100,1)
  res[,4] <- round(count2/nSim*100,1)
  res[,5] <- round(count3/nSim*100,1)
  colnames(res) <- c("ESS", "H0", "H00", "H01", "H02")
  row.names(res) <- "%"
  cat("\n")
  cat("The expected sample size and overall power: \n")
  print(res)
  cat("\n")
  cat("Pr(Early stopping for futility) = ", round(ESF/nSim*100,1), "%.\n")
  cat("Pr(Early stopping for efficacy) = ", round(ESE/nSim*100,1), "%.\n")
  cat("Pr(Enrich subgroup 1) = ", round(trigger01/nSim*100,1), "%.\n")
  cat("Pr(Enrich subgroup 2) = ", round(trigger02/nSim*100,1), "%.\n")

  return(list(nTotal = ceiling(nTotal),
              H00 = round((count1 + count6)/nSim*100,1),
              H01 = round(count2/nSim*100,1),
              H02 = round(count3/nSim*100,1),
              H0 = round((count1 + count2 + count3 + count6)/nSim*100,1),
              ESF = round(ESF/nSim*100,1),
              ESE = round(ESE/nSim*100,1),
              Enrich01 = round(trigger01/nSim*100,1),
              Enrich02 = round(trigger02/nSim*100,1),
              Trigger03 = round(trigger03/nSim*100,1)))
}


#' @title Calculate the \eqn{N2} and the critical value \eqn{C} in the Adaptive
#'   Enrichment Design (Strategy 2) with Sample Size Re-estimation Procedure
#'
#' @description The \code{AED2_SSR.CP()} is used to calculate the sample size required
#'    at the second stage and the critical value used at the final analysis in the
#'    Adaptive Enrichment Design with Sample Size Re-estimation Procedure. In
#'    addition, this function can also used to conduct the conditional power
#'    analysis in terms of \eqn{N2}
#'
#' @param Z1 The test statistic obtained at the interim analysis
#' @param delta The standardized size of treatment effect, which can be estimated
#'   by using \eqn{(\mu_{X} - \mu_{Y})/\sqrt{\sigma^2}}.
#' @param N1 The sample size used at the first stage
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param rho The proportion of subgroup 1
#' @param epsilon The threshold of the difference between subgroup-specific test
#'   statistics.
#' @param alpha The overall Type I error rate
#' @param beta The \code{(1 - Power)}
#' @param N2 The pre-specified sample size used at the second stage, which is used
#'   to conduct the conditional power analysis
#'
#' @return A list contains
#' \itemize{
#'   \item upper.boundary The efficacy stopping boundary
#'   \item lower.boundary The futility stopping boundary
#'   \item N2 The pre-specified sample size used at the second stage, which is
#'           used to implement the conditional power analysis
#'   \item Conditional.Power The value of conditional power given the value of \code{N2} in the
#'           conditional power analysis
#'   \item P.Value The corresponding P-Value used at the final analysis in the conditional
#'           power analysis
#'   \item N2.CP The re-estimated sample size of \code{N2} to ensure an adequate
#'           conditional power
#'   \item c.CP The estimated the critical value used at the final analysis based
#'           the conditional power
#' }
#'
#' @import stats
#'
#' @examples
#' Z1 <- 1.974
#' delta <- 0.355
#' N1 <- 248
#' pstar <- 0.15
#' alpha <- 0.05
#' rho <- 0.5
#' epsilon <- 0.5
#' beta <- 0.20
#' N2 <- 104
#' res <- AED2_SSR.CP(Z1 = Z1, delta = delta, N1 = N1, pstar = pstar,
#'            alpha = alpha, rho = rho, epsilon = epsilon,
#'            beta = beta, N2 = N2)
#'
#' @export
AED2_SSR.CP <- function(Z1 = NULL, delta = NULL, N1 = NULL,
                        pstar, rho, epsilon, alpha, beta, N2 = NULL){
  zp <- qnorm(1 - pstar)
  k <- AED2_SSR.boundary(rho = rho, alpha = alpha, pstar = pstar, epsilon = epsilon)$upper.boundary
  z1 <- Z1
  n1 <- N1*0.5
  A <- 0
  if (z1 >= zp & z1 < k) {
    A <- 1 - pnorm(sqrt(k^2 - z1^2))
  } else if (z1 >= k) {
    A <- 1
  }
  za <- qnorm(1 - A)
  if (!is.null(N2)) {
    cp <- 1 - pnorm(za - sqrt(N2/4)*delta)
    c <- (sqrt(N1)*z1 + sqrt(N2)*za)/(sqrt(N1 + N2))
    p <- 1 - pnorm(c)
    cat("\n")
    cat("Conditional power analysis of AED2-SSR: \n")
    cat("The conditional power = ", round(cp*100,1), "% given N2 = ", N2, ".\n")
    cat("The exact critical value used at the final analysis = ", round(c, 3),
        ", with corresponding P-Value = ", round(p,3), ".\n")
    cat("\n")
  }

  zb <- qnorm(1 - beta)
  N2.CP <- n1/2 * (za + zb)^2/z1^2
  c.CP <- (z1^2 + za*(za + zb))/(sqrt(z1^2 + (za + zb)^2))
  cat("\n")
  cat("The N2 required at the second stage = ", ceiling(N2.CP)*2, ".\n")
  cat("The estimated critical value used at the final analysis = ", round(c.CP,3),
      ".\n")
  cat("\n")
  return(list(upper.boundary = k, lower.boundary = zp,
              N2 = N2, Conditional.Power = cp, P.Value = p,
              N2.CP = 2*ceiling(N2.CP), c.CP = c.CP))
}


#' @title Calculate the critical value used at the final analysis of the
#'   Adaptive Enrichment Design (Strategy 1) with Sample Size Re-estimation Procedure
#'
#' @description The \code{AED1_SSR.boundary()} is used to calculate the critical
#'   value required at the final analysis of the Adaptive Enrichment Design
#'   (Strategy 1) with sample size
#'   re-estimation procedure. In the AED1-SSR design, the adaptive enrichment
#'   strategy is guided by a pre-specified futility stopping boundary and a
#'   threshold of the difference between the subgroup-specific test statistics.
#'
#' @param rho The proportion of subgroup 1.
#' @param alpha The overall Type I error rate.
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param Info The observation information, which is commonly calculated through
#'   the sample size used at each stage of the trial.
#' @param epsilon The threshold of the difference between subgroup-specific test
#'   statistics.
#'
#' @import stats
#'
#' @references
#' \itemize{
#'   \item Lin, R., Yang, Z., Yuan, Y. and Yin, G., 2021. Sample size re-estimation
#'    in adaptive enrichment design. Contemporary Clinical Trials, 100, p.106216.
#'    <doi: 10.1016/j.cct.2020.106216>
#' }
#'
#' @examples
#' AED1_SSR.boundary(rho = 0.5, alpha = 0.05, pstar = 0.2, Info = 0.5, epsilon = 0.5)
#' @export
#'
AED1_SSR.boundary <- function(rho, alpha, pstar, Info, epsilon) {
  N <- 500
  n11 <- rho*N
  n12 <- (1 - rho)*N
  zp <- qnorm(1 - pstar)
  a <- sqrt(n11/(n11 + n12))
  b <- sqrt(n12/(n11 + n12))

  Func <- function(C) {
    f1 <- function(a, b, C, epsilon, z11, z12, Info) {
      z1 <- a*z11 + b*z12
      (1 - pnorm( (C - sqrt(Info)*z1)/sqrt(1 - Info)))*dnorm(z11)*dnorm(z12)
    }
    inte1Func <- function(z12) {
      sapply(z12, function(z12){
        integrate(f1, zp, z12 + epsilon,
                  z12 = z12, epsilon = epsilon, a = a, b = b,
                  C = C, Info = Info)$value
      })
    }
    f12 <- function(a, b, C, epsilon, z11, z12, Info) {
      z1 <- a*z11 + b*z12
      (1 - pnorm( (C - sqrt(Info)*z1)/sqrt(1 - Info)))*dnorm(z12)*dnorm(z11)
    }
    inte2Func <- function(z11) {
      sapply(z11, function(z11) {
        integrate(f12, zp, z11 + epsilon,
                  z11 = z11, epsilon = epsilon, a = a, b = b,
                  C = C, Info = Info)$value
      })
    }

    f2 <- function(z11, z12, Info, epsilon, C) {
      (1 - pnorm( (C - sqrt(Info)*z11)/sqrt(1 - Info)))*dnorm(z11)*dnorm(z12)
    }
    inte3Func <- function(z12) {
      sapply(z12, function(z12){
        integrate(f2, pmin(z12 + epsilon, zp ), C,
                  z12 = z12, Info = Info,
                  epsilon = epsilon, C = C)$value
      })
    }

    f3 <- function(z11, z12, Info, epsilon, C) {
      (1 - pnorm( (C - sqrt(Info)*z12)/sqrt(1 - Info)))*dnorm(z12)*dnorm(z11)
    }
    inte4Func <- function(z11) {
      sapply(z11, function(z11){
        integrate(f3, pmin(z11 + epsilon,zp), C,
                  z11 = z11, Info = Info,
                  epsilon = epsilon, C = C)$value
      })
    }

    (integrate(inte1Func, zp, Inf)$value +
        integrate(inte2Func, zp, Inf)$value +
        integrate(inte3Func, -Inf, zp)$value +
        integrate(inte4Func, -Inf, zp)$value + 2*(1 - pnorm(C))) - alpha
  }
  c <- uniroot(Func, c(zp,4))$root
  return(C = c)
}

#' @title Calculate the sample size required at the second stage of the adaptive
#'   enrichment design (Strategy1) with Sample Size Re-estimation Procedure
#'
#' @description The \code{AED1_SSR.N2()} is used to calculated the sample size
#'   required at the second stage of the Adaptive Enrichment Design (Strategy 1)
#'   with Sample Size Re-estimation Procedure.
#'
#' @param c The critical value used at the final analysis
#' @param z1 The test statistic obtained at the interim analysis
#' @param N1 The sample size used at the first stage
#' @param beta The (1 - power)
#'
#' @return The Value of the re-estimated sample size
#'
#' @references
#' \itemize{
#'   \item Lin, R., Yang, Z., Yuan, Y. and Yin, G., 2021. Sample size re-estimation
#'    in adaptive enrichment design. Contemporary Clinical Trials, 100, p.106216.
#'    <doi: 10.1016/j.cct.2020.106216>
#' }
#'
#' @examples
#' c <- 2.258
#' z1 <- 1.974
#' N1 <- 248
#' beta <- 0.2
#' AED1_SSR.N2(c = c, z1 = z1, N1 = N1, beta = beta)
#'
#' @export
AED1_SSR.N2 <- function(c, z1, N1, beta) {
  func <- function(N2) {
    (c - z1/(sqrt(N1/(N1 + N2))))/sqrt(N2/(N1 + N2)) - qnorm(beta)
  }
  ratio <- uniroot(func, c(0,10000*N1))$root
  return(ratio)
}


#' @title Conduct the simulation studies of the Adaptive Enrichment Design
#'   (Strategy 1) with Sample Size Re-estimation Procedure
#'
#' @description The \code{AED1_SSR.sim()} is used to conduct the simulation study
#'   of the Adaptive Enrichment Design (Strategy 1) with Sample Size Re-estimation
#'   procedure
#'
#' @param N1 The sample size used at the first stage
#' @param rho The proportion of subgroup 1 among the overall patients
#' @param alpha The overall Type I error rate
#' @param beta The (1 - Power)
#' @param pstar The \code{(1 - power)} of accepting the null hypothesis at the
#'    interim analysis.
#' @param theta The sizes of the treatment effect in subgroups 1 and 2 with the
#'   experimental arm
#' @param theta0 The size of the treatment effect in standard arm
#' @param Info The observation information
#' @param K The number of subgroups. The default value is \code{K = 2}
#' @param epsilon The threshold of the difference between the subgroup-specific
#'   test statistic
#' @param sigma0 The variance of the treatment effect
#' @param nSim The number of simulated studies
#' @param Seed The random seed
#'
#' @references
#' \itemize{
#'   \item Lin, R., Yang, Z., Yuan, Y. and Yin, G., 2021. Sample size re-estimation
#'    in adaptive enrichment design. Contemporary Clinical Trials, 100, p.106216.
#'    <doi: 10.1016/j.cct.2020.106216>
#' }
#'
#' @return A list contains
#' \itemize{
#'   \item nTotal The average expected sample size
#'   \item H00 The probability of rejecting the null hypothesis of \eqn{H_{00}}
#'   \item H01 The probability of rejecting the null hypothesis of \eqn{H_{01}}
#'   \item H02 The probability of rejecting the null hypothesis of \eqn{H_{02}}
#'   \item H0  The probabilities of rejecting at least one of the null hypothesis
#'   \item ESF The probability of early stopping for futility
#'   \item ESE The probability of early stopping for efficacy
#'   \item Enrich01 The prevalence of adaptive enrichment of subgroup 1
#'   \item Enrich02 The prevalence of adaptive enrichment of subgroup 2
#' }
#'
#' @examples
#' res <- AED1_SSR.sim(
#'   N1 = 310, rho = 0.5,
#'   alpha = 0.05, beta = 0.2, pstar = 0.2,
#'   theta = c(0,0), theta0 = 0, Info = 0.5,
#'   epsilon = 0.5, sigma0 = 1, nSim = 1000, Seed = 6)
#' @export
AED1_SSR.sim <- function(N1, rho, alpha, beta, pstar,
                         theta, theta0, Info, K = 2, epsilon, sigma0,
                         nSim, Seed) {
  set.seed(Seed)
  cat("\n")
  cat("## Summary of AED1-SSR:    \n")
  cat("\n")
  zalpha <- AED1_SSR.boundary(rho = rho, alpha = alpha, pstar = pstar,
                              Info = Info, epsilon = epsilon)
  c <- round(zalpha,3)
  zb <- pnorm(1-beta)
  n11 <- rho * N1
  n12 <- (1 - rho) * N1
  n1 <- 0.5*N1

  count1 <- count2 <- count3 <-  nTotal <- 0
  Trigger01 <- Trigger02 <- ESE <- ESF <- 0
  for (i in 1:nSim) {
    x1 <- rnorm(n11/2, theta[1], sigma0)
    y1 <- rnorm(n11/2, theta0, sigma0)
    x2 <- rnorm(n12/2, theta[2], sigma0)
    y2 <- rnorm(n12/2, theta0, sigma0)
    x <- c(x1, x2)
    y <- c(y1, y2)
    sigma1 <- var(c(x,y))
    z1 <- (mean(x) - mean(y))/sqrt(2 * sigma1/n1)
    z11 <- (mean(x1) - mean(y1))/sqrt(4 * sigma1/n11)
    z12 <- (mean(x2) - mean(y2))/sqrt(4 * sigma1/n12)

    nTotal <- nTotal + 2 * n1/nSim
    opt <- which.max(c(z11, z12))
    if (opt == 1) {
      if (z11 < zb) {
        ESF <- ESF + 1/nSim
      } else if (z11 >= zb & (z11 - z12) >= epsilon) {
        if (z11 >= c ) {
          ESE <- ESE + 1/nSim
          count1 <- count1 + 1
        } else {
          n2 <- 0.5*AED1_SSR.N2(c = zalpha, z1 = z11, N1 = n11, beta = beta)

          if (is.na(n2)) {
            next
          } else {
            nTotal <- nTotal + 2 * n2/nSim
            x1 <- c(x1, rnorm(ceiling(n2), theta[1], sigma0))
            y1 <- c(y1, rnorm(ceiling(n2), theta0, sigma0))
            sigma1 <- var(c(x1,x2, y1, y2))
            z <- (mean(x1) - mean(y1))/sqrt(2 * sigma1/(n11/2 + n2))
            if (z > c) {
              count1 <- count1 + 1
            }
            Trigger01 <- Trigger01 + 1/nSim
          }
        }


      } else {
        if (z1 >= c) {
          ESE <- ESE + 1/nSim
          count3 <- count3 + 1
        } else {
          n2 <- 0.5*AED1_SSR.N2(c = zalpha, z1 = z1, N1 = N1, beta = beta)

          if (is.na(n2)) {
            next
          } else {
            nTotal <- nTotal + 2 * n2/nSim
            n21 <- 2 * ceiling((n2 * rho))
            n22 <- 2 * ceiling((n2 * (1 - rho)))
            x1 <- c(x1, rnorm(n21/2, theta[1], sigma0))
            y1 <- c(y1, rnorm(n21/2, theta0, sigma0))
            x2 <- c(x2, rnorm(n22/2, theta[2], sigma0))
            y2 <- c(y2, rnorm(n22/2, theta0, sigma0))
            x <- c(x1, x2)
            y <- c(y1, y2)
            sigma1 <- var(c(x, y))
            z <- (mean(x) - mean(y))/sqrt(2 * sigma1/(n1 + n2))
            if (z >= c) {
              count3 <- count3 + 1
            }
          }
        }
      }
    }
    if (opt == 2) {
      if (z12 < zb) {
        ESF <- ESF + 1/nSim
      } else if (z12 >= zb & (z12 - z11) >= epsilon) {
        if (z12 >= c) {
          ESE <- ESE + 1/nSim
          count2 <- count2 + 1
        } else {
          Trigger02 <- Trigger02 + 1/nSim
          n2 <- 0.5*AED1_SSR.N2(c = zalpha, z1 = z12, N1 = n12, beta = beta)

          if (is.na(n2)) {
            next
          } else {
            nTotal <- nTotal + 2 * n2/nSim
            x2 <- c(x2, rnorm(ceiling(n2), theta[2], sigma0))
            y2 <- c(y2, rnorm(ceiling(n2), theta0, sigma0))
            sigma1 <- var(c(x1,x2, y1, y2))
            z <- (mean(x2) - mean(y2))/sqrt(2 * sigma1/(n12/2 + n2))
            if (z > c) {
              count2 <- count2 + 1
            }
          }

        }

      } else {
        if (z1 >= c) {
          ESE <- ESE + 1/nSim
          count3 <- count3 + 1
        } else {
          n2 <- 0.5*AED1_SSR.N2(c = zalpha, z1 = z1, N1 = N1, beta = beta)

          if (is.na(n2)) {
            next
          } else {
            nTotal <- nTotal + 2 * n2/nSim
            n21 <- 2 * ceiling(n2 * rho)
            n22 <- 2 * ceiling(n2 * (1 - rho))

            x1 <- c(x1, rnorm(n21/2, theta[1], sigma0))
            y1 <- c(y1, rnorm(n21/2, theta0, sigma0))
            x2 <- c(x2, rnorm(n22/2, theta[2], sigma0))
            y2 <- c(y2, rnorm(n22/2, theta0, sigma0))
            x <- c(x1, x2)
            y <- c(y1, y2)
            sigma1 <- var(c(x, y))
            z <- (mean(x) - mean(y))/sqrt(2 * sigma1/(n1 + n2))
            if (z >= c) {
              count3 <- count3 + 1
            }
          }

        }
      }
    }
  } ## nSim
  res <- matrix(NA, ncol = 5, nrow = 1)
  res[,1] <- ceiling(nTotal)
  res[,2] <- round((count1 + count2 + count3)/nSim*100,1)
  res[,3] <- round(count3/nSim*100,1)
  res[,4] <- round(count1/nSim*100,1)
  res[,5] <- round(count2/nSim*100,1)
  colnames(res) <- c("ESS", "H0", "H00", "H01", "H02")
  row.names(res) <- "%"
  cat("\n")
  cat("The expected sample size and overall power: \n")
  print(res)
  cat("\n")
  cat("Pr(Early stopping for futility) = ", round(ESF*100,1), "%.\n")
  cat("Pr(Early stopping for efficacy) = ", round(ESE*100,1), "%.\n")
  cat("Pr(Enrich subgroup 1) = ", round(Trigger01*100,1), "%.\n")
  cat("Pr(Enrich subgroup 2) = ", round(Trigger02*100,1), "%.\n")
  cat("\n")
  return(list(nTotal = ceiling(nTotal),
              H00 = round(count3/nSim*100,1),
              H01 = round(count1/nSim*100,1),
              H02 = round(count2/nSim*100,1),
              H0  = round((count1 + count2 + count3)/nSim*100,1),
              ESF = round(ESF*100,1), ESE = round(ESE*100,1),
              Enrich01 = round(Trigger01*100,1),
              Enrich02 = round(Trigger02*100,1)))
}

#' @title Calculate the conditional power of the Adaptive Enrichment Design with
#'   (Strategy 1) Sample Size Re-estimation Procedure
#'
#' @description The \code{AED1_SSR.CP()} is used to calculate the conditional
#'    power of the Adaptive Enrichment Design (Strategy 1) with sample size
#'    re-estimation procedure
#'
#' @param c The critical value used at the final analysis
#' @param Z1 The test statistic obtained at the interim analysis
#' @param N1 The sample size used at the first stage
#' @param N2 The sample size used at the second stage
#'
#' @return A list contains
#' \itemize{
#'   \item Critical.Value The critical value used at the final analysis
#'   \item Conditional.Power The value of conditional power given the observed data
#' }
#'
#' @references
#' \itemize{
#'   \item Lin, R., Yang, Z., Yuan, Y. and Yin, G., 2021. Sample size re-estimation
#'    in adaptive enrichment design. Contemporary Clinical Trials, 100, p.106216.
#'    <doi: 10.1016/j.cct.2020.106216>
#' }
#'
#' @examples
#' c <- 2.258
#' Z1 <- 1.975
#' N1 <- 248
#' N2 <- 200
#' AED1_SSR.CP(c = 2.258, Z1 = 1.974, N1 = 248, N2 = 200)
#' @export
AED1_SSR.CP <- function(c, Z1, N1, N2) {
  t <- N1/(N1 + N2)
  cp <- 1 - pnorm((c - Z1/sqrt(t))/(sqrt(1 - t)))
  return(list(
    Critical.Value = c,
    Conditional.Power = cp))
}


#' @title Conduct the simulation studies of the Marker Sequential Test design
#'
#' @description The \code{MaST.sim()} is used to conduct the simulation studies
#'   of the marker sequential test design (MaST).
#'
#' @param N The total sample size used at the trial
#' @param rho The proportion of subgroup 1 among the overall patients
#' @param alpha The overall Type I error rate
#' @param beta The (1 - Power)
#' @param theta The sizes of treatment effect in subgroups 1 and 2 with the
#'    experimental arm
#' @param theta0 The size of treatment effect in the standard arm
#' @param sigma0 The variance of the treatment effect
#' @param nSim The number of simulated studies
#' @param Seed The random seed
#'
#' @return A list contains
#' \itemize{
#'   \item nTotal The average expected sample size
#'   \item H00 The probability of rejecting the null hypothesis of \eqn{H_{00}}
#'   \item H01 The probability of rejecting the null hypothesis of \eqn{H_{01}}
#'   \item H02 The probability of rejecting the null hypothesis of \eqn{H_{02}}
#'   \item H0  The probabilities of rejecting at least one of the null hypothesis
#' }
#'
#' @references
#' \itemize{
#'   \item Freidlin, B., Korn, E. L., and Gray, R. (2014). Marker sequential
#'   test (MaST) design. Clinical trials, 11(1), 19-27. <doi:10.1177/1740774513503739>
#' }
#'
#' @examples
#' N <- 310
#' rho <- 0.5
#' alpha <- 0.05
#' beta <- 0.20
#' theta <- c(0,0)
#' theta0 <- 0
#' sigma0 <- 1
#' nSim <- 1000
#' Seed <- 6
#' MaST.sim(N = N, rho = rho, alpha = alpha, beta = beta,
#'          theta = theta, theta0 = theta0, sigma0 = sigma0,
#'          nSim = nSim, Seed = Seed)
#' @export
MaST.sim <- function(N, rho, alpha, beta,
                     theta, theta0, sigma0, nSim, Seed) {
  set.seed(Seed)
  c <- qnorm(1 - alpha)
  if (alpha >= 0.04) {
    alpha1 <- 0.04;         c1 <- qnorm(1 - alpha1)
    alpha2 <- alpha - 0.04; c2 <- qnorm(1 - alpha2)
  } else if (alpha >= 0.025 & alpha < 0.04) {
    alpha1 <- 0.022;        c1 <- qnorm(1 - alpha1)
    alpha2 <- alpha - 0.022;c2 <- qnorm(1 - alpha2)
  }

  count1 <- count2 <- count3 <- nTotal <- 0
  for (i in 1:nSim) {
    n11 <- N*rho
    n12 <- N*(1 - rho)
    n1 <- N/2
    nTotal <- nTotal + N
    x1 <- rnorm(n11/2, theta[1], sigma0)
    y1 <- rnorm(n11/2, theta0, sigma0)
    x2 <- rnorm(n12/2, theta[2], sigma0)
    y2 <- rnorm(n12/2, theta0, sigma0)
    sigma1 <- var(c(x1,x2,y1,y2))
    z1 <- (mean(c(x1,x2)) - mean(c(y1,y2)))/sqrt(2 * sigma1/n1); #print(z1)
    z11 <- (mean(x1) - mean(y1))/sqrt(4 * sigma1/n11); #print(z11) ## Test statistics for subgroup 1
    z12 <- (mean(x2) - mean(y2))/sqrt(4 * sigma1/n12); #print(z12) ## Test statistics for subgroup 2

    if (z11 > c1) {
      if (z12 > c) {
        count3 <- count3 + 1
      } else {
        count1 <- count1 + 1
      }
    } else {
      if (z1 > c2) {
        count2 <- count2 + 1
      }
    }

  }
  return(list(nTotal = ceiling(nTotal/nSim),
              H01 = round(count1/nSim*100,1),
              H02 = round(count2/nSim*100,1),
              H00 = round(count3/nSim*100,1),
              H0  = round((count1 + count2 + count3)/nSim*100,1)))
}


#' @title Calculate the critical value used at the final analysis in AED
#'
#' @description \code{AED.boundary()} is used to calculate the critical value
#'     used at the final analysis in AED design, meanwhile preserving the overall
#'     type I error rate at \eqn{\alpha} level
#'
#' @param rho The proportion of subgroup 1
#' @param alpha The overall type I error rate
#' @param Info The infromation fraction
#' @param epsilon The threshold of difference between the subgroup-specific test
#'   statistics
#'
#' @return The critical value used at the final analysis
#'
#' @references
#' \itemize{
#'   \item Lin, R., Yang, Z., Yuan, Y. and Yin, G., 2021. Sample size re-estimation
#'    in adaptive enrichment design. Contemporary Clinical Trials, 100, p.106216.
#'    <doi: 10.1016/j.cct.2020.106216>
#' }
#'
#' @examples
#' AED.boundary(rho = 0.5, alpha = 0.05, Info = 0.5, epsilon = 0.5)
#' @export
AED.boundary <- function(rho, alpha, Info, epsilon) {
  N <- 500
  n11 <- N * rho  ## total sample size used in stage 1
  n12 <- N * (1 - rho) ## total sample size used to subgroup 2 in stage 2
  #zp <- qnorm(1 - pstar)    ## quantile of power (conditional power)
  a <- sqrt(n11/(n11 + n12))  ## weight for subgroup 1
  b <- sqrt(n12/(n11 + n12))  ## weight for subgroup 2
  ##
  Func <- function(C) {
    f1 <- function(C, z11, z12, a, b, Info, epsilon) {
      z1 <- a*z11 + b*z12
      (1 - pnorm((C - sqrt(Info)*z1)/sqrt(1 - Info)))*dnorm(z11)*dnorm(z12)
    }
    inte1Func <- function(z12) {
      sapply(z12, function(z12) {
        integrate(f1, z12 - epsilon, z12 + epsilon,
                  a = a, b = b, epsilon = epsilon, z12 = z12,
                  C = C, Info = Info)$value
      })
    }

    f3 <- function(z11, z12, C, Info, epsilon) {
      (1 - pnorm((C - sqrt(Info)*z11)/sqrt(1 - Info)))*dnorm(z11)*dnorm(z12)
    }
    inte3Func <- function(z12) {
      sapply(z12, function(z12){
        integrate(f3, z12 + epsilon, Inf,
                  z12 = z12, epsilon = epsilon,
                  Info = Info, C = C)$value
      })
    }
    f4 <- function(z11, z12, C, Info, epsilon) {
      (1 - pnorm((C - sqrt(Info)*z12)/sqrt(1 - Info)))*dnorm(z12)*dnorm(z11)
    }
    inte4Func <- function(z11) {
      sapply(z11, function(z11){
        integrate(f4, z11 + epsilon,Inf, epsilon = epsilon,
                  C = C, Info = Info, z11 = z11)$value
      })
    }
    integrate(inte1Func, -Inf, Inf)$value +
      integrate(inte3Func, -Inf, Inf)$value +
      integrate(inte4Func, -Inf, Inf)$value  - alpha
  }

  c <- uniroot(Func, c(-400,400))$root
  return(C = c)
}


#' @title Conduct the simulation studies of the Adaptive Enrichment Design without
#'   early stopping boundary
#'
#' @description The \code{AED.sim()} is used to conduct the simulation studies
#'    of the Adaptive Enrichment Design without early stopping boundary. The AED
#'    design is quite similar with the AED1_SSR design. But, in the AED design,
#'    the futility stopping boundary and the Sample Size Re-estimation Procedure
#'    are removed. On the contrary, a fixed sample size is used to replace the
#'    sample size re-estimated procedure. In addition, an \eqn{\epsilon}-rule is
#'    also introduced to select the subgroup with larger subgroup-specific test
#'    statistic.
#'
#' @param N1 The sample size used at the first stage
#' @param N2 The sample size used at the second stage
#' @param rho The proportion of the subgroup 1
#' @param alpha The overall Type I error rate
#' @param beta The (1 - Power)
#' @param theta The sizes of treatment effects in subgroups 1 and 2 among the
#'     experimental arm
#' @param theta0 The size of treatment effect in standard arm
#' @param K The number of subgroups
#' @param Info The observed information
#' @param epsilon The threshold of difference between the subgroup-specific test
#'   statistics
#' @param sigma0 The variance of the treatment effect
#' @param nSim The number of simulated studies
#' @param Seed The random Seed
#'
#' @return A list contains
#' \itemize{
#'   \item nTotal The average expected sample size
#'   \item H00 The probability of rejecting the null hypothesis of \eqn{H_{00}}
#'   \item H01 The probability of rejecting the null hypothesis of \eqn{H_{01}}
#'   \item H02 The probability of rejecting the null hypothesis of \eqn{H_{02}}
#'   \item H0  The probabilities of rejecting at least one of the null hypothesis
#'   \item Enrich01 The prevalence of adaptive enrichment of subgroup 1
#'   \item Enrich02 The prevalence of adaptive enrichment of subgroup 2
#' }
#'
#' @references
#' \itemize{
#'   \item Lin, R., Yang, Z., Yuan, Y. and Yin, G., 2021. Sample size re-estimation
#'    in adaptive enrichment design. Contemporary Clinical Trials, 100, p.106216.
#'    <doi: 10.1016/j.cct.2020.106216>
#' }
#'
#' @examples
#' N1 <- 310
#' N2 <- 310
#' rho <- 0.5
#' alpha <- 0.05
#' beta <- 0.20
#' theta <- c(0,0)
#' theta0 <- 0
#' K <- 2
#' Info <- 0.5
#' epsilon <- 0.5
#' sigma0 <- 1
#' nSim <- 1000
#' Seed <- 6
#' AED.sim(N1 = N1, N2 = N2, rho = rho, alpha = alpha,
#'         beta = beta, theta = theta, theta0 = theta0,
#'         K  = K, Info = Info, epsilon = epsilon,
#'         sigma0 = sigma0, nSim = nSim, Seed = Seed)
#' @export
AED.sim <- function(N1, N2, rho, alpha, beta,
                    theta, theta0, K, Info,epsilon, sigma0,
                    nSim, Seed) {
  set.seed(Seed)
  c <- round(AED.boundary(rho = rho, alpha = alpha, Info = Info, epsilon = epsilon),3)
  n11 <- rho * N1
  n12 <- (1 - rho) * N1
  n1 <- 0.5*N1

  count1 <- count2 <- count3 <-  nTotal <- 0
  Trigger01 <- Trigger02 <- 0
  for (i in 1:nSim) {
    x1 <- rnorm(n11/2, theta[1], sigma0)
    y1 <- rnorm(n11/2, theta0, sigma0)
    x2 <- rnorm(n12/2, theta[2], sigma0)
    y2 <- rnorm(n12/2, theta0, sigma0)
    x <- c(x1, x2)
    y <- c(y1, y2)
    sigma1 <- var(c(x,y))
    z11 <- (mean(x1) - mean(y1))/sqrt(4 * sigma1/n11); #print(z11) ## Test statistics for subgroup 1
    z12 <- (mean(x2) - mean(y2))/sqrt(4 * sigma1/n12); #print(z12) ## Test statistics for subgroup 2

    nTotal <- nTotal + 2 * n1/nSim
    opt <- which.max(c(z11, z12))
    if (opt == 1) {
      if ((z11 - z12) >= epsilon) {
        n2 <- rho*N2*0.5
        if (is.na(n2)) {
          next
        } else {
          nTotal <- nTotal + 2 * n2/nSim
          x1 <- c(x1, rnorm(ceiling(n2), theta[1], sigma0))
          y1 <- c(y1, rnorm(ceiling(n2), theta0, sigma0))
          sigma1 <- var(c(x1,x2, y1, y2))
          z <- (mean(x1) - mean(y1))/sqrt(2 * sigma1/(n11/2 + n2))
          if (z > c) {
            count1 <- count1 + 1
          }
          Trigger01 <- Trigger01 + 1/nSim
        }
      } else {
        n2 <- 0.5*N2
        if (is.na(n2)) {
          next
        } else {
          nTotal <- nTotal + 2 * n2/nSim
          n21 <- 2 * ceiling((n2 * rho))
          n22 <- 2 * ceiling((n2 * (1 - rho)))
          x1 <- c(x1, rnorm(n21/2, theta[1], sigma0))
          y1 <- c(y1, rnorm(n21/2, theta0, sigma0))
          x2 <- c(x2, rnorm(n22/2, theta[2], sigma0))
          y2 <- c(y2, rnorm(n22/2, theta0, sigma0))
          x <- c(x1, x2)
          y <- c(y1, y2)
          sigma1 <- var(c(x, y))
          z <- (mean(x) - mean(y))/sqrt(2 * sigma1/(n1 + n2))
          if (z >= c) {
            count3 <- count3 + 1
          }
        }
      }
    }

    if (opt == 2) {
      if ((z12 - z11) >= epsilon) {
        n2 <- (1-rho)*N2*0.5
        if (is.na(n2)) {
          next
        } else {
          nTotal <- nTotal + 2 * n2/nSim
          x2 <- c(x2, rnorm(ceiling(n2), theta[2], sigma0))
          y2 <- c(y2, rnorm(ceiling(n2), theta0, sigma0))
          sigma1 <- var(c(x1,x2, y1, y2))
          z <- (mean(x2) - mean(y2))/sqrt(2 * sigma1/(n12/2 + n2))
          if (z > c) {
            count2 <- count2 + 1
          }
          Trigger02 <- Trigger02 + 1/nSim
        }
      } else {
        n2 <- N2*0.5
        if (is.na(n2)) {
          next
        } else {
          nTotal <- nTotal + 2 * n2/nSim
          n21 <- 2 * ceiling((n2 * rho))
          n22 <- 2 * ceiling((n2 * (1 - rho)))
          x1 <- c(x1, rnorm(n21/2, theta[1], sigma0))
          y1 <- c(y1, rnorm(n21/2, theta0, sigma0))
          x2 <- c(x2, rnorm(n22/2, theta[2], sigma0))
          y2 <- c(y2, rnorm(n22/2, theta0, sigma0))
          x <- c(x1, x2)
          y <- c(y1, y2)
          sigma1 <- var(c(x, y))
          z <- (mean(x) - mean(y))/sqrt(2 * sigma1/(n1 + n2))
          if (z >= c) {
            count3 <- count3 + 1
          }
        }
      }
    }



  } ## nSim
  return(list(nTotal = ceiling(nTotal),
              H00 = round(count3/nSim*100,1),
              H01 = round(count1/nSim*100,1),
              H02 = round(count2/nSim*100,1),
              H0  = round((count1 + count2 + count3)/nSim*100,1),
              Enrich01 = round(Trigger01*100,1),
              Enrich02 = round(Trigger02*100,1)
  )
  )
}

