
#' EM algorithm for Maximum Likelihood estimation in Linear Mixed Models
#'
#' Fits the model: $$y = X\beta + Zu + e$$
#' where
#' $$u \sim N(0, \sigma_u^2 K), \quad e \sim N(0, \sigma_e^2 I_n)$$
#' using the EM algorithm for ML (not REML).
#'
#' @param y Response vector (n × 1)
#' @param X Design matrix for fixed effects (n × p)
#' @param Z Design matrix for random effects (n × q)
#' @param K Covariance matrix for random effects (q × q)
#' @param niter Maximum number of EM iterations (default 100)
#' @param tol Convergence threshold for relative change of variance components (default 1e-6)
#' @param silent Logical, if TRUE suppress iteration messages (default TRUE)
#' @param verbose Logical, if TRUE prints convergence message (default TRUE)
#'
#' @return A list containing:
#' \item{beta}{Estimated fixed effects vector $$\hat{\beta}$$}
#' \item{u}{Estimated random effects vector $$\hat{u}$$ (BLUP)}
#' \item{Vu}{Estimated variance of random effects $$\sigma_u^2$$}
#' \item{Ve}{Estimated residual variance $$\sigma_e^2$$}
#' \item{loglik}{Vector of log-likelihood values across iterations}
#' \item{iterations}{Number of iterations performed}
#'
#' @references
#' Henderson, C. R. (1975).
#' Best Linear Unbiased Estimation and Prediction under a Selection Model.
#' \emph{Biometrics}, 31(2), 423–447.
#'
#' Gumedze, F. N., & Dunne, T. T. (2011).
#' Parameter estimation and inference in the linear mixed model.
#' \emph{Statistical Methods in Medical Research}, 20(4), 397–415..
#'
#' @author Jemay SALOMON
EM_ML <- function(y, X, Z, K,
                  niter = 100,
                  tol = 1e-6,
                  silent = TRUE,
                  verbose = TRUE) {


  y <- as.matrix(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  K <- as.matrix(K)

  n <- nrow(y)
  p <- ncol(X)
  q <- ncol(Z)


  ZK <- Z %*% K
  ZKZt <- tcrossprod(ZK, Z)
  K_inv <- solve(K)
  beta <- solve(crossprod(X)) %*% crossprod(X, y)

  # Initial residuals
  e <- as.numeric(y - X %*% beta)
  Ve <- var(e) * 0.7
  Vu <- var(e) * 0.3
  loglik <- numeric(niter)


  for (iter in 1:niter) {
    # V^(t) = σ_u^2 Z K Z' + σ_e^2 I_n
    V <- Vu * ZKZt + diag(Ve, n)
    Vinv <- solve(V)

    # GLS estimator for β:
    # β̂^(t) = (X' V^{-1} X)^{-1} X' V^{-1} y
    XtVinv <- crossprod(X, Vinv)
    XtVinvX <- XtVinv %*% X
    XtVinvX_inv <- solve(XtVinvX)
    beta <- XtVinvX_inv %*% XtVinv %*% y

    # BLUP for u:
    # û^(t) = σ_u^2 K Z' V^{-1} (y - X β̂^(t))
    u <- Vu * K %*% t(Z) %*% Vinv %*% (y - X %*% beta)

    # Conditional variance of u:
    # Var(u|y) = σ_u^2 K - σ_u^4 K Z' V^{-1} Z K
    ZtVinvZ <- t(Z) %*% Vinv %*% Z
    Varu <- Vu * K - Vu^2 * K %*% ZtVinvZ %*% K
    e_hat <- y - X %*% beta - Z %*% u

    # σ_e^2 update (ML):
    # σ_e^2^(t+1) = (e_hat' e_hat + tr(Z Var(u|y) Z')) / n
    tr_e <- sum(diag(Z %*% Varu %*% t(Z)))
    Ve_new <- as.numeric((t(e_hat) %*% e_hat + tr_e) / n)

    # σ_u^2 update (ML):
    # σ_u^2^(t+1) = (û' K^{-1} û + tr(K^{-1} Var(u|y))) / q
    tr_u <- sum(diag(K_inv %*% Varu))
    Vu_new <- as.numeric((t(u) %*% K_inv %*% u + tr_u) / q)

    Ve_new <- max(Ve_new, 1e-10)
    Vu_new <- max(Vu_new, 1e-10)

    V_new <- Vu_new * ZKZt + diag(Ve_new, n)
    cholV <- chol(V_new)
    logdetV <- 2 * sum(log(diag(cholV)))
    r <- y - X %*% beta
    loglik[iter] <- as.numeric(
      -0.5 * (n * log(2 * pi) + logdetV + t(r) %*% solve(V_new) %*% r)
    )
    delta <- max(
      abs(Ve_new - Ve) / Ve,
      abs(Vu_new - Vu) / Vu
    )

    Ve <- Ve_new
    Vu <- Vu_new

    if (!silent && iter %% 10 == 0) {
      cat(sprintf(
        "Iter %3d | σ²u = %.4f | σ²e = %.4f | logL = %.2f\n",
        iter, Vu, Ve, loglik[iter]
      ))
    }

    if (delta < tol && iter > 5) {
      loglik <- loglik[1:iter]
      if (verbose)
        cat(sprintf("\nConverged at iteration %d\n", iter))
      break
    }
  }
  V <- Vu * ZKZt + diag(Ve, n)
  Vinv <- solve(V)
  XtVinv <- t(X) %*% Vinv
  XtVinvX <- XtVinv %*% X
  XtVinvX_inv <- solve(XtVinvX)
  beta <- XtVinvX_inv %*% XtVinv %*% y
  u <- Vu * K %*% t(Z) %*% Vinv %*% (y - X %*% beta)

  return(list(
    beta = beta,
    u = u,
    blue =  as.numeric(beta[1])+u,
    Vu = Vu,
    Ve = Ve,
    loglik = loglik,
    iterations = length(loglik)
  ))
}
