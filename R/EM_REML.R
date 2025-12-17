#' EM–REML for Linear Mixed Models
#'
#' Estimate variance components in a linear mixed model using
#' the Expectation–Maximization (EM) algorithm under
#' Restricted Maximum Likelihood (REML).
#'
#' The model is:
#' $$
#' y = X\beta + Zu + e
#' $$
#'
#' with
#' $$
#' u \sim \mathcal{N}(0, \sigma_u^2 K), \qquad
#' e \sim \mathcal{N}(0, \sigma_e^2 I_n)
#' $$
#'
#' Estimation is performed using Henderson’s Mixed Model Equations (MME),
#' avoiding explicit inversion of the marginal covariance matrix.
#'
#' @param y Numeric vector or matrix of observations (n × 1).
#' @param X Fixed-effects design matrix (n × p).
#' @param Z Random-effects design matrix (n × q).
#' @param K Relationship (covariance) matrix for random effects (q × q).
#' @param initVe Optional initial value for residual variance \eqn{\sigma_e^2}.
#'   Defaults to half the empirical variance of \code{y}.
#' @param initVu Optional initial value for random-effect variance
#'   \eqn{\sigma_u^2}. Defaults to half the empirical variance of \code{y}.
#' @param silent Logical. If \code{FALSE}, intermediate estimates are printed.
#' @param tol Convergence tolerance on relative change of variance components.
#' @param niter Maximum number of EM iterations.
#' @param verbose Logical. If \code{TRUE}, prints a convergence message.
#'
#' @details
#' The EM algorithm alternates between:
#'
#' \strong{E-step}
#' \itemize{
#'   \item Solve Henderson's MME to obtain the BLUE of \eqn{\beta}
#'         and BLUP of \eqn{u}
#'   \item Compute the conditional variance
#'         \eqn{\mathrm{Var}(u \mid y)}
#' }
#'
#' \strong{M-step}
#' \itemize{
#'   \item Update residual variance:
#'   \eqn{\sigma_e^2 = E(e'e \mid y)/(n - p)}
#'   \item Update random-effect variance:
#'   \eqn{\sigma_u^2 =
#'   [\hat{u}'K^{-1}\hat{u} + \sigma_e^2 \mathrm{tr}(K^{-1}C_{22})]/\mathrm{rank}(K)}
#' }
#'
#' REML corrects for the loss of degrees of freedom due to estimation of
#' fixed effects.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{beta}: Estimated fixed effects (BLUE).
#'   \item \code{u}: Estimated random effects (BLUP).
#'   \item \code{Vu}: Estimated random-effect variance \eqn{\sigma_u^2}.
#'   \item \code{Ve}: Estimated residual variance \eqn{\sigma_e^2}.
#'   \item \code{iterations}: Number of EM iterations performed.
#' }
#'
#' @references
#' Henderson, C. R. (1975).
#' Best Linear Unbiased Estimation and Prediction under a Selection Model.
#' \emph{Biometrics}, 31(2), 423–447.
#'
#' Gumedze, F. N., & Dunne, T. T. (2011).
#' Parameter estimation and inference in the linear mixed model.
#' \emph{Statistical Methods in Medical Research}, 20(4), 397–415.
#'
#' @author Jemay SALOMON

EM_REML <- function(y, X, Z, K,
                    initVe = NULL,
                    initVu = NULL,
                    silent = TRUE,
                    tol = 1e-6,
                    niter = 200,
                    verbose = TRUE) {

  y <- as.matrix(y)   # (n × 1)
  X <- as.matrix(X)   # (n × p)
  Z <- as.matrix(Z)   # (n × q)
  K <- as.matrix(K)   # (q × q)

  n <- nrow(y)
  p <- qr(X)$rank
  q <- ncol(Z)
  rankK <- qr(K)$rank # rang effectif de K

  Kinv <- solve(K + diag(1e-8, q))

  if (is.null(initVe)) initVe <- var(y) * 0.5
  if (is.null(initVu)) initVu <- var(y) * 0.5

  Ve <- as.numeric(initVe)  # σ²_e
  Vu <- as.numeric(initVu)  # σ²_u

  XtX <- crossprod(X)
  XtZ <- crossprod(X, Z)
  ZtX <- crossprod(Z, X)
  ZtZ <- crossprod(Z, Z)
  Xty <- crossprod(X, y)
  Zty <- crossprod(Z, y)

  for (iter in 1:niter) {

    lambda <- Ve / Vu

    lhs <- rbind(
      cbind(XtX, XtZ),
      cbind(ZtX, ZtZ + lambda * Kinv)
    )
    rhs <- rbind(Xty, Zty)
    sol <- solve(lhs, rhs)

    beta <- sol[1:p]           # BLUE de β
    u    <- sol[(p + 1):(p + q)]  # BLUP de u
    e <- y - X %*% beta - Z %*% u
    lhs_inv <- solve(lhs)
    C22 <- lhs_inv[(p + 1):(p + q), (p + 1):(p + q)]
    Ve_new <- as.numeric(t(e) %*% y) / (n - p)
    Vu_new <- (
      t(u) %*% Kinv %*% u +
        Ve * sum(Kinv * C22)
    ) / rankK

    Ve_new <- max(Ve_new, 1e-10)
    Vu_new <- max(Vu_new, 1e-10)
    delta <- max(
      abs(Ve_new - Ve) / Ve,
      abs(Vu_new - Vu) / Vu
    )

    Ve <- Ve_new
    Vu <- Vu_new

    if (!silent && iter %% 10 == 0) {
      cat(sprintf(
        "Iter %3d | σ²u = %.4f | σ²e = %.4f | h² = %.3f\n",
        iter, Vu, Ve
      ))
    }

    if (delta < tol && iter > 5) {
      if (verbose) cat("Converged\n")
      break
    }
  }
  y_hat <- X %*% beta + Z %*% u
  genos <- colnames(Z)
  y_hat_genos <- tapply(as.numeric(y_hat), c(genos, genos), mean)
  list(
    beta = beta,
    blue = y_hat_genos, # BLUE
    u = u,         # BLUP
    Vu = Vu,       # σ²_u (REML)
    Ve = Ve,       # σ²_e (REML)
    iterations = iter
  )
}
