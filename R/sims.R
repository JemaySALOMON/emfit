simulPhenos <- R6::R6Class(
  "simulPhenos",
  public = list(
    panelGenos = NULL,
    K = NULL,
    nbBlocks = NULL,
    nbSimuls = NULL,
    mainSeed = NULL,
    simulSeeds = NULL,
    simulIds = NULL,
    levSpecies = NULL,
    fEffs = NULL,
    sigma2_u = NULL,
    sigma2_err = NULL,
    h2 = NULL,
    sims = list(),
    contrasts = "contr.sum",
    initialize = function(panelGenos,
                          nbBlocks,
                          h2,
                          sigma_u,
                          K,
                          nbSimuls = 2,
                          mainSeed = 123,
                          fEffs) {

      self$panelGenos <- panelGenos
      self$nbBlocks <- nbBlocks
      self$nbSimuls <- nbSimuls
      self$mainSeed <- mainSeed
      self$fEffs <- fEffs
      self$h2 <- h2
      self$K = K
      self$sigma2_u <- sigma2_u
      set.seed(mainSeed)
      self$simulSeeds <- sample.int(1e6, nbSimuls)
      self$simulIds <- sprintf("simul%03d", seq_len(nbSimuls))
      self$sigma2_err <- ((1 - h2) / h2) * (nbBlocks * self$sigma2_u)

    },
    simulateOne = function(simulId) {
      stopifnot(simulId %in% seq_len(self$nbSimuls))
      set.seed(self$simulSeeds[simulId])
      dat <- data.frame(
        id  = rep(self$panelGenos, times = self$nbBlocks),
        block  = factor(rep(seq_len(self$nbBlocks),
                            each = length(self$panelGenos)))
      )
      X <- model.matrix(~ 1 + block, dat,
                        contrasts.arg = list(block =contrasts))
      Z <- model.matrix(~ 0 + id, dat)
      colnames(Z) <- self$panelGenos
      # simulate u=>using K
      Q <- self$sigma2_u * self$K
      u =  MASS::mvrnorm(n=1, mu=rep(0,length(panelGenos)), Sigma=Q)
      names(u) <- self$panelGenos
      R <- diag(self$sigma2_err, nrow(dat))
      e <- MASS::mvrnorm(1, rep(0, nrow(dat)), R)
      stopifnot(identical(names(u), colnames(Z)))
      y <- X %*% fEffs + Z %*% u + e
      dat$yield <- as.numeric(y)
      list(
        data = dat,
        y = dat$yield,
        X = X,
        Z = Z,
        u = u
      )
    },
    runSimul = function() {
      for (i in seq_len(self$nbSimuls)) {
        self$sims[[self$simulIds[i]]] <- self$simulateOne(i)
      }
      invisible(self)
    }

  )
)
