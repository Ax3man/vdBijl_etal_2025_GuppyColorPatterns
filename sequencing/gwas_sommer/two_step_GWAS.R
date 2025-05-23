## This function is a modified version of sommer::GWAS
# It makes two changes:
# 1. It can now estimate the variance components based on a larger population than is present in the
#    marker matrix. It is assumed one uses a combined GRM and A matrix (the H matrix).
# 2. It can combine multiple random effect matrices into 1 for the marker evaluation.

# for testing
# fixed = car_6 ~ 1
# random = ~ vsr(fish_idA, Gu = H2) + vsr(fish_idX, Gu = X2) + vsr(dsr(patriline))
# rcov = ~ units
# #gTerm = "u:fish_idA"
# data = ornaments
# M = genomat[, 1:10]
# min.MAF = 0.1
# getPEV = FALSE
#
# weights <- NULL
# W <- NULL
# nIters=20; tolParConvLL = 1e-03; tolParInv = 1e-06
# init=NULL; constraints=NULL; method="NR"
#
# naMethodX="exclude"
# naMethodY="exclude"
# returnParam=FALSE
# dateWarning=TRUE; date.warning=TRUE
# verbose=TRUE
# stepWeight=NULL; emWeight=NULL
# n.PC = 0;
# P3D = TRUE

# data <- left_join(ornaments, dplyr::select(pd, fish_id, generation) |> distinct()) |>
#  filter(generation == 'gen_4')
# H3 <- H2[data$fish_id, data$fish_id]
# X3 <- X2[data$fish_id, data$fish_id]
# Y3 <- Y2[data$fish_id, data$fish_id]
#
# two_step_GWAS(
#   fixed = car_6 ~ selection2,
#   random = ~ vsr(fish_idA, Gu = H3) + vsr(fish_idX, Gu = X3) + vsr(fish_idY, Gu = Y3),
#   data = data,
#   M = geno$mat,
#   ind_names = data$fish_id
# )

two_step_GWAS <- function(
    fixed, random, rcov, data,
    M = NULL, min.MAF = 0.05,
    ind_names, # a vector with the names of the individuals, for correct ordering of all the matrices

    nIters = 20, tolParConvLL = 1e-03, tolParInv = 1e-06,
    init = NULL, constraints = NULL, method = "NR",
    getPEV = FALSE,
    naMethodX = "exclude",
    naMethodY = "exclude",
    verbose = TRUE,
    stepWeight = NULL, emWeight = NULL
) {

  if (length(which(is.na(M))) > 0) {
    stop("Please provide an imputed marker matrix (M).", call. = FALSE)
  }
  if (length(setdiff(rownames(M), ind_names)) > 0) {
    stop("Please make sure all rownames of `M` are present in `ind_names`.")
  }
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }

  ## return all parameters for a mixed model (NOT ACTUAL FITTING)
  res <- mmer(
    fixed = fixed, random = random, rcov = rcov, data = data,
    returnParam = TRUE,
    verbose = TRUE
  )

  cat('Estimating genetic variances\n')

  lastmodel <- .Call(
    "_sommer_MNR", PACKAGE = "sommer",
    res$yvar, res$X,res$Gx,res$Z, res$K, res$R, res$GES, res$GESI, res$W, res$isInvW,
    res$nIters, res$tolParConvLL, res$tolParInv, res$selected, res$getPEV, res$verbose,
    TRUE, res$stepWeight, res$emWeight
  )

  cat("Processing matrices\n")

  ## get input matrices
  # the individuals with marker data:
  ind_marker <- levels(factor(rownames(M)))
  # some matrices are sorted by lexical order...
  ind_names2 <- levels(factor(ind_names))

  Y <- scale(res$yvar)
  rownames(Y) <- ind_names2
  Y <- Y[ind_marker, , drop = FALSE]

  # Z is the random effect incidence matrix, and should be the same for the different genetic effects
  Z <- res$Z[[1]]
  rownames(Z) <- ind_names
  # subsetting the full matrix to only have the individuals present in the marker matrix
  Z <- Z[ind_names %in% ind_marker, ind_marker]

  # fixed effect matrix X
  X <- do.call(cbind, res$X) |> make.full()
  rownames(X) <- ind_names
  X <- X[ind_marker, , drop = FALSE]

  # Vinv
  Vinv <- lastmodel$Vi
  rownames(Vinv) <- ind_names2
  colnames(Vinv) <- colnames(res$Z[[1]])
  Vinv <- Vinv[ind_marker, ind_marker]

  if (nrow(M) != ncol(Z)) {
    stop(paste(
      "Marker matrix M needs to have same numbers of rows(", nrow(M),
      ") than columns of the gTerm incidence matrix(",ncol(Z),")."
    ), call. = FALSE)
  }

  cat("Performing GWAS evaluation\n")
  results <- .Call("_sommer_gwasForLoop",PACKAGE = "sommer",
                     M,Y,as.matrix(Z),X,Vinv,min.MAF,TRUE
  )
  preScores <- matrix(results[, , 1], nrow(results), ncol(results))
  bs <- matrix(results[, , 2], nrow(results), ncol(results))
  bs.se <- matrix(results[, , 3], nrow(results), ncol(results))

  v2 <- length(Y) - ((ncol(X) + 1) * ncol(Y)) # ncol(XZMi)
  logpvals <- pbeta(preScores, v2/2, 1/2, log.p = TRUE)
  scores <- -logpvals / log(10)
  ########################
  rownames(scores) <- colnames(M)
  colnames(scores) <- colnames(Y)

  lastmodel$effects <- bs
  lastmodel$effects.se <- bs.se
  lastmodel$preScores <- preScores
  lastmodel$pvals <- exp(logpvals)
  lastmodel$logpvals <- logpvals
  lastmodel$scores <- scores
  lastmodel$shape1 <- v2/2
  lastmodel$shape2 <- 1/2
  lastmodel$method <- method
  lastmodel$constraints <- res[[8]]
  class(lastmodel) <- c("mmergwas")

  return(lastmodel)
}
