# Some of these functions have been significantly edited by me, Wouter van der Bijl.
# Truncated_TestScore, or SHet, now support parallel computation using {future}

library("MASS")
library("Matrix")
library("compiler")

solar_cormat <- function(sol) {
  s <- sol$solar$files$model$polygenic.out
  s <- s[grep('RhoG', s)]
  s <- str_split_i(s, ' is ', 2)[c(TRUE, FALSE)]
  r <- readr::parse_double(s)

  n <- length(r)
  k <- (1 + sqrt(1 + 8 * n)) / 2

  m <- matrix(nrow = k, ncol = k)
  m[lower.tri(m)] <- r
  m <- t(m)
  m[lower.tri(m)] <- r
  diag(m) <- 1

  return(m)
}

combine_assoc <- function(files) {
  stopifnot(length(files) > 1)
  d <- data.table::fread(files[1], data.table = FALSE)
  # detect if the standard errors are in column "se", or "beta.se"
  if ('beta.se' %in% names(d)) {
    se_name <- 'beta.se'
    d$Z1 <- d$beta / d[[se_name]]
    d <- dplyr::select(d, chr:alt, Z1, beta1 = beta)
  } else {
    se_name <- 'se'
    d$Z1 <- d$beta / d[[se_name]]
    d <- dplyr::select(d, chr:af, Z1, beta1 = beta)
  }

  Zs <- map_dfc(2:length(files), \(i) {
    D <- data.table::fread(files[i], data.table = FALSE)
    D$Z <- D$beta / D[[se_name]]
    D <- D[c('Z', 'beta')]
    names(D) <- paste0(c('Z', 'beta'), i)
    return(D)
  })

  bind_cols(d, Zs) %>%
    dplyr::relocate(starts_with('Z'), starts_with('beta'), .after = last_col()) %>%
    as_tibble()
}

Non_Trucated_TestScore <- function(X, SampleSize, CorrMatrix) {
	Wi = matrix(SampleSize, nrow = 1);
	sumW = sqrt(sum(Wi^2));
	W = Wi / sumW;

	Sigma = ginv(CorrMatrix);
	XX = apply(X, 1, function(x) {
		x1 <- matrix(x, ncol = length(x), nrow = 1);
		Tmat = W %*% Sigma %*% t(x1);
		Tmat = (Tmat*Tmat) / (W %*% Sigma %*% t(W));
		return(Tmat[1,1]);
		}
	);
	return(XX);
}
SHom <- cmpfun(Non_Trucated_TestScore);

Trucated_TestScore <- function(
    X, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05,
    isAllpossible = TRUE
) {
	N <- nrow(X)

	Wi <- matrix(SampleSize, nrow = 1)
	sumW <- sqrt(sum(Wi^2))
	W <- Wi / sumW

	calc_TTT <- function(x) {
	  if (isTRUE(all.equal(var(x), 0))) return(0)
	  TTT = -1

	  if (isAllpossible == T ) {
	    cutoff <- sort(unique(abs(x)))	  ## it will filter out any of them.
	  } else {
	    cutoff = seq(startCutoff, endCutoff, CutoffStep)
	  }

	  for (threshold in cutoff) {
	    x1 <- x
	    index <- which(abs(x1) < threshold)

	    if (length(index) == N) break;

	    A <- CorrMatrix

	    W1 <- W
	    if (length(index) != 0) {
	      x1 <- x1[-index]
	      A  <- A[-index, -index]   ## update the matrix
	      W1 <- W[-index]
	    }

	    if (correct == 1)	{
	      W1 <- W1 * sign(x1)
	    }

	    A <- ginv(A)
	    x1 <- matrix(x1, nrow = 1)
	    W1 <- matrix(W1, nrow = 1)
	    Tmat <- W1 %*% A %*% t(x1)
	    Tmat <- (Tmat * Tmat) / (W1 %*% A %*% t(W1))

	    if (TTT < Tmat[1, 1]) TTT <- Tmat[1, 1]
	  }
	  return(TTT)
	}
	res <- future.apply::future_apply(X, 1, calc_TTT, future.seed = TRUE)

	return(res)
}
SHet <- cmpfun(Trucated_TestScore);

EstimateGamma <- function (N = 1E6, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = T) {

	Wi = matrix(SampleSize, nrow = 1);
	sumW = sqrt(sum(Wi^2));
	W = Wi / sumW;

	Permutation = mvrnorm(n = N, mu = c(rep(0, length(SampleSize))), Sigma = CorrMatrix, tol = 1e-8, empirical = F);

	Stat =  Trucated_TestScore(X = Permutation, SampleSize = SampleSize, CorrMatrix = CorrMatrix,
					correct = correct, startCutoff = startCutoff, endCutoff = endCutoff,
					CutoffStep = CutoffStep, isAllpossible = isAllpossible);
    	a = min(Stat)*3/4
	ex3 = mean(Stat*Stat*Stat)
	V =	var(Stat);

	for (i in 1:100){
		E = mean(Stat)-a;
		k = E^2/V
		theta = V/E
		a = (-3*k*(k+1)*theta**2+sqrt(9*k**2*(k+1)**2*theta**4-12*k*theta*(k*(k+1)*(k+2)*theta**3-ex3)))/6/k/theta
	}

	para = c(k,theta,a);
	return(para);
}

EmpDist <- function (N = 1E6, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = T) {

	Wi = matrix(SampleSize, nrow = 1);
	sumW = sqrt(sum(Wi^2));
	W = Wi / sumW;

	Permutation = mvrnorm(n = N, mu = c(rep(0, length(SampleSize))), Sigma = CorrMatrix, tol = 1e-8, empirical = F);

	Stat =  Trucated_TestScore(X = Permutation, SampleSize = SampleSize, CorrMatrix = CorrMatrix, correct = correct, startCutoff = startCutoff, endCutoff = endCutoff, CutoffStep = CutoffStep, isAllpossible = isAllpossible);


	return(Stat);
}




