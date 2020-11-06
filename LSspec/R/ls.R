# ls.R
# computing the lomb-scargle periodogram
# using the nfft

LSspec <- function( x, t, os = 1, tRange = range(t) ){

  # compute nfft of irregular time series
  sp <- nfft( t, x, os = os, setN = TRUE, tRange = tRange )
  # compute nfft of 1's at frequencies 2w
  spf <- nfft( t, rep(1, length(t)), N = 4*length(t)*os, os = os, tRange = tRange )
  spf <- spf[ rep_len(c(FALSE,TRUE),length(spf)) ]
  
  # compute components of LS periodogram from Press & Rybicki
  Sh <- -Im(sp)    # \sum_i x_i sin(2\pi f t_i)
  Ch <- Re(sp)     # \sum_i x_i cos(2\pi f t_i)
  S2 <- -Im(spf)   # \sum_i 1 sin(2\pi 2f)
  C2 <- Re(spf)    # \sum_i 1 cos(2\pi 2f)
  
  cos2wt <- C2/Mod(spf)
  sin2wt <- S2/Mod(spf)
  coswt <- sign(sin2wt)*sqrt(0.5*(1+cos2wt))
  sinwt <- sqrt(0.5*(1-cos2wt))
  
  
  A1 <- Ch*coswt + Sh*sinwt
  B1 <- Sh*coswt - Ch*sinwt
  A2 <- length(t)/2 + 0.5*C2*cos2wt + 0.5*S2*sin2wt
  B2 <- length(t) - A2
  
  P <- 0.5 * ( A1*A1/A2 + B1*B1/B2 )
  P <- P * 0.5 * diff(tRange) / length(t)
  freq <- seq(from = 0, by = 1/os/diff(tRange), along.with = P )
#browser()
# list( P, Sh, Ch, S2, C2, cos2wt, sin2wt, coswt, sinwt, A1, B1, A2, B2 )  
#data.frame( freq = seq(from = 0, by = 1/os/length(t), along.with = P ), P = P )
data.frame( freq = freq, P = P )
}

dpssApproxFun <- function(t, w, k, tRange = range(t), returnEigenvalues = TRUE){
  # t are the times
  # w is the bandwith
  # k is the number of tapers to return
  # tRange is the total range of times (possibly beyond range(t))
  require(multitaper)
  
  # number of samples
  n <- length(t)
  
  # sample the dpss, then interpolate to needed values
  # The dpss is sampled regularly for accuracy (tridiagonal and whatnot)
  tp <- dpss(n = n, k = k, nw = n * w, returnEigenvalues = returnEigenvalues)
  tpVal <- as.data.frame(tp$v)
  splF <- function(taper, t) splinefun(t, taper)
  tpI <- lapply(tpVal, splF, t = seq(from = tRange[1], to = tRange[2],
                                     length.out = n))
  out <- list(tpI = tpI)
  
  # Return eigenvalues
  if(returnEigenvalues){
    out$eigen <- tp$eigen
  }
  attr(out, "w") <- w
  attr(out, "k") <- k
  attr(out, "n") <- n
  attr(out, "tRange") <- tRange
  out
}

Ftestfun <- function(t, k, x=NULL, w=NULL, tpI = NULL, preprocess = NULL, tpIVal = NULL, freq = NULL, powerK = NULL, tRange = range(t)){
  # Calculate Delta t
  N <- length(t)
  T <- diff(tRange)
  deltat <- T/N
  
  if(is.null(preprocess)){
    # Compute Interpolated DPSS tapers if not provided
    if(is.null(tpI))
      tpI <- dpssApproxFun(t = t, w = w, k = k, tRange = tRange)
    
    tpIVal <- lapply(tpI, do.call, list(t))
    
    # The tapers * data
    pf <- function(taper, t, x, T, dt){
      taper/sum(taper^2) * x * sqrt(T/dt)
    }
    taperedData <- lapply(tpIVal, pf, t = t, x = x, T = T, dt = deltat)
    
    # Compute spectra
    spc <- lapply(taperedData, LSspec, t = t, tRange = tRange)
    
    # Frequency and power
    freq <- spc[[1]]$freq
    powerK <- sapply(spc, getElement, name = "P")
  }
  else{
    if(is.null(tpIVal) | is.null(freq) | is.null(powerK)){
      stop("Values of tapers and Spectrum Estimate must be provided.")
    }
  }
  
  tpIVal <- data.frame(tpIVal)
  tpIVal <- tpIVal*sqrt(deltat)
  Uk0 <- apply(tpIVal, 2, sum)    ## Percival and Walden H0
  
  ## Even tapers are symmetric and hence their summation goes to 0
  if(k >= 2) {
    Uk0[seq(2,k,2)] <- 0
  }
  
  # H0 sum of squares
  Uk0SumSq <- as.numeric(t(Uk0)%*%Uk0)
  
  # Mean estimate of the amplitude of the periodic signal
  MuEstimate <- (powerK %*% Uk0) /Uk0SumSq
  # Variance explained by the presence of the periodic signal in the data
  SigVar <-  (Mod(MuEstimate)^2)*Uk0SumSq
  
  # Variance of residual coloured noise
  Uk0 <- as.matrix(Uk0)
  NoiseVar <- apply(Mod(powerK - (MuEstimate %*% t(Uk0)))^2,
                    1, sum)
  # k-1 degrees of freedom
  F_ <- (k-1)*SigVar/NoiseVar
  out <- data.frame(F_ = F_, freq = freq)
  out
}

LSspecMT <- function(t, x, os = 1, tRange = range(t), w, k, tpI = NULL, subtract.mean = TRUE, niterations = NULL, Ftest = TRUE){
  # Number of samples
  N <- length(t)
  # Total time window
  T <- diff(tRange)
  
  if(subtract.mean) x <- x - mean(x)
  
  # Compute Interpolated DPSS tapers if not provided
  if(is.null(tpI)){
    dpssIN <- dpssApproxFun(t, w, k, tRange = tRange, returnEigenvalues = TRUE)
    tpI <- dpssIN$tpI
    lambdaK <- dpssIN$eigen    # Eigenvalues are required only if adaptive weighting is performed
  }
  
  tpIVal <- lapply(tpI, do.call, list(t))
  
  # The tapers * data
  pf <- function(taper, t, x, T, dt){
    taper/sum(taper^2) * x * sqrt(T/dt) 
  }
  taperedData <- lapply(tpIVal, pf, t = t, x = x, T = T, dt = T/N)
  
  # Compute spectra - you can use lsp in place of LSspec for comparison
  spc <- lapply(taperedData, LSspec, t = t, os = os, tRange = tRange)
  # Frequency and power
  freq <- spc[[1]]$freq    # For lsp, this becomes spc$V1$scanned
  powerK <- sapply(spc, getElement, name = "P") # For lsp, this becomes name = "power"
  
  # Averaging to get multitaper statistic
  Sf <- rowMeans(powerK)
  
  # Adaptive weighting
  if(!is.null(niterations)){
    for(i in 1:niterations){
      # Weighting using local "signal" and broad-band "noise"
      dK <- Sf/(Sf%*%t(lambdaK) + (1 - lambdaK)*var(x))
      # Spectrum Estimate
      Sf <- rowMeans(dK^2*lambdaK*powerK)/sum(dK^2*lambdaK)
    }
  }
  out <- data.frame(freq = freq, P = Sf)
  
  if(Ftest){
    ftestRes <- Ftestfun(t, k, preprocess = TRUE, tpIVal = tpIVal, freq = freq, powerK = powerK)
    attr(out, "Ftest") <- ftestRes
  }
  
  attr( out, "w" ) <- w
  attr( out, "k" ) <- k
  attr( out, "n" ) <- N
  attr( out, "tRange" ) <- tRange
  attr( out, "os" ) <- os
  attr( out, "subtract.mean" ) <- subtract.mean
  out
}

PspecMT <- function( t, x, os = 1, tRange = range(t), w, k, tpI = NULL, subtract.mean = TRUE, niterations = NULL, Ftest = TRUE){
  # Number of samples
  N <- length(t)
  # Total time window
  T <- diff(tRange)
  
  if(subtract.mean) x <- x - mean(x)
  
  # Compute Interpolated DPSS tapers if not provided
  if(is.null(tpI)){
    dpssIN <- dpssApproxFun(t, w, k, tRange = tRange, returnEigenvalues = TRUE)
    tpI <- dpssIN$tpI
    lambdaK <- dpssIN$eigen    # Eigenvalues are required only if adaptive weighting is performed
  }
  
  tpIVal <- lapply(tpI, do.call, list(t))
  
  # The tapers * data
  pf <- function(taper, t, x, T, dt){
    taper/sum(taper^2) * x * sqrt(T/dt) 
  }
  taperedData <- lapply(tpIVal, pf, t = t, x = x, T = T, dt = T/N)
  
  # Redefine nfft function
  nfft.2 <- function(x, t, os, setN, tRange) nfft(t, x, os = os, setN = TRUE, tRange = tRange)
  # Compute spectra
  spc <- lapply(taperedData, nfft.2, t = t, os = os, tRange = tRange)
  # Frequency and Power
  freq <- seq(from = 0, by = 1/os/diff(tRange), along.with = spc[[1]])
  powerK <- sapply(spc, function(x) abs(x)^2)
  
  # Averaging to get multitaper statistic
  # Note that we multiply by deltat * oversample
  Sf <- rowMeans(powerK) * 0.5 * T / N / N
  
  if(!is.null(niterations)){
    for(i in 1:niterations){
      # Weighting using local "signal" and broad-band "noise"
      dK <- Sf/(Sf%*%t(lambdaK) + (1 - lambdaK)*var(x))
      # Spectrum Estimate
      Sf <- rowMeans(dK^2*lambdaK*powerK)/sum(dK^2*lambdaK)
    }
  }
  out <- data.frame(freq = freq, P = Sf)
  
  if(Ftest){
    ftestRes <- Ftestfun(t, k, preprocess = TRUE, tpIVal = tpIVal, freq = freq, powerK = powerK)
    attr(out, "Ftest") <- ftestRes
  }
  
  attr( out, "w" ) <- w
  attr( out, "k" ) <- k
  attr( out, "n" ) <- N
  attr( out, "tRange" ) <- tRange
  attr( out, "os" ) <- os
  attr( out, "subtract.mean" ) <- subtract.mean
  out
}

## Slow Non-uniform Discrete Fourier Transform (NUDFT)
nudft <- function(z, t, inverse=FALSE) {
  N <- length(t)
  if(N == 0) return(z)
  T <- diff(range(t))
  ff <- (if(inverse) 1 else -1) * 1i * 2*pi/T * (t - t[1])
  vapply(1:N, function(h) sum(z * exp(ff*(h-1))), complex(1))
}

NUspecMT <- function(t, x, os = 1, tRange = range(t), w, k, tpI = NULL, subtract.mean = TRUE, niterations = NULL, Ftest = TRUE){
  # Number of samples
  N <- length(t)
  # Total time window
  T <- diff(tRange)
  
  if(subtract.mean) x <- x - mean(x)
  
  # Compute Interpolated DPSS tapers if not provided
  if( is.null(tpI) ){
    dpssIN <- dpssApproxFun(t, w, k, tRange = tRange, returnEigenvalues = TRUE)
    tpI <- dpssIN$tpI
    lambdaK <- dpssIN$eigen
  }
  
  tpIVal <- lapply(tpI, do.call, list(t))
  
  # The tapers * data
  pf <- function(taper, t, x, T, dt){
    taper/sum(taper^2) * x * sqrt(T/dt) 
  }
  taperedData <- lapply(tpIVal, pf, t = t, x = x, T = T, dt = T/N)
  
  # Compute spectra
  spc <- lapply(taperedData, nudft, t = t)
  # Frequency and Power
  freq <- seq(from = 0, by = 1/os/diff(tRange), along.with = spc[[1]])
  powerK <- sapply(spc, function(x) abs(x)^2)
  
  # Averaging to get multitaper statistic
  # Note that we multiply by deltat * oversample
  Sf <- rowMeans(powerK) * 0.5 * T / N / N
  
  if(!is.null(niterations)){
    for(i in 1:niterations){
      # Weighting using local "signal" and broad-band "noise"
      dK <- Sf/(Sf%*%t(lambdaK) + (1 - lambdaK)*var(x))
      # Spectrum Estimate
      Sf <- rowMeans(dK^2*lambdaK*powerK)/sum(dK^2*lambdaK)
    }
  }
  out <- data.frame(freq = freq, P = Sf)
  
  if(Ftest){
    ftestRes <- Ftestfun(t, k, preprocess = TRUE, tpIVal = tpIVal, freq = freq, powerK = powerK)
    attr(out, "Ftest") <- ftestRes
  }
  
  attr( out, "w" ) <- w
  attr( out, "k" ) <- k
  attr( out, "n" ) <- N
  attr( out, "tRange" ) <- tRange
  attr( out, "os" ) <- os
  attr( out, "subtract.mean" ) <- subtract.mean
  out
}


