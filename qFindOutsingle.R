#Sys.setLanguage('en')

library(OUTRIDER)
library(pracma)


#' 
#' Outsingle method to find the latent space dimension q
#' 

qFindOutsingleOpt <- function(ods, type = "ratio"){
  N <- ncol(ods) # number of samples
  
  # control for sequencing depth (2)
  ods <- estimateSizeFactors(ods) 
  controlled_counts <- t(t(counts(ods, normalized=FALSE)) / sizeFactors(ods)) 
  
  # log-transform controlled counts (3)
  row_means <- rowMeans(controlled_counts)
  log_control_counts <- log2((controlled_counts +1) / (row_means +1))
  
  # Compute Z-scores based on mean and standard deviation (5) using the scaling factor c4
  # still unsure where this is coming from and whether we need it
  # formula used in order to normalize data based on sample size
  tryCatch({
    c4 <- sqrt(2 / (N - 1)) * gamma(N / 2) / gamma((N - 1) / 2)
  }, error = function(e) {
    # if OverflowError occurs, use an alternative formula for c4
    c4 <- 1 - 1 / (4 * N) - 7 / (32 * (N ^ 2)) - 19 / (128 * (N ^ 3))
  })
  z_scores <- apply(log_control_counts, 1, calcZScores, c4)
  
  # perform Singular Value Decomposition (SVD) on the matrix of Z-scores 
  # and extract singular values
  sv <- svd(z_scores)$d
  
  # aspect ratio of the count matrix, 0<beta<=1
  beta <- ncol(ods) / nrow(ods) 
  
  # compute the optimal w(beta)
  coef = (optimal_SVHT_coef(beta) / sqrt(median_marcenko_pastur(beta))) # (13)
  print("optimal coefficient: ")
  print(coef)
  
  # compute cutoff (12)
  cutoff = coef * median(sv)
  print("cutoff: ")
  print(cutoff)
  
  # compute and return rank
  greater_than_cutoff = which(sv > cutoff)
  if (length(greater_than_cutoff) > 0){
    k = max(greater_than_cutoff)}
  else {
    k = 0}
  print("Target rank: ")
  return(k)
}

calcZScores <- function(log_control_counts, c4){
  gene_mean <- mean(log_control_counts) # Compute mean of log-controlled counts
  gene_std <- sd(log_control_counts) / c4 # Compute standard deviation with scaling factor c4
  (log_control_counts - gene_mean) / gene_std # Calculate Z-scores
}

# not necessary for our purposes as we can calculate w(beta) numerically
approximation_SVHT_coef <- function(beta){
  0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta + 1.43
}

optimal_SVHT_coef <- function(beta){ 
  # calculate lambda(beta) (14)
  sqrt(2 * (beta + 1) + (8 * beta) / (beta + 1 + sqrt(beta^2 + 14 * beta + 1)))
}


# The following formulas are derived from the robust estimator for the noise 
# parameter gamma in the model Z~ = Z + gamma * E.
# Gamma(Z~) = sigma / sqrt(n * mu)) with mu: median of the Marcenko-Pastur distribution.
# More detailed explanations can be found in Gavish and Donoho 2014.

mar_pas <- function(x, topSpec, botSpec, beta){ 
  # implement Marcenko-Pastur distribution
  if (((topSpec - x) * (x - botSpec)) > 0){
    return(sqrt((topSpec - x) * (x - botSpec)) / (beta * x) / (2 * pi))}
  
  else{
    return(0)}
}

median_marcenko_pastur <- function(beta){ 
  # compute mu: median of Marcenko-Pastur distribution (15)
  botSpec <- (1 - sqrt(beta))^2
  lobnd <- copy(botSpec)
  topSpec <- (1 + sqrt(beta))^2
  hibnd <- copy(topSpec)
  change <- 1
  
  while (change & ((hibnd - lobnd) > 0.001)){ # iterate until convergence
    change <- 0
    x <- seq(lobnd, hibnd, length.out = 10) # range of values for mu
    y <- rep(0, length(x)) # save corresponding results of integral approximation
    for (i in 1:length(x)){
      # approximation of 1 - Integral(15) with Adaptive Gauss-Kronrod Quadrature
      yi <- quadgk(Vectorize(mar_pas), a=x[i], b=topSpec, topSpec=topSpec, botSpec=botSpec, beta=beta)
      # calculation of Integral(15)
      y[i] = 1.0 - yi
    }  
    
    # choose new boundaries for x that yield the closest results to 0.5
    if (any(y < 0.5)){
      lobnd = max(x[y < 0.5])
      change = 1
    }
    
    if (any(y > 0.5)){
      hibnd = min(x[y > 0.5])
      change = 1
    }
  }
  # if hibnd and lobnd are similar enough, return their mean as value for mu
  return((hibnd + lobnd) / 2.0)
}
