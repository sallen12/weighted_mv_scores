################################################################################
# threshold-weighted energy score

twes_sample <- function(y, dat, v = identity, w = NULL){
  # y: realised values (numeric vector of length d)
  # dat: numeric matrix of data (columns are simulation draws from multivariate forecast distribution)
  # v: chaining function to transform the forecasts and observations (function)
  # w: numeric vector of weights for forecast draws (length equal to number of columns of dat)
  
  v_y <- v(y)
  check.chaining(y, v_y)
  v_dat <- apply(dat, 2, v)
  check.chaining(dat, v_dat)
  
  out <- es_sample(y = v_y, dat = v_dat, w = w)
  
  return(out)
}


################################################################################
# threshold-weighted variogram score

twvs_sample <- function(y, dat, v = identity, w_vs = NULL, p = 0.5){
  # y: realised values (numeric vector of length d)
  # dat: numeric matrix of data (columns are simulation draws from multivariate forecast distribution)
  # v: chaining function to transform the forecasts and observations (function)
  # w_vs: numeric matrix of weights for dat used in the variogram score. This matrix must be square and symmetric, with all elements being non-negative.
  # p: order of variogram score. Standard choices include p = 1 and p = 0.5.
  
  v_y <- v(y)
  check.chaining(y, v_y)
  v_dat <- apply(dat, 2, v)
  check.chaining(dat, v_dat)
  
  out <- vs_sample(y = v_y, dat = v_dat, w_vs = w_vs, p = p)
  
  return(out)
}


################################################################################
# threshold-weighted inverse multiquadric score

twims_sample_mv <- function(y, dat, v = identity, w = NULL){
  # y: realised values (numeric vector of length d)
  # dat: numeric matrix of data (columns are simulation draws from multivariate forecast distribution)
  # v: chaining function to transform the forecasts and observations (function)
  # w: numeric vector of weights for forecast draws (length equal to number of columns of dat)
  
  v_y <- v(y)
  check.chaining(y, v_y)
  v_dat <- apply(dat, 2, v)
  check.chaining(dat, v_dat)
  
  out <- ims_sample_mv(y = v_y, dat = v_dat, w = w)
  
  return(out)
}
