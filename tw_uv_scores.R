################################################################################
# threshold-weighted energy score

twcrps_sample <- function(y, dat, v = identity, w = NULL){
  # y: vector of realised values
  # dat: vector or matrix (depending on y) of simulation draws from forecast distribution
  # v: chaining function to transform the forecasts and observations (function)
  # w: vector or matrix (matching dat) of weights
  
  v_y <- v(y)
  check.chaining(y, v_y)
  v_dat <- v(dat)
  check.chaining(dat, v_dat)
  
  out <- crps_sample(y = v_y, dat = v_dat, w = w)
  
  return(out)
}

################################################################################
# threshold-weighted inverse multiquadric score

twims_sample_uv <- function(y, dat, v = identity, w = NULL){
  # y: vector of realised values
  # dat: vector or matrix (depending on y) of simulation draws from forecast distribution
  # v: chaining function to transform the forecasts and observations (function)
  # w: vector or matrix (matching dat) of weights
  
  v_y <- v(y)
  check.chaining(y, v_y)
  v_dat <- v(dat)
  check.chaining(dat, v_dat)
  
  out <- ims_sample_uv(y = v_y, dat = v_dat, w = w)
  
  return(out)
}

