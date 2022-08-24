################################################################################
# vertically re-scaled energy score

vres_sample <- function(y, dat, w_vr_y = NULL, w_vr_dat = NULL, x0 = NULL, w = NULL){
  # y: realised values (numeric vector of length d)
  # dat: numeric matrix of data (columns are simulation draws from multivariate forecast distribution)
  # w_vr_y: weight function evaluated at y (double)
  # w_vr_dat: weight function evaluated at each ensemble member (length equal to number of columns of dat)
  # x0: centring parameter for the vertical re-scaling
  # w: numeric vector of weights for forecast draws (length equal to number of columns of dat)
  
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  w_vr_y <- check.mv.weight(y, w_vr_y)
  w_vr_dat <- check.mv.weight(dat, w_vr_dat)
  w <- w.helper.multiv(dat, w)
  x0 <- x0.helper.multiv(dat, x0)
  out <- vresC(y, dat, w_vr_y, w_vr_dat, x0, w)
  
  return(out)
}


################################################################################
# vertically re-scaled variogram score

vrvs_sample <- function(y, dat, w_vr_y = NULL, w_vr_dat = NULL, x0 = NULL, w_vs = NULL, w = NULL, p = 0.5){
  # y: realised values (numeric vector of length d)
  # dat: numeric matrix of data (columns are simulation draws from multivariate forecast distribution)
  # w_vr_y: weight function evaluated at y (double)
  # w_vr_dat: weight function evaluated at each ensemble member (length equal to number of columns of dat)
  # x0: centring parameter for the vertical re-scaling
  # w_vs: numeric matrix of weights for dat used in the variogram score. This matrix must be square and symmetric, with all elements being non-negative.
  # w: numeric vector of weights for forecast draws (length equal to number of columns of dat)
  # p: order of variogram score. Standard choices include p = 1 and p = 0.5.
  
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  w_vr_y <- check.mv.weight(y, w_vr_y)
  w_vr_dat <- check.mv.weight(dat, w_vr_dat)
  w <- w.helper.multiv(dat, w)
  w_vs <- w_vs.helper(dat, w_vs)
  x0 <- x0.helper.multiv(dat, x0)
  out <- vrvsC(y, dat, w_vr_y, w_vr_dat, x0, w_vs, w, p)
  
  return(out)
}


################################################################################
# vertically re-scaled inverse multiquadric score

vrims_sample_mv <- function(y, dat, w_vr_y = NULL, w_vr_dat = NULL, w = NULL){
  # y: realised values (numeric vector of length d)
  # dat: numeric matrix of data (columns are simulation draws from multivariate forecast distribution)
  # w_vr_y: weight function evaluated at y (double)
  # w_vr_dat: weight function evaluated at each ensemble member (length equal to number of columns of dat)
  # w: numeric vector of weights for forecast draws (length equal to number of columns of dat)
  
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  w_vr_y <- check.mv.weight(y, w_vr_y)
  w_vr_dat <- check.mv.weight(dat, w_vr_dat)
  w <- w.helper.multiv(dat, w)
  out <- vrimsC(y, dat, w_vr_y, w_vr_dat, w)
  
  return(out)
}