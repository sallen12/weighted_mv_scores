################################################################################
# vertically re-scaled CRPS

vrcrps_sample <- function(y, dat, w_vr_y = NULL, w_vr_dat = NULL, x0 = NULL, w = NULL){
  # y: vector of realised values
  # dat: vector or matrix (depending on y) of simulation draws from forecast distribution
  # w_vr_y: weight function evaluated at y (vector)
  # w_vr_dat: weight function evaluated at each ensemble member (length equal to number of columns of dat)
  # x0: centring parameter for the vertical re-scaling
  # w: vector or matrix (matching dat) of weights
  
  w_vr_y <- check.uv.weight(y, w_vr_y)
  w_vr_dat <- check.uv.weight(dat, w_vr_dat)
  if (is.null(x0)) {
    x0 <- 0
  }else if (!is.numeric(x0) || length(x0) != 1 ){
    stop("x0 must be numeric of length 1")
  }
  
  input <- list(y = y, dat = dat)
  if (identical(length(y), 1L) && is.vector(dat)) {
    check_sample(input)
    out <- vrcrps_func(y, dat, w_vr_y, w_vr_dat, w)
  } else {
    check_sample2(input)
    out <- sapply(seq_along(y), function(i) vrcrps_func(y[i], dat[i, ], w_vr_y[i], w_vr_dat[i, ], w[i, ]))
  }
  
  if (is.null(w)) {
    s3 <- rowMeans(abs(dat - x0)*w_vr_dat) - abs(y - x0)*w_vr_y
    s4 <- rowMeans(w_vr_dat) - w_vr_y
  } else{
    s3 <- rowSums(w*abs(dat - x0)*w_vr_dat) - abs(y - x0)*w_vr_y
    s4 <- rowSums(w*w_vr_dat) - w_vr_y
  }

  out <- out + s3*s4
  
  return(out)
}

# single vrcrps kernel score
vrcrps_func <- function(y, dat, w_vr_y, w_vr_dat, w = NULL){
  if (is.null(w)) {
    w <- rep(1/length(dat), length(dat))
  } else{
    if (abs(sum(w) - 1) > 1e-5) {
      message("Weights in w don't add to one - they have been re-scaled to one")
    }
    w <- w/sum(w)
  }
  s1 <- sum(vrcrps_kernel(y, dat, w_vr_y, w_vr_dat)*w)
  s2 <- sum(sapply(seq_along(dat), function(m) sum(vrcrps_kernel(dat[m], dat, w_vr_dat[m], w_vr_dat)*w))*w)/2
  return(s1 - s2)
}

# vrcrps_kernel
vrcrps_kernel <- function(x1, x2, w1, w2) abs(x1 - x2)*w1*w2


################################################################################
# vertically re-scaled inverse multiquadric score

vrims_sample_uv <- function(y, dat, w_vr_y = NULL, w_vr_dat = NULL, w = NULL){
  # y: vector of realised values
  # dat: vector or matrix (depending on y) of simulation draws from forecast distribution
  # w_vr_y: weight function evaluated at y (double)
  # w_vr_dat: weight function evaluated at each ensemble member (length equal to number of columns of dat)
  # w: vector or matrix (matching dat) of weights
  
  w_vr_y <- check.uv.weight(y, w_vr_y)
  w_vr_dat <- check.uv.weight(dat, w_vr_dat)
  
  input <- list(y = y, dat = dat)
  if (identical(length(y), 1L) && is.vector(dat)) {
    check_sample(input)
    out <- vrims_func(y, dat, w_vr_y, w_vr_dat, w)
  } else {
    check_sample2(input)
    out <- sapply(seq_along(y), function(i) vrims_func(y[i], dat[i, ], w_vr_y[i], w_vr_dat[i, ], w[i, ]))
  }
  
  return(out)
}

# single vrims kernel score
vrims_func <- function(y, dat, w_vr_y, w_vr_dat, w = NULL){
  if (is.null(w)){
    w <- 1/length(dat)
  } else{
    if (abs(sum(w) - 1) > 1e-5){
      message("Weights in w don't add to one - they have been re-scaled to one")
    }
    w <- w/sum(w)
  }
  s1 <- sum(vrims_kernel_uv(y, dat, w_vr_y, w_vr_dat)*w)
  s2 <- sum(sapply(seq_along(dat), function(m) sum(vrims_kernel_uv(dat[m], dat, w_vr_dat[m], w_vr_dat)*w))*w)/2
  s3 <- vrims_kernel_uv(y, y, w_vr_y, w_vr_y)/2
  out <- s1 - s2 - s3
  return(out)
}

# vrims_kernel
vrims_kernel_uv <- function(x1, x2, w1, w2) (-1/sqrt(1 + abs(x1 - x2)^2))*w1*w2