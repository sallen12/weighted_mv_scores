################################################################################
# inverse multiquadric score

# multivariate 
ims_sample_mv <- function(y, dat, w = NULL) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  w <- w.helper.multiv(dat, w)
  out <- imsC(y, dat, w)
  return(out)
}

# univariate
ims_sample_uv <- function(y, dat, w = NULL) {
  input <- list(y = y, dat = dat)
  if (!is.null(w)) input$w <- w
  if (identical(length(y), 1L) && is.vector(dat)) {
    check_sample(input)
    out <- ims_func(y, dat, w)
  } else {
    check_sample2(input)
    out <- sapply(seq_along(y), function(i) ims_func(y[i], dat[i, ], w[i, ]))
  }
  return(out)
}

# single ims kernel score
ims_func <- function(y, dat, w = NULL){
  if (is.null(w)){
    w <- 1/length(dat)
  } else{
    if (abs(sum(w) - 1) > 1e-5){
      message("Weights in w don't add to one - they have been re-scaled to one")
    }
    w <- w/sum(w)
  }
  s1 <- sum(ims_kernel_uv(y, dat)*w)
  s2 <- sum(sapply(seq_along(dat), function(m) sum(ims_kernel_uv(dat[m], dat)*w))*w)/2
  s3 <- ims_kernel_uv(y, y)/2
  out <- s1 - s2 - s3
  return(out)
}

# ims kernel
ims_kernel_uv <- function(x1, x2) -1/sqrt(1 + abs(x1 - x2)^2)
