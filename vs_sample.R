################################################################################
# variogram score (with weights)

vs_sample_w <- function(y, dat, w_vs = NULL, w = NULL, p = 0.5) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  if (!is.numeric(p) || length(p) != 1 ){
    stop("Order 'p' must be numeric of length 1")
  } else if (p < 0) {
    stop("Order 'p' must be positive")
  }
  w <- w.helper.multiv(dat, w)
  w_vs <- w_vs.helper(dat, w_vs)
  vs <- vsC_w(y, dat, w_vs, w, p)
  return(vs)
}
