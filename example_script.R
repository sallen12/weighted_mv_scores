################################################################################
# example of weighted score function usage

require(scoringRules)

source("vs_sample.R") # load variogram score (with weights)
source("ims_sample.R") # load inverse multiquadric score
source("tw_uv_scores.R") # load univariate threshold-weighted scores
source("vr_uv_scores.R") # load univariate vertically re-scaled scores
source("tw_mv_scores.R") # load multivariate threshold-weighted scores
source("vr_mv_scores.R") # load multivariate vertically re-scaled scores
source("helper_functions.R") # load helper_functions
Rcpp::sourceCpp("cpp_functions.cpp") # load required Cpp functions

n <- 10 # number of observations
d <- 3 # number of dimensions
m <- 20 # number of ensemble members


################################################################################
# univariate scores

# reinsert checkNumeric

y <- rnorm(n)
dat <- matrix(rnorm(n*m), ncol = m)

crps_sample(y, dat)
ims_sample_uv(y, dat)


################################################################################
# multivariate scores

y <- rnorm(d)
dat <- matrix(rnorm(d*m), ncol = m)

es_sample(y, dat)
vs_sample(y, dat)
ims_sample_mv(y, dat)

vs_sample(y, dat)
vs_sample_w(y, dat) # kernel score construction of the variogram score that allows weights
vs_sample_w(y, dat, w = runif(m))


################################################################################
# threshold-weighted univariate scores

t <- 0.5 # threshold in chaining function v(x) = max(x, t) (i.e. w(x) = 1{x > t})

y <- rnorm(n)
dat <- matrix(rnorm(n*m), ncol = m)

crps_sample(y, dat)
twcrps_sample(y, dat)
twcrps_sample(y, dat, v = function(x) pmax(x, -Inf))
twcrps_sample(y, dat, v = function(x) pmax(x, t))

ims_sample_uv(y, dat)
twims_sample_uv(y, dat)
twims_sample_uv(y, dat, v = function(x) pmax(x, -Inf))
twims_sample_uv(y, dat, v = function(x) pmax(x, t))

# check vs numerical integration
dat <- qnorm(seq(1, 1e6)/(1e6 + 1))
twcrps_sample(y[1], dat, v = function(x) pmax(x, t))
integrate(f = function(x) (pnorm(x) - as.numeric(y[1] <= x))^2, lower = t, upper = Inf)


################################################################################
# vertically re-scaled univariate scores

t <- 0.5 # threshold in weight function w(x) = = 1{x > t}
w_vr <- function(x, t) as.numeric(x > t)

y <- rnorm(n)
dat <- matrix(rnorm(n*m), ncol = m)

crps_sample(y, dat)
vrcrps_sample(y, dat)
vrcrps_sample(y, dat, w_vr_y = w_vr(y, -Inf), w_vr_dat = apply(dat, 2, w_vr, -Inf))
vrcrps_sample(y, dat, w_vr_y = w_vr(y, t), w_vr_dat = apply(dat, 2, w_vr, t))
vrcrps_sample(y, dat, w_vr_y = w_vr(y, t), w_vr_dat = apply(dat, 2, w_vr, t), x0 = t)
twcrps_sample(y, dat, v = function(x) pmax(x, t))

ims_sample_uv(y, dat)
vrims_sample_uv(y, dat)
vrims_sample_uv(y, dat, w_vr_y = w_vr(y, -Inf), w_vr_dat = apply(dat, 2, w_vr, -Inf))
vrims_sample_uv(y, dat, w_vr_y = w_vr(y, t), w_vr_dat = apply(dat, 2, w_vr, t))


################################################################################
# outcome-weighted univariate scores

y <- rnorm(n)
dat <- matrix(rnorm(n*m), ncol = m)

t <- 0.5 # threshold in weight function w(x) = = 1{x > t}
w_y <- as.numeric(y > t)
w_dat <- array(as.numeric(dat > t), dim(dat))

crps_sample(y, dat)
crps_sample(y, dat, w = w_dat)*w_y # outcome weighting can be achieved using the inbuilt weight function
crps_sample(y, dat, w = array(0, dim(dat)))*w_y # score is undefined when the weight is zero

ims_sample_uv(y, dat)
ims_sample_uv(y, dat, w = w_dat)*w_y
ims_sample_uv(y, dat, w = array(0, dim(dat)))*w_y


################################################################################
# threshold-weighted multivariate scores

y <- rnorm(d)
dat <- matrix(rnorm(d*m), ncol = m)

t <- rep(0.5, d) # threshold in chaining function v(x) = max(x, t)

es_sample(y, dat)
twes_sample(y, dat)
twes_sample(y, dat, v = function(x) pmax(x, -Inf))
twes_sample(y, dat, v = function(x) pmax(x, t))

vs_sample(y, dat)
twvs_sample(y, dat)
twvs_sample(y, dat, v = function(x) pmax(x, -Inf))
twvs_sample(y, dat, v = function(x) pmax(x, t))

ims_sample_mv(y, dat)
twims_sample_mv(y, dat)
twims_sample_mv(y, dat, v = function(x) pmax(x, -Inf))
twims_sample_mv(y, dat, v = function(x) pmax(x, t))


################################################################################
# vertically re-scaled multivariate scores

y <- rnorm(d)
dat <- matrix(rnorm(d*m), ncol = m)

t <- rep(-0.5, d) # threshold in weight function w(x) = 1{x1 > t1, x2 > t2, ..., xd > td}
w_vr <- function(x, t) all(x > t)

es_sample(y, dat)
vres_sample(y, dat)
vres_sample(y, dat, w_vr_y = 1, w_vr_dat = rep(1, m))
vres_sample(y, dat, w_vr_y = w_vr(y, -Inf), w_vr_dat = apply(dat, 2, w_vr, -Inf))
vres_sample(y, dat, w_vr(y, t), apply(dat, 2, w_vr, t))
vres_sample(y, dat, w_vr(y, t), apply(dat, 2, w_vr, t), x0 = t)
twes_sample(y, dat, v = function(x) x*all(x > t) + t*(1 - all(x > t)))

vs_sample(y, dat)
vrvs_sample(y, dat)
vrvs_sample(y, dat, w_vr_y = 1, w_vr_dat = rep(1, m))
vrvs_sample(y, dat, w_vr_y = w_vr(y, -Inf), w_vr_dat = apply(dat, 2, w_vr, -Inf))
vrvs_sample(y, dat, w_vr(y, t), apply(dat, 2, w_vr, t))
vrvs_sample(y, dat, w_vr(y, t), apply(dat, 2, w_vr, t), x0 = t)
twvs_sample(y, dat, v = function(x) x*all(x > t) + t*(1 - all(x > t)))

ims_sample_mv(y, dat)
vrims_sample_mv(y, dat)
vrims_sample_mv(y, dat, w_vr_y = 1, w_vr_dat = rep(1, m))
vrims_sample_mv(y, dat, w_vr_y = w_vr(y, -Inf), w_vr_dat = apply(dat, 2, w_vr, -Inf))
vrims_sample_mv(y, dat, w_vr(y, t), apply(dat, 2, w_vr, t))


################################################################################
# outcome-weighted multivariate scores

y <- rnorm(d)
dat <- matrix(rnorm(d*m), ncol = m)

t <- rep(-0.5, d) # threshold in weight function w(x) = 1{x1 > t1, x2 > t2, ..., xd > td}
w_y <- as.numeric(all(y > t))
w_dat <- apply(dat > t, 2, prod)

es_sample(y, dat)
es_sample(y, dat, w = w_dat)*w_y 
es_sample(y, dat, w = numeric(ncol(dat)))*w_y

vs_sample_w(y, dat)
vs_sample_w(y, dat, w = w_dat)*w_y
vs_sample_w(y, dat, w = numeric(ncol(dat)))*w_y

ims_sample_mv(y, dat)*w_y
ims_sample_mv(y, dat, w = w_dat)*w_y
ims_sample_mv(y, dat, w = numeric(ncol(dat)))*w_y

