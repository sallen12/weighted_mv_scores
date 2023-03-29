# weighted_mv_scores

This repository contains R code accompanying the paper  

> Allen, S., Ginsbourger, D. and Ziegel, J. (2022). 
> Evaluating forecasts for high-impact events using transformed kernel scores.
> Arxiv pre-print: [arxiv.org/abs/2202.12732](https://arxiv.org/abs/2202.12732).

## Weighted scoring rules

Proper scoring rules provide a well-established framework with which to assess probabilistic forecasts. Weighted scoring rules have been developed that incorporate a weight function into conventional scoring rules, thereby allowing particular outcomes to be emphasised when evaluating forecast performance. Several classes of weighted scoring rules have been proposed, both in the univariate and multivariate case.

This repository contains weighted versions of the **continuous ranked probability score (CRPS)**, **inverse multiquadric score** (IMS), **energy score** (ES), and **variogram score** (VS).

In particular, the following univariate weighted scores can be calculated:
  - Threshold-weighted CRPS (twcrps_sample)
  - Threshold-weighted IMS (twims_sample_uv)
  - Vertically re-scaled CRPS (vrcrps_sample)
  - Vertically re-scaled IMS (vrims_sample_uv)
  
as well as the following multivariate scores:
  - Threshold-weighted ES (twes_sample)
  - Threshold-weighted VS (twvs_sample)
  - Threshold-weighted IMS (twims_sample_mv)
  - Vertically re-scaled ES (vres_sample)
  - Vertically re-scaled VS (vrvs_sample)
  - Vertically re-scaled IMS (vrims_sample_mv)

The weighted scores assume that the forecasts are finite samples from a predictive distribution (i.e. ensemble forecasts). Arbitrary weight and chaining functions can be employed to emphasise particular outcomes.

## scoringRules

Threshold-weighted and outcome-weighted versions of the CRPS, ES, VS, and a Gaussian kernel score are now available within the [scoringRules R package](https://github.com/FK83/scoringRules).
