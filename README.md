# weighted_mv_scores

Code to calculate weighted scoring rules in R. The weighted scores are discussed in [this paper](https://arxiv.org/abs/2202.12732).

The weighted scores assume the forecasts are finite samples from a predictive distribution (i.e. ensemble forecasts). Arbitrary weight and chaining functions can be employed to emphasise particular outcomes.

Weighted versions of the **continuous ranked probability score (CRPS)**, **inverse multiquadric score** (IMS), **energy score** (ES), and **variogram score** (VS) are included.

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

Outcome-weighted versions of these scores can be calculated using the existing functionality within the [scoringRules package](https://github.com/FK83/scoringRules). 

