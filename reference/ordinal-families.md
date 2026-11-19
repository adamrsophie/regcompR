# Cumulative link ordinal families

`P(Y <= k) = F((theta_k - mu) / scale)` with `F` the inverse of the
chosen link. Registered as `cumulative_logit` (proportional odds),
`cumulative_probit`, `cumulative_cloglog` (proportional hazards),
`cumulative_loglog` and `cumulative_cauchit`. A regression on `scale`
gives the location-scale model of McCullagh (1980).

## References

McCullagh (1980), Regression models for ordinal data. JRSS B, 42,
109-142.
