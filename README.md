Longitudinal Disease Progression Modeling
================
Adam Lang

-   [Overview](#overview)
-   [Model fitting process](#model-fitting-process)

Overview
========

The trajectory of a given AD biomarker may take decades, however we are
usually limited to much shorter follow up data in ADNI and other AD
focused studies. This can impose a challenge when estimating the full
course of AD disease progression. Using methodology described in Budgeon
et. al \[1\], we can estimate a full-term disease pathology curve from
short term follow up data.

Model fitting process
=====================

In order to estimate the full term disease pathology, we can fit linear
models to each subjects follow up data. We can then estimate the mean
*μ̂*
(
*t*<sub>•</sub>
) and slope
*μ̂*
‘(
*t*<sub>•</sub>
) of each subjects curve where
*t*<sub>•</sub>
is a subjects midpoint. We then fit a polynomial to the points (
*μ̂*
(
*t*<sub>•</sub>
),
*μ̂*
’(
*t*<sub>•</sub>
)). Denoting the fitting values of the polynomial
*μ̂*
’(t) and inputs as
*μ̂*
(t), we can take the reciprocal of the curve, resulting in a new set of
points
$$\\frac{a}{b}$$
