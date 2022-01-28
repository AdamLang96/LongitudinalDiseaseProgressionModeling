Longitudinal Disease Progression Modeling
================
Adam Lang

-   [Overview](#overview)
-   [Installation](#installation)
-   [Modeling Process](#modeling-process)
-   [FitDiseaseProgressionCurve](#fitdiseaseprogressioncurve)
    -   [Arguments](#arguments)
    -   [Returns](#returns)
-   [Example](#example)
    -   [Data](#data)
    -   [Model Fitting](#model-fitting)
-   [Plots](#plots)
-   [Comparison](#comparison)
-   [Notes](#notes)
-   [Citations](#citations)
    -   [Packages](#packages)

Overview
========

The trajectory of a given AD biomarker may take decades, however we are
usually limited to much shorter follow up data in ADNI and other AD
focused studies. This can impose a challenge when estimating the full
course of AD disease progression. Using methodology described in Budgeon
et. al [\[1\]](#1), we can estimate a full-term disease pathology curve
from short term follow up data.

<br>

Installation
============

Both ***FitDiseaseProgressionCurve.R*** and
***FitDiseaseProgressionCurveHelpers.R*** must be downloaded. The
following packages must be installed:

``` r
install.packages("purrr") 
install.packages("plyr") 
install.packages("dplyr")
install.packages("stringr")
install.packages("nlme")
install.packages("matrixStats")
install.packages("ggplot2")
install.packages("RConics")
install.packages("zoo")
```

<br>

Modeling Process
================

The modeling process to is outlined in detail in [Budgeon et.
al](https://pubmed.ncbi.nlm.nih.gov/28444781/)

<br>

FitDiseaseProgressionCurve
==========================

Arguments
---------

**FitDiseaseProgressionCurve** requires the following arguments:

***data*** (data.frame) a data frame containing:

   Outcome of Interest <code>(numeric)</code>  
   Column <code>ID</code> subject ids <code>(factor)</code>  
   Column <code>Time\_Since\_Baseline</code> <code>(numeric)</code> time
since baseline obs

***formula.fixed*** <code>(character)</code> fixed effect formula
<code>(should be “Outcome \~ Time\_Since\_Baseline”) </code>

***formula.random*** <code>(character)</code> random effect formula
<code>(default is “\~1+Time\_Since\_Baseline\|ID”)</code>

***n\_iter*** <code>(numeric)</code> number of iterations for
bootstrapping, <code>(default is 1)</code>

***n\_sample*** <code>(numeric)</code> proportion of data preserved for
bootstrapping <code>(default is 1)</code>

***lme.control*** <code>(list)</code> lmeControl list for mixed effect
modeling <code>(optional)</code>

***individual.lm*** <code>(logical)</code> use linear models for each
subject instead of one mixed effects model <code>(default is
FALSE)</code>

***seq.by*** <code>(numeric)</code> \# of times to integrate along
polynomial domain <code>(default is 1000)</code>

Returns
-------

**FitDiseaseProgressionCurve** Returns a list with the following
objects:

***Function\_Arguments*** <code>(list)</code> A list of the arguments
specified in the model fitting

***Model\_Output*** <code>(list)</code> A list with:

   <code>Model\_Data (data.frame)</code> data of estimated population
model  
   <code>Model\_Plot (ggplot)</code> plot of estimated population model

***Mean\_Slope\_Output*** <code>(list)</code> A list with:

   <code>Mean\_Slope\_Data (data.frame)</code> data of estimated means
vs. slopes  
   <code>Mean\_Slope\_Plot (ggplot)</code> plot of estimated means
vs. slopes

***IntegrationBounds*** <code>(list)</code> A list with:

   <code>direction (character)</code> whether the model is increasing or
decreasing (positive or negative)  
   <code>roots.frame (data.frame)</code> a data.frame with the
integration start and end points

Example
=======

Here we will generate simulated data to show model fitting process

Data
----

``` r
#Functions for simulating data in SimulationFunctions.R
set.seed(123)

#generate data based on sigmoid curve
example.data.list <- construct.simulated.dataset(0.4, 13, .01, .98, 20, seq(0,25,.1),
                                            length.subj = 3,  start.sim = 51, eps = 0.5, id.start = 1)
```

Example of data frame for model

``` r
data <- example.data.list$data
knitr::kable(data[1:14,])
```

| Time\_Since\_Baseline | Simulated\_Response | ID  |
|----------------------:|--------------------:|:----|
|                   0.0 |           0.0099764 | 1   |
|                   0.5 |           0.0099843 | 1   |
|                   1.0 |           0.0099923 | 1   |
|                   1.5 |           0.0100002 | 1   |
|                   2.0 |           0.0100082 | 1   |
|                   2.5 |           0.0100161 | 1   |
|                   3.0 |           0.0100240 | 1   |
|                   0.0 |           0.0099753 | 2   |
|                   0.5 |           0.0099876 | 2   |
|                   1.0 |           0.0099999 | 2   |
|                   1.5 |           0.0100123 | 2   |
|                   2.0 |           0.0100246 | 2   |
|                   2.5 |           0.0100369 | 2   |
|                   3.0 |           0.0100492 | 2   |

Model Fitting
-------------

``` r
set.seed(123)
model.fit  <-  FitDiseaseProgressionCurve(data           = data, 
                                          formula.fixed  = "Simulated_Response ~ Time_Since_Baseline", 
                                          formula.random = "~1 + Time_Since_Baseline|ID", 
                                          n_iter         = 200, 
                                          n_sample       = .75, 
                                          seq.by         = 1000, 
                                          verbose        = FALSE, 
                                          individual.lm  = FALSE) 
```

Estimated population model

``` r
model.fit$Model_Output$Model_Plot
```

[![modelfit.png](https://i.postimg.cc/yN2n25cG/modelfit.png)](https://postimg.cc/wtchmwW5)
<br>

Plots
=====

<br>

Population sigmoidal curve that we will base our data generating process
on

[![init-curve.png](https://i.postimg.cc/tJ3r351M/init-curve.png)](https://postimg.cc/K3YrFtDP)
<br>

Longitudinal subject data generated based on our sigmoidal curve with
noise introduced

[![curve-lines.png](https://i.postimg.cc/j5G9Q2ML/curve-lines.png)](https://postimg.cc/1nrJymhQ)

<br>

Subjects realigned so that baseline is T=0

[![time-bline.png](https://i.postimg.cc/vHy2QpkF/time-bline.png)](https://postimg.cc/xqs35F9t)

<br>

Comparison
==========

Here the red curve is the population sigmoid curve we generated and the
black curve is the population curve estimated in the model.

[![curvewithlines.png](https://i.postimg.cc/c1RT5FHD/curvewithlines.png)](https://postimg.cc/ThPrKqnm)

Notes
=====

The length of the domain is data dependant, as the origin of the curve
begins at the smallest observed mean in the data. However, the main
utility of the disease progression curve is to estimate the time between
disease states, so the origin of the model and the total time span of
the model is unimportant.

Citations
=========

<a id="1">\[1\]</a> Budgeon, Charley et al. “Constructing longitudinal
disease progression curves using sparse, short-term individual data with
an application to Alzheimer’s disease.” Statistics in medicine 36 17
(2017): 2720-2734 .

Packages
--------

\[2\] Lionel Henry and Hadley Wickham (2019). purrr: Functional
Programming Tools. R package version 0.3.2.
<a href="https://CRAN.R-project.org/package=purrr" class="uri">https://CRAN.R-project.org/package=purrr</a>

\[3\] Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data
Analysis. Journal of Statistical Software, 40(1), 1-29. URL
<a href="http://www.jstatsoft.org/v40/i01/" class="uri">http://www.jstatsoft.org/v40/i01/</a>.

\[4\] Hadley Wickham, Romain François, Lionel Henry and Kirill Müller
(2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7.
<a href="https://CRAN.R-project.org/package=dplyr" class="uri">https://CRAN.R-project.org/package=dplyr</a>

\[5\] Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2021). nlme:
Linear and Nonlinear Mixed Effects Models. R package version 3.1-153,
&lt;URL:
<a href="https://CRAN.R-project.org/package=nlme" class="uri">https://CRAN.R-project.org/package=nlme</a>&gt;.

\[6\] Hadley Wickham (2019). stringr: Simple, Consistent Wrappers for
Common String Operations. R package version 1.4.0.
<a href="https://CRAN.R-project.org/package=stringr" class="uri">https://CRAN.R-project.org/package=stringr</a>

\[7\] Henrik Bengtsson (2020). matrixStats: Functions that Apply to Rows
and Columns of Matrices (and to Vectors). R package version 0.57.0.
<a href="https://CRAN.R-project.org/package=matrixStats" class="uri">https://CRAN.R-project.org/package=matrixStats</a>

\[8\] H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
Springer-Verlag New York, 2016.

\[9\] Emanuel Huber (2014). RConics: Computations on Conics. R package
version 1.0.
<a href="https://CRAN.R-project.org/package=RConics" class="uri">https://CRAN.R-project.org/package=RConics</a>

\[10\] Achim Zeileis and Gabor Grothendieck (2005). zoo: S3
Infrastructure for Regular and Irregular Time Series. Journal of
Statistical Software, 14(6), 1-27.
<a href="doi:10.18637/jss.v014.i06" class="uri">doi:10.18637/jss.v014.i06</a>
