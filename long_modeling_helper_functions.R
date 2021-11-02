#--------------------------------------------------------------------------------------------- Libraries

library(purrr) 
library(plyr) 
library(dplyr)
library(stringr)
library(nlme)
library(matrixStats)
library(ggplot2)
library(RConics)

#--------------------------------------------------------------------------------------------- Functions

#' Negation of %in% function
#' 
`%notin%` <- Negate(`%in%`)

#--------------------------------------------------------------------------------------------- 

#' Creates data.frame of Midpoints vs Rates for
#'  individual subjects curves
#'  
#' @param BuildDictionary (list) A list of arguments from  \emph{FitDiseaseProgressionCurve} function
#' @return (data.frame) A data frame with columns ID, Midpoints, Rates
#' 

EstimateMeanSlope <- function(BuildDictionary, test.average, test.new) {
  method.logical  <- BuildDictionary[["individual.lm"]]
  data            <- BuildDictionary[["data"]]
  formula.fixed   <- BuildDictionary[["formula.fixed"]]
  formula.random  <- BuildDictionary[["formula.random"]]
  model.control   <- BuildDictionary[["lmeControl"]]
  rate.vec        <- c()
  midpoint.vec.av <- midpoint.vec.rate <- c()

  if(method.logical) {
    split.data      <- split(data, data$ID)
    ids             <- names(split.data)
    for(i in 1:length(split.data)) {
      subj          <- split.data[[i]]
      subj.lm       <- lm(formula = as.formula(formula.fixed), data = subj)
      rate          <- coef(subj.lm)[["Time_Since_Baseline"]]
      pred.model    <- predict(subj.lm)
      midpoint      <- (max(pred.model) + min(pred.model)) / 2
      rate.vec      <- append(rate.vec, rate)
      midpoint.vec  <- append(midpoint.vec, midpoint)
    }
    mean.slope.data <- data.frame("ID"        = ids,
                                  "Rates"     = rate.vec,
                                  "Midpoints" = midpoint.vec)
    
  } else {
    model      <- lme(as.formula(formula.fixed), 
                      random  = as.formula(formula.random), 
                      data    = data,
                      control = model.control)
    re         <- nlme :: ranef(model)
    fe         <- nlme :: fixef(model)
    splitdata  <- split(data, data$ID)
    ids        <- names(splitdata)
    for(i in ids) {
      subj       <- splitdata[[i]]
      int        <- re["(Intercept)"][i,]
      int        <- int + fe[["(Intercept)"]]
      rate       <- re["Time_Since_Baseline"][i, ]
      rate       <- rate + fe[["Time_Since_Baseline"]]
      rate.vec   <- append(rate.vec, rate)
      pred.val   <- stats::predict(object = model)
      keeps      <- which(names(pred.val) == i)
      pred.val   <- pred.val[keeps]
      
      midpoint.av   <- (max(pred.val) + min(pred.val)) / 2
      midpoint.vec.av  <- append(midpoint.vec.av, midpoint.av)
      
      midpoint.rate <- int + (rate * (.5 * max(subj$Time_Since_Baseline)))
      midpoint.vec.rate  <- append(midpoint.vec.rate, midpoint.rate)
      
      }
    mean.slope.data <- data.frame("ID"          = ids,
                                  "Rates"       = rate.vec,
                                  "Midpoints" = midpoint.vec.rate)
    
  }
  
  return(mean.slope.data)
}
#---------------------------------------------------------------------------------------------

#' Fits 3rd degree poly over Midpoints/Rates data
#' 
#' @param EstimateMeanSlopeOutput (data.frame) Output from \emph{EstimateMeanSlope} function
#' @return (list) Coefficients of fitted 3rd degree polynomial
#' 
FitPolynomial <- function(EstimateMeanSlopeOutput) {
  data <- EstimateMeanSlopeOutput
  third.degree.poly <- lm(Rates ~ poly(Midpoints, 3, raw = TRUE), 
                                        data = data)
  poly.coefs <- coef(third.degree.poly)
  poly.coefs <- list("coef.x3"  = poly.coefs[["poly(Midpoints, 3, raw = TRUE)3"]], 
                     "coef.x2"  = poly.coefs[["poly(Midpoints, 3, raw = TRUE)2"]], 
                     "coef.x1"  = poly.coefs[["poly(Midpoints, 3, raw = TRUE)1"]], 
                     "int"      = poly.coefs[["(Intercept)"]])
  
  return(poly.coefs)
}

#---------------------------------------------------------------------------------------------
#' Calculates the roots (real) of the fitted polynomial
#' 
#' @param FitPolynomialOutput (list) Output from \emph{FitPolynomial} function
#' @return (vector) Vector of real roots of fitted polynomial
#' 
FindRealRoots <- function(FitPolynomialOutput) {
  poly.coefs  <- FitPolynomialOutput
  poly.coefs  <- c(poly.coefs[["coef.x3"]], poly.coefs[["coef.x2"]], 
                   poly.coefs[["coef.x1"]], poly.coefs[["int"]])
  roots        <- RConics::cubic(poly.coefs)
  roots        <- as.complex(roots)
  roots.keep   <- which(Im(roots) == 0)
  roots        <- roots[roots.keep]
  roots        <- Re(roots)
  return(roots)
}

#---------------------------------------------------------------------------------------------
#' Checks to see whether the roots satisfy criterion for integration
#'  
#' @param FindRealRootsOutput (vector) Output from \emph{FindRealRoots} function
#' @param EstimateMeanSlopeOutput (data.frame) Output from \emph{EstimateMeanSlope} function
#' @return (list) List of all real roots with logical value indicating
#'  whether they satisfy integration conditions
#'  
CheckRealRoots <- function(FindRealRootsOutput, EstimateMeanSlopeOutput) {
  all.roots     <- c()
  roots.satisfy <- list()
  roots         <- FindRealRootsOutput
  midpoints     <- EstimateMeanSlopeOutput$Midpoints
  min.midpoint  <- min(midpoints)
  max.midpoint  <- max(midpoints)
  if(length(roots) == 0) {
    roots.satisfy[[1]] <- list("root" = NA,
                          "satisfies_condition" = TRUE)
    all.roots <- TRUE
  } else {
    for(i in 1:length(roots)) {
      if(roots[i] >= min.midpoint & roots[i] <= max.midpoint) {
        logi               <- FALSE
        all.roots          <- append(all.roots, logi)
        roots.satisfy[[i]] <- list("root" = roots[i],
                                   "meets_condition" = logi)
      } else {
        logi               <- TRUE
        all.roots          <- append(all.roots, logi)
        roots.satisfy[[i]] <- list("root" = roots[i],
                                   "meets_condition" = logi)
      }
    }
  }
  all.roots              <- all(all.roots)
  roots.satisfy[["all_roots_satisfy?"]] <- all.roots
  return(roots.satisfy)
}

#---------------------------------------------------------------------------------------------
#' Defines functions for 3rd degree polynomial and reciprocal of
#'  3rd degree polynomial from fitted polynomial coefficients
#'  
#'  @param FitPolynomialOutput (list) Output from \emph{FitPolynomial} function
#'  @return (list) List of functions for checking polynomial fit visually 
#'   and integrating polynomial reciprocal
#'   
DefinePolynomialCurveAndReciprocal <- function(FitPolynomialOutput) {
  poly.coefs      <- FitPolynomialOutput
  PolynomialCurve <- function(x) {
    (poly.coefs[["coef.x3"]] * (x^3)) +
    (poly.coefs[["coef.x2"]] * (x^2)) +
    (poly.coefs[["coef.x1"]] * (x))   +
     poly.coefs[["int"]]
  }
  ReciprocalCurve <- function(x) {
     1 / ((poly.coefs[["coef.x3"]] * (x^3)) +
          (poly.coefs[["coef.x2"]] * (x^2)) +
          (poly.coefs[["coef.x1"]] * (x))   +
           poly.coefs[["int"]])
  }
    
 return(list("Polynomial_Function" = PolynomialCurve,
             "Reciprocal_Function" = ReciprocalCurve))
}

#---------------------------------------------------------------------------------------------
#' Calculates the window of integration based on min and max midpoints,
#'  roots of the fitted polynomial, and the direction of the curve.
#'  
#'  @param CheckRealRootsOutput (list) Output from \emph{CheckRealRoots} function
#'  @param EstimateMeanSlopeOutput (data.frame) Output from \emph{EstimateMeanSlope} function
#'  @param DefinePolynomialCurveAndReciprocalOutput (list) Output from \emph{DefinePolynomialCurveAndReciprocal} function
#'  @param seq.by (numeric) # of times to integrate along curve (default is 1000)
#'  @return (list) a list of values containing the integration domain,
#'   whether the domain was subsetted due to roots of the polynomial, and
#'   the direction of the curve (positive/negative)
#'   

CalculateBoundsofIntegration <- function(CheckRealRootsOutput, EstimateMeanSlopeOutput, DefinePolynomialCurveAndReciprocalOutput, seq.by) {
  polynomial.curve <- DefinePolynomialCurveAndReciprocalOutput[["Polynomial_Function"]]
  min.midpoint.row <- which.min(EstimateMeanSlopeOutput$Midpoints)
  max.midpoint.row <- which.max(EstimateMeanSlopeOutput$Midpoints)
  min.midpoint     <- EstimateMeanSlopeOutput["Midpoints"][min.midpoint.row, ]
  max.midpoint     <- EstimateMeanSlopeOutput["Midpoints"][max.midpoint.row, ]
  min.midpoint.y   <- polynomial.curve(min.midpoint)
  max.midpoint.y   <- polynomial.curve(max.midpoint)
  number.positive.rates <- length(which(EstimateMeanSlopeOutput$Rates > 0))
  number.negative.rates <- length(which(EstimateMeanSlopeOutput$Rates < 0))
  a <- min.midpoint
  b <- max.midpoint
  head_subset <- FALSE
  tail_subset <- FALSE

  if(number.positive.rates < number.negative.rates) {
    direction <- "decreasing"
  } else {
    direction <- "increasing"
  }
  
  if(!CheckRealRootsOutput[["all_roots_satisfy?"]]) {
    roots <- CheckRealRootsOutput
    roots[["all_roots_satisfy?"]] <- NULL
    which.fails    <- unlist(map(roots, pluck, 2))
    which.fails    <- !which.fails
    roots.fails    <- unlist(map(roots[which.fails], pluck, 1))
  if(direction == "increasing") {    # positive rates
    if(length(roots.fails) == 1) {
      if(min.midpoint.y < 0 & max.midpoint.y > 0) {
         a           <- roots.fails
         head_subset <- TRUE
      } else if(min.midpoint.y > 0 & max.midpoint.y < 0) {
         b           <- roots.fails
         tail_subset <- TRUE
       }
    } else if(length(roots.fails == 2)) {
        a <- min(roots.fails)
        b <- max(roots.fails)
        head_subset <- TRUE
        tail_subset <- TRUE
    } else if(length(roots.fails) == 3) {
        return("3_FAILS")
      }
  } else {                         # negative rates 
      if(length(roots.fails) == 1) {
        if(min.midpoint.y < 0 & max.midpoint.y > 0) {
          b           <- roots.fails
          tail_subset <- TRUE
        } else if(min.midpoint.y > 0 & max.midpoint.y < 0) {
          a           <- roots.fails
          head_subset <- TRUE
        }
      } else if(length(roots.fails == 2)) {
        a <- min(roots.fails)
        b <- max(roots.fails)
        head_subset <- TRUE
        tail_subset <- TRUE
      } else if(length(roots.fails) == 3) {
        return("3_FAILS")
      }
    }
  }
  seq.by = (b - a) / seq.by
  integration.domain <- seq(a, b, by = seq.by)
  if(head_subset & tail_subset) {
  integration.domain <- integration.domain[2 : (length(integration.domain) - 1)]
  } else if (head_subset & !tail_subset) {
    integration.domain <- integration.domain[2 : length(integration.domain)]
  } else if(!head_subset & tail_subset) {
    integration.domain <- integration.domain[1 : (length(integration.domain) - 1)]
  } else {
    integration.domain <- integration.domain
  }
  return(list("integration_start"    = a,
              "integration_end"      = b,
              "integration_domain"   = integration.domain,
              "head_subset"          = head_subset,
              "tail_subset"          = tail_subset,
              "direction"            = direction))
  }
  





#---------------------------------------------------------------------------------------------
#' Calculates integral of polynomial along polynomial domain
#' 
#' @param DefinePolynomialCurveAndReciprocalOutput (list) Output of \emph{DefinePolynomialCurveAndReciprocal} function
#' @param CalculateBoundsofIntegrationOutput (list) Output of \emph{CalculateBoundsofIntegration} function
#' @return (vector) Vector of integrated values 
#' 
IntegratePolynomial <- function(DefinePolynomialCurveAndReciprocalOutput, CalculateBoundsofIntegrationOutput) {
  integration.function      <- DefinePolynomialCurveAndReciprocalOutput[["Reciprocal_Function"]]
  integration.domain        <- CalculateBoundsofIntegrationOutput[["integration_domain"]]
  integrated.values.vector  <- c(0)
  integration.domain <- as.vector(integration.domain)
  for(k in 2:length(integration.domain)) {
    integration.val <- try(integrate(integration.function,
                                 lower = min(integration.domain),
                                 upper = integration.domain[k])$val, silent = TRUE)
    integration.val <- suppressWarnings(as.numeric(integration.val))
    integrated.values.vector <- suppressWarnings(append(integrated.values.vector, integration.val))
    if(is.na(integration.val)) {
      warning("NA value generated during integration")
    }
  }

  return(integrated.values.vector)
}

#---------------------------------------------------------------------------------------------
#' Calculates mean and standard error at each point along the domain of the
#'  curve from boostrap iterations
#'  
#' @param CalculateBoundsofIntegrationOutput (list) Output of \emph{CalculateBoundsofIntegration} function
#' @param BootStrappedDF (data.frame) Data frame of bootstrapped values
#' @return (data.frame) Disease progression model data
#' 

CalculateSE <- function(CalculateBoundsofIntegrationOutput, BootStrappedDF) {
  BootStrappedDF <- as.matrix(BootStrappedDF)
  response       <- CalculateBoundsofIntegrationOutput[["integration_domain"]]
  Conf.Low       <- rep(NA, nrow(BootStrappedDF))
  Conf.Hi        <- rep(NA, nrow(BootStrappedDF))
  domain         <- rep(NA, nrow(BootStrappedDF))
  for(i in 1 : nrow(BootStrappedDF)) {
    variance      <- quantile(BootStrappedDF[i, ], 
                              probs = c(0.05, .95), na.rm = TRUE)
    Conf.Low[i] <- variance[[1]]
    Conf.Hi[i]  <- variance[[2]]
    row.mean    <- mean(BootStrappedDF[i,])
    domain[i]   <- row.mean
    
  }
  disease.progression.data <- cbind(response, domain, Conf.Low, Conf.Hi)
  colnames(disease.progression.data) <- c("Response", "Domain", "CI_Low", "CI_Hi")
  return(disease.progression.data)
}

#---------------------------------------------------------------------------------------------
#' Determines whether to reorder the disease progression model data depending on if the curve
#'  is decreasing
#'  
#' @param CalculateSEOutput (matrix) Output from \emph{CalculateSE} function
#' @param CalculateBoundsofIntegrationOutput (list) Output from \emph{CalculateBoundsofIntegration} function
#' @return (data.frame) Disease progression data reordered if necessary
#' 
ReorderIfDecreasing <- function(CalculateSEOutput, CalculateBoundsofIntegrationOutput) {
  direction <- CalculateBoundsofIntegrationOutput[["direction"]]
  disease.progression.data <- as.data.frame(CalculateSEOutput)
  if(direction == "decreasing") {
    disease.progression.data["Domain"]              <- disease.progression.data["Domain"] * -1
    disease.progression.data["CI_Low"]              <- disease.progression.data["CI_Low"] * -1
    disease.progression.data["CI_Hi"]               <- disease.progression.data["CI_Hi"]  * -1
    disease.progression.data["Response"]            <- disease.progression.data["Response"][nrow(disease.progression.data):1,]
  }
 return(disease.progression.data)  
}

#---------------------------------------------------------------------------------------------
#' Recalculate bounds of integration for each bootstrap iteration
#'  to ensure integration step does not fail
#'  
#' @param IntegrationBoundsFull (list) CalculateBoundsofIntegrationOutput list
#' @param IntegrationBoundsSample (list) output from \emph{CalculateBoundsofIntegration} on bootstrapped data
#' @return (list) list with integration domain, start/end points of integration, index of start/end points
#' 
BootStrapCurves <- function(IntegrationBoundsFull, IntegrationBoundsSample) {
  integration.domain.full         <- IntegrationBoundsFull[["integration_domain"]]
  integration.domain.full.start   <- IntegrationBoundsFull[["integration_start"]]
  integration.domain.full.end     <- IntegrationBoundsFull[["integration_end"]]
  integration.domain.sample.start <- IntegrationBoundsSample[["integration_start"]]
  integration.domain.sample.end   <- IntegrationBoundsSample[["integration_end"]]
  start.values                    <- max(c(integration.domain.full.start, integration.domain.sample.start))
  end.values                      <- min(c(integration.domain.full.end, integration.domain.sample.end))
  integration.domain              <- integration.domain.full[integration.domain.full >= start.values & integration.domain.full <= end.values]
  integration.start               <- integration.domain[1]
  integration.end                 <- integration.domain[length(integration.domain)]
  integration.start.index         <- which(integration.domain == integration.start)
  integration.end.index           <- which(integration.domain == integration.end)
  bootstrap.out.list              <- list("integration_domain" = integration.domain,
                                          "integration_start"  = integration.start,
                                          "integration_end"    = integration.end,
                                          "start_index"        = integration.start.index,
                                          "end_index"          = integration.end.index)
                                          
                          
  return(bootstrap.out.list)
}

BootStrapCurves.test <- function(IntegrationBoundsFull, IntegrationBoundsSample) {
  integration.domain.full         <- IntegrationBoundsFull[["integration_domain"]]
  integration.domain.full.start   <- IntegrationBoundsFull[["integration_start"]]
  integration.domain.full.end     <- IntegrationBoundsFull[["integration_end"]]
  integration.domain.sample.start <- IntegrationBoundsSample[["integration_start"]]
  integration.domain.sample.end   <- IntegrationBoundsSample[["integration_end"]]
  start.values                    <- max(c(integration.domain.full.start, integration.domain.sample.start))
  end.values                      <- min(c(integration.domain.full.end, integration.domain.sample.end))
  return(list("start.values" = start.values,
         "end.values" = end.values,
         "integration.domain" = integration.domain.full))
  integration.domain              <- integration.domain.full[integration.domain.full >= start.values & integration.domain.full <= end.values]
  integration.start               <- integration.domain[1]
  integration.end                 <- integration.domain[length(integration.domain)]
  integration.start.index         <- which(integration.domain == integration.start)
  integration.end.index           <- which(integration.domain == integration.end)
  bootstrap.out.list              <- list("integration_domain" = integration.domain,
                                          "integration_start"  = integration.start,
                                          "integration_end"    = integration.end,
                                          "start_index"        = integration.start.index,
                                          "end_index"          = integration.end.index)
  
  
  return(bootstrap.out.list)
}


#--------------------------------------------------------------------------------------------- Generate plots
#' Plots disease progression curve
#' 
#' @param data (data.frame) Disease progression data frame
#' @return (ggplot) Plot of disease progression curve
#' 
PlotCurve <- function(data) {
  na.remove       <- which(is.na(data$Domain))
  if(length(na.remove) >= 1) {
   warning("NA values were removed generating plot")
   data            <- data[-na.remove,]
    }
  init.plot       <- ggplot(data = data, aes(x = Domain,  y = Response)) + geom_point()
  plot.curve      <- init.plot  + geom_line(aes(x=CI_Low, y = Response), linetype  = "dashed")
  plot.curve      <- plot.curve + geom_line(aes(x=CI_Hi,  y = Response), linetype  = "dashed")
  return(plot.curve)
 }

#---------------------------------------------------------------------------------------------
#' Plots Midpoints vs Rates
#' 
#' @param data (data.frame) Midpoints vs Rates data frame (EstimateMeanSlopeOutput)
#' @return (ggplot) Plot of Midpoints vs Rates
#' 
PlotMeanSlope <- function(data) {
  plot        <- ggplot(data = data, aes(x = Midpoints, y = Rates)) + geom_point()
  return(plot)
}

#---------------------------------------------------------------------------------------------


