#' A function to model full term disease pathogenisis 
#'  based on short term longitudinal data
#'  
#'  @param data (data.frame) a data frame containing:
#'   Outcome of interest (numeric)
#'   Column ID with subject id's (factor)
#'   Column Time_Since_Baseline amount of time since baseline observation, can be days, months, years etc where baseline value = 0
#'  @param formula.fixed (character) String with fixed effect formula (should be "outcome ~ Time_Since_Baseline")
#'  @param formula.random (character) String with random effect formula, default reccomended (default is "~1+Time_Since_Baseline|ID")
#'  @param n_iter (numeric) Number of iterations for bootstrapping, (default is 1)
#'  @param n_sample (numeric) Proportion of data preserved for bootstrapping (default is 1)
#'  @param lme.control (list) lmeControl list for mixed effect modeling
#'  @param individual.lm (logical) Logical indicating whether to use multiple linear models fit
#'   to individual subjects in order to calculate Midpoints/Rates instead of Mixed Effects model. Default is reccomended (default is FALSE )
#'  @param seq.by (numeric) # of times to integrate along polynomial reciprocal domain, as the curve is generated pointwise. Upping the value increases fitting time.
#'   (default is 1000)
#'  @param verbose (logical) Prints information during fitting process (default is TRUE)
#'  @returns (list) List with disease progression model data and other related information
#'  
FitDiseaseProgressionCurve <- function(data, formula.fixed, 
                                       formula.random    = "~1+Time_Since_Baseline|ID", 
                                       n_iter            = 1,
                                       n_sample          = 1, 
                                       lme.control       = lmeControl(msTol = 1e-06, msVerbose = FALSE, opt = "optim"), 
                                       individual.lm     = FALSE,
                                       seq.by            = 1000,
                                       verbose           = TRUE) {
  if(!is.data.frame(data)) {
    stop("Data must be of type data.frame")
  }
  
  if("ID" %notin% colnames(data)) {
    stop("Data must include column ID")
  }
  
  if("Time_Since_Baseline" %notin% colnames(data)) {
    stop("Data must include column Time_Since_Baseline")
  }
  
  if(!is.factor(data$ID)) {
    stop("ID must be of class factor")
  }
  
  if(!is.numeric(data$Time_Since_Baseline)) {
    stop("Time_Since_Baseline must be of class numeric")
  }
  
  if(!all(levels(data$ID) %in% unique(data$ID))) {
    stop("Remove missing ID factor level(s) before modeling")
  }
  
  init.bootstrap.dataframe <- list()
  bootstrap.list           <- list()
  build.dictionary         <- list("data"              = data,
                                  "formula.fixed"      = formula.fixed,
                                  "formula.random"     = formula.random,
                                  "n.iterations"       = n_iter,
                                  "n.sample"           = n_sample,
                                  "lmeControl"         = lme.control,
                                  "individual.lm"      = individual.lm,
                                  "verbose"            = verbose)
  
  if(verbose) {
    cat("beginning model fit")
    cat("\n")
  }
  
  # Fits model on entire data before bootstrapping
  EstimateMeanSlopeOutput            <- EstimateMeanSlope(build.dictionary)
  PlotEstimateMeanSlopeOutput        <- PlotMeanSlope(EstimateMeanSlopeOutput)
  FitPolynomialOutput                <- FitPolynomial(EstimateMeanSlopeOutput)
  FindRealRootsOutput                <- FindRealRoots(FitPolynomialOutput)
  CheckRealRootsOutput               <- CheckRealRoots(FindRealRootsOutput                   = FindRealRootsOutput,
                                                       EstimateMeanSlopeOutput               = EstimateMeanSlopeOutput)
  PolynomialCurveOutput              <- DefinePolynomialCurveAndReciprocal(FitPolynomialOutput)
  CalculateBoundsofIntegrationOutput <- CalculateBoundsofIntegration(CheckRealRootsOutput                     = CheckRealRootsOutput,
                                                                     EstimateMeanSlopeOutput                  = EstimateMeanSlopeOutput,
                                                                     DefinePolynomialCurveAndReciprocalOutput = PolynomialCurveOutput,
                                                                     seq.by                                   = seq.by)
  
  init.bootstrap.vector    <- rep(NA, length(CalculateBoundsofIntegrationOutput[["integration_domain"]]))
  for(i in 1:n_iter) {
     # Bootstrapping
    sample.mean.slope             <- sample_n(EstimateMeanSlopeOutput, 
                                              round(n_sample * nrow(EstimateMeanSlopeOutput)))
    sample.polynomial.output      <- FitPolynomial(sample.mean.slope)
    sample.real.roots.output      <- FindRealRoots(sample.polynomial.output)
    sample.check.roots.output     <- CheckRealRoots(sample.real.roots.output,
                                                    sample.mean.slope)
    sample.curve.output           <- DefinePolynomialCurveAndReciprocal(sample.polynomial.output)
    sample.calculate.bounds       <- CalculateBoundsofIntegration(sample.check.roots.output,
                                                                  sample.mean.slope,
                                                                  sample.curve.output,
                                                                  seq.by)

 
    check.bounds                  <- BootStrapCurves(CalculateBoundsofIntegrationOutput,
                                                     sample.calculate.bounds)
    bootstrap.subset              <- IntegratePolynomial(sample.curve.output,
                                                         check.bounds)
    start.index                   <- check.bounds[["start_index"]]
    end.index                     <- check.bounds[["end_index"]]
    bootstrap.vector              <- init.bootstrap.vector
    bootstrap.vector[start.index : end.index] <- bootstrap.subset
    init.bootstrap.dataframe[[i]] <- as.data.frame(bootstrap.vector)
    iter.list                     <- list("iter_polynomial_coefs"   = sample.polynomial.output,
                                          "iter_root_check"         = sample.check.roots.output,
                                          "iter_integration_bounds" = check.bounds)
    
    if(verbose) {
    cat(paste("fit iteration", i, "out of", n_iter, sep = " "))
    cat("\n")
    }
    bootstrap.list[[i]] <- iter.list
   }
  bootstrap.dataframe             <- do.call(bind_cols, 
                                             init.bootstrap.dataframe)
  colnames(bootstrap.dataframe)   <- names(bootstrap.list) <-  paste("iter_", 1 : n_iter, sep = "")
                                           
  CalculateSEOutput               <- CalculateSE(CalculateBoundsofIntegrationOutput,
                                                 bootstrap.dataframe)
  ReorderIfDecreasingOutput       <- ReorderIfDecreasing(CalculateSEOutput,
                                                         CalculateBoundsofIntegrationOutput)
  PlotCurveOutput                 <- PlotCurve(ReorderIfDecreasingOutput)
  
  if(verbose) {
    cat("finished")
    cat("\n")
  }
  
  final.list <- list("Function_Arguments"      = build.dictionary,
                     "Model_Output"            = list("Model_Data" = ReorderIfDecreasingOutput,
                                                      "Model_Plot" = PlotCurveOutput),
                     "Mean_Slope_Output"       = list("Mean_Slope_Data" = EstimateMeanSlopeOutput,
                                                      "Mean_Slope_Plot" = PlotEstimateMeanSlopeOutput),
                     "Polynomial_Coefficients" = FitPolynomialOutput,
                     "Polynomial_Roots"        = CheckRealRootsOutput,
                     "Polynomial/Reciprocal_Functions" = PolynomialCurveOutput,
                     "Integration_Bounds"      = CalculateBoundsofIntegrationOutput,
                     "Bootstrapped_Data"       = bootstrap.dataframe,
                     "Misc_Bootstrap_info"     = bootstrap.list)
  
  return(final.list)
}
