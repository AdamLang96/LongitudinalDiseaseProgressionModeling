simulated.set.1 <- construct.simulated.dataset(0.4, 15, 10, 20, 1, seq(0,25,.05),
                                               length.subj = 5,  start.sim = 51, eps=0.5, id.start = 1, slope.error = 0, b.error = 0.2)
simulated.set.2 <- construct.simulated.dataset(-0.4, 15, 10, 20, 1, seq(0,25,.05),
                                               length.subj = 5,  start.sim = 51, eps=0.5, id.start = 1, slope.error = .7, b.error = 1)

if(source.plot) {
  plot(simulated.set.1$curve_lines)
  plot(simulated.set.2$curve_lines)
}
n_sample <- seq(.4, 1, by = .05)
n_iter   <- seq(1, 2001, by = 100)
b.error  <- seq(0, .2, by =.1)
ind.lm   <- c(FALSE, TRUE)
parameter_search <- expand.grid(n_sample, n_iter, b.error, ind.lm)
colnames(parameter_search) <- c("n_sample", "n_iter", "b.error", "ind.lm")
init.error.vec <- c()
init.time.vec <- c()
for(i in 1:10) {
  
  simulated.set <- construct.simulated.dataset(0.4, 15, 10, 20, 1, seq(0,25,.05),
                                                 length.subj = 5,  start.sim = 51, 
                                                 eps=0.5, id.start = 1, slope.error = 0, 
                                                 b.error = parameter_search$b.error[i])
  sigmoid.df <- simulated.set$sigmoid_df
  model.data <- simulated.set[["data"]]
  start.time <- Sys.time()
  model.fit  <-  FitDiseaseProgressionCurve(data = model.data, formula.fixed = "Simulated_Response ~ Time_Since_Baseline", 
                                              formula.random = "~1 + Time_Since_Baseline|ID", n_iter = parameter_search$n_iter[i], 
                                              n_sample =parameter_search$n_sample[i], seq.by = 1000, verbose = TRUE, individual.lm = parameter_search$ind.lm[i], 
                                              lme.control =  lmeControl(msTol = 1e-06, msVerbose = FALSE, opt = "optim", msMaxIter = 50))
  end.time   <- Sys.time()
  time.diff <- difftime(start.time, end.time, units = "secs")
  model.fit <- model.fit$Model_Output$Model_Data
  model.fit <- model.fit[complete.cases(model.fit$Domain), ]
  perc.error <- area.between.curves(sigmoid.df, model.fit, perc.error = TRUE)
  init.error.vec <- append(init.error.vec, perc.error)
  init.time.vec <- append(init.time.vec, as.numeric(time.diff))
}


