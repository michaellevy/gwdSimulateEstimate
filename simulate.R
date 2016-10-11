if(!require("netUtils")) {
  install.packages("devtools")
  devtools::install_github("michaellevy/netUtils")
}
if (!require("pacman")) install.packages("pacman")
pacman::p_load(statnet, tidyverse, netUtils, broom)
set.seed(80112)

nNodes = 100
meanDegs = 1:2 # c(1, 3, 10)
coefs = c(-1, 1)  # seq(-.5, .5, len = 4)
decays = c(.5, 2)  # 10^(seq(-1, .5, len = 5))  # May want better resolution here. 

# Organized as models/generating_mechanism/constrained-or-density/fixed-or-CEF
directory = "models/degPop/constrained/fixed"
dir.create(directory, recursive = TRUE)

estimates = 
  lapply(meanDegs, function(meanDeg) {
    
    randomGraph = makeNetwork(nNodes, meanDegree = meanDeg, directed = FALSE)    
    
    coefWithinMeandeg = 
      lapply(coefs, function(dp) {
        
        # Simulate a network with the given parameter value (dp) for this iteration
        simNet = 
          simulate.formula(
            randomGraph ~ degreepopularity
            , coef = dp
            , constraints = ~ edges
            , nsim = 1
            # , control = control.simulate(MCMC.burnin = 1e4, MCMC.interval = 1e4))
          )
        
        # lapply over theta_s values and estimate an ergm for each
        decaysWithinCoef = 
          lapply(decays, function(decay) {
            
            # Open connection to write warnings and note parameter values
            lg = file(file.path(directory, "EstimationLog.txt"), open = "at")
            sink(lg, type = "message", append = TRUE)
            message("\nAt ", Sys.time(), ", with coef = ", dp, " & decay = ", decay)
            
            # Estimate
            m = try(ergm(simNet ~ gwdegree(decay, fixed = TRUE),
                         control = control.ergm(MCMLE.maxit = 50,
                                                MCMC.samplesize = 5e3,
                                                seed = 475)
                         , constraints = ~ edges
            ))
            
            # Close connection
            sink(type = "message")
            close(lg)
            
            if(class(m) == "try-error")  
              return(NULL)
            
            tidy(m) %>%
              mutate(theta_s = decay)
            
          }) 
        
        if(all(sapply(decaysWithinCoef, is.null)))
          return(NULL)
        
        do.call(rbind, decaysWithinCoef) %>%
          mutate(degPop = dp)
      })
    
    if(all(sapply(coefWithinMeandeg, is.null)))
      return(NULL)
    
    do.call(rbind, coefWithinMeandeg) %>%
      mutate(meanDegree = meanDeg)
    
  })

write_csv(do.call(rbind, estimates), file.path(directory, "estimates.csv"))
