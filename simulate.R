# On aws had to install libssl-dev, libcurl4-gnutls-dev, and libxml2-dev and maybe run the following line
# Sys.setenv(PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig")
if (!require("devtools")) install.packages("devtools", dependencies = TRUE)  
if (!require("pacman")) install.packages("pacman")
if (!require("netUtils")) {
  install.packages("devtools")
  devtools::install_github("michaellevy/netUtils")
}
pacman::p_load(statnet, netUtils, tidyverse, broom, stringr, parallel)
set.seed(80112)

cores = detectCores()
nNodes = 100
meanDegs = c(1, 3, 10)
coefs = seq(-.5, .5, len = 4)
decays = 10^(seq(-1, .7, len = 9))

for (simParam in c("twopath", "gwd1.0", "degreepopularity", "twopath", "gwd0.25", "gwd1.0", "gwd3.0")) {
  for (constrained in c(TRUE, FALSE)) {
    for (fixed in c(TRUE, FALSE)) {
      
      start = Sys.time()
      
      # Organized as models/generating_mechanism/constrained-or-density/fixed-or-CEF
      directory = paste0("models/", simParam, "/constrained-", constrained, "/fixed-", fixed)
      if(!dir.exists(directory)) dir.create(directory, recursive = TRUE)
      
      estimates =
        lapply(meanDegs, function(meanDeg) {

          randomGraph = makeNetwork(nNodes, meanDegree = meanDeg, directed = FALSE)

          clust = makeCluster(cores, type = 'FORK')
          coefWithinMeandeg =
            parLapply(clust, coefs, function(dp) {
              
              # Simulate a network with the given parameter value (dp) for this iteration
              simForm = if(str_detect(simParam, "gwd")) {
                simForm = formula(paste0("randomGraph ~ gwdegree(",
                                         str_extract(simParam, "[0-9]+\\.[0-9]+"),
                                         ", fixed = TRUE)"))
              } else {
                simForm = formula(paste0("randomGraph ~ ", simParam))
              }
              simNet =
                simulate.formula(
                  object = simForm
                  , coef = dp
                  , constraints = ~ edges
                  , nsim = 1
                  , control = control.simulate(MCMC.burnin = 1e6, MCMC.interval = 1e6)
                )

              # lapply over theta_s values and estimate an ergm for each
              decaysWithinCoef =
                lapply(decays, function(decay) {

                  # Open connection to write warnings and note parameter values
                  lg = file(file.path(directory, "EstimationLog.txt"), open = "at")
                  sink(lg, type = "message", append = TRUE)
                  message("\nAt ", Sys.time(),
                          "  Sim'ing on ", simParam,
                          ", edges constrained: ", constrained,
                          ", theta fixed: ", fixed,
                          ", mean degree: ", meanDeg,
                          ", coef: ", round(dp, 2),
                          ", decay = ", round(decay, 2))


                  # Estimate
                  m =
                    if(constrained) {
                      try(ergm(simNet ~ gwdegree(decay, fixed = fixed),
                               control = control.ergm(MCMLE.maxit = 50,
                                                      MCMC.samplesize = 5e3,
                                                      seed = 475)
                               , constraints = ~ edges))
                    } else {
                      try(ergm(simNet ~ gwdegree(decay, fixed = fixed) + edges,
                               control = control.ergm(MCMLE.maxit = 50,
                                                      MCMC.samplesize = 5e3,
                                                      seed = 475)))
                    }
                  
                  # Close connection
                  sink(type = "message")
                  close(lg)

                  # If model didn't converge, get rid of the estimates
                  if(class(m) != "ergm" | m$iterations == 50)
                    return(NULL)

                  # Something is causing failures here (I think it's here) with
                  # "Error in svd(X) : a dimension is zero". Rather than 
                  # chasing it down, just return NULL
                  
                  out = 
                    try({
                      tidy(m) %>%
                        mutate(theta_s = decay)
                    })
                  
                  if(!is.data.frame(out)) return(NULL) else return(out)
              
                })

              if(all(sapply(decaysWithinCoef, is.null)))
                return(NULL)

              do.call(rbind, decaysWithinCoef) %>%
                mutate(simCoef = dp)
            })

          stopCluster(clust)

          if(all(sapply(coefWithinMeandeg, is.null)))
            return(NULL)

          do.call(rbind, coefWithinMeandeg) %>%
            mutate(meanDegree = meanDeg)
          

        })

      tab =
        mutate(do.call(rbind, estimates),
               sim_parameter = simParam,
               edges_constrained = constrained,
               fixed_decay = fixed)

      write_csv(tab, file.path(directory, "estimates.csv"))
      
      finish = Sys.time()
      message("At ", Sys.time(), 
              "\nDone with ", simParam, 
              ", edges constrained: ", constrained, 
              ", theta fixed: ", fixed)
      print(finish - start)
      
      
      
    }
  }
}
