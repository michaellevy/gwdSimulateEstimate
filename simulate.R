if(!require("netUtils")) {
  install.packages("devtools")
  devtools::install_github("michaellevy/netUtils")
}
if (!require("pacman")) install.packages("pacman")
pacman::p_load(statnet, dplyr, ggplot2, netUtils, stargazer)
set.seed(80112)

n = makeNetwork(100, meanDegree = 1, directed = FALSE)    
nSim = simulate.formula(
  n ~ gwdegree(.1, TRUE)
  , coef = -3
  , constraints = ~ edges
  , nsim = 1
  , control = control.simulate(MCMC.burnin = 1e4, MCMC.interval = 1e4))
decays = c(.1, .5, 1.5, 3)
# 
# ergms = lapply(decays, function(decay) 
#   try(ergm(nSim ~ edges + gwdegree(decay, fixed = TRUE),
#            control = control.ergm(MCMLE.maxit = 50,
#                                   MCMC.samplesize = 5e3)))
# )

ergmsC = lapply(decays, function(decay) 
  structure(
    try(ergm(nSim ~ gwdegree(decay, fixed = TRUE),
             control = control.ergm(MCMLE.maxit = 50,
                                    MCMC.samplesize = 5e3),
             constraints = ~ edges)),
    "decay" = decay)
)

# stargazer(ergms[sapply(ergms, class) == "ergm"], type = "text")
ergmsC = ergmsC[sapply(ergmsC, class) == "ergm"]
stargazer(ergmsC, type = "text")

sims = lapply(ergmsC, function(model)
  simulate.ergm(model, contstraints = ~ edges, nsim = 500))

degDists = 
  do.call(rbind,
          lapply(seq_along(sims), function(i) {
            deg = getDegreeDist(sims[[i]])
            mutate(deg,
                   decay = attr(ergmsC[[i]], 'decay'))
          }) 
  )

# Calculate variance of degree distribution
vars = lapply(seq_along(sims), function(i) {
  c(decay = attr(ergmsC[[i]], "decay"), 
    variance = mean(sapply(sims[[i]], function(n) var(degree(n, gmode = 'graph')))))
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

vars
var(degree(nSim, gmode = "graph"))  # yikes

png("results/varyDecayAndEst_gofs.png", width = 800, height = 1000)
par(mfrow = c(4, 3))
for(i in seq_along(ergmsC)) {
  decay = attr(ergmsC[[i]], "decay")
  plot(gof(ergmsC[[i]]))
  mtext(decay)
}
dev.off()

png("results/sampleGraphs.png", width = 800, height = 1000)
par(mfrow = c(3, 2))
plot(nSim, main = "basis")
for(i in 1:length(sims))
  plot(sims[[i]][[1]], main = paste0("Decay = ",
                                     attr(ergmsC[[i]], "decay"),
                                     "\nVariance = ",
                                     round(vars$variance[i], 2)))
dev.off()
