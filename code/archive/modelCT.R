model {
  
  #Data model for pCO2 observations
  for(i in 1:ns){
    pco2[i] ~ dnorm(pco2_m[a.r[i]], pco2.pre[i])
  }
  
  #Process model
  for(i in 2:nall){
    pco2_m[i] = pco2_m[i-1] + pco2_m.eps[i]
    
    pco2_m.eps[i] ~ dnorm(pco2_m.eps[i-1] * pco2_m.eps.ac ^ (2 * tau[i]), 
                          pco2_m.eps.num / (1.00001 - pco2_m.eps.ac ^ (2 * tau[i]))) 
    
    tau[i] = (ages[i-1] - ages[i])
  }
  
  #Age model
  ##rank and sort
  a.r = rank(-ai)
  ages = -sort(-ai)
  ##combine with fixed time points
  ai = c(pco2.ai, ts)
  ##pick data ages
  for(i in 1:length(lp)){
    a.off[i] ~ dnorm(0, lp[i])
    for(j in lc[i,1]:lc[i,2]){
      pco2.ai[j] = pco2.age[j] + a.off[i]
    }
  }
  #pco2.ai ~ dmnorm(pco2.age, pco2.age.pre)
  
  #Initial condition
  pco2_m.eps[1] ~ dnorm(0, pco2_m.pre)
  pco2_m[1] ~ dunif(6, 8)

  #TS model constant
  pco2_m.eps.num = (1 - pco2_m.eps.ac ^ 2) * pco2_m.pre

  #Priors
  ##time series autocorrelation
  pco2_m.eps.ac ~ dunif(0.01, 0.99)
  ##time series error
  pco2_m.pre ~ dgamma(pco2_m.pre.shp, pco2_m.pre.rate)
  pco2_m.pre.shp = 1
  pco2_m.pre.rate = 0.01
  
}