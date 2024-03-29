model {
  
  #Data model for pCO2 observations
  
  for(i in 1:length(pco2)){
    pco2[i] ~ dnorm(pco2.off[i], pco2.pre[i])
    pco2.off[i] ~ dnorm(pco2_m[pco2.aii[i]], pco2.off.pre)
  }

  #Age model
  for(i in 1:length(pco2)){
    pco2.aii[i] = max(min(round((70 - pco2.ai[i]) / ages.bin), al), 1)
  }
  for(i in 1:length(lp)){
    a.off[i] ~ dnorm(0, lp[i])
    for(j in lc[i,1]:lc[i,2]){
      pco2.ai[j] = pco2.age[j] + a.off[i]
    }
  }  
  
  #Process model
  for(i in 2:al){
    pco2_m[i] = pco2_m[i-1] + pco2_m.eps[i]
    
    pco2_m.eps[i] ~ dnorm(pco2_m.eps[i-1] * pco2_m.eps.ac, pco2_m.pre) 
  }
  
  pco2_m.eps[1] ~ dnorm(0, pco2_m.pre)
  pco2_m[1] ~ dunif(6, 8)

  #Priors on model parameters  
  pco2.off.pre ~ dgamma(2, 0.1)
  
  pco2_m.eps.ac ~ dunif(0.01, 0.99)
  
  pco2_m.pre ~ dgamma(pco2_m.pre.shp, pco2_m.pre.rate)
  pco2_m.pre.shp = 1
  pco2_m.pre.rate = 0.01
  
}