model {
  
  #Data model for pCO2 observations
  
  for(i in 1:length(pco2)){
    pco2[i] ~ dnorm(pco2_m[pco2.age.ind[i]], 1 / pco2.sd[i] ^ 2)
    
  }
  
  #Process model for MgCa_sw timeseries
  
  for(i in 2:length(pco2.ages)){
    pco2_m[i] = pco2_m[i-1] + pco2_m.eps[i]
    
    pco2_m.eps[i] ~ dnorm(pco2_m.eps[i-1] * pco2_m.eps.ac ^ pco2.tau[i], 
                             pco2_m.eps.num / 
                               (1 - pco2_m.eps.ac ^ (2 * pco2.tau[i]))) 
 
    pco2.tau[i] = pco2.ages[i-1] - pco2.ages[i]   
  }
  
  pco2_m.eps.num = (1 - pco2_m.eps.ac ^ 2) * pco2_m.pre
  
  pco2_m.eps[1] ~ dnorm(0, pco2_m.pre)
  pco2_m[1] ~ dunif(4, 8)

  #Priors on MgCa_sw model parameters  
  
  pco2_m.eps.ac ~ dunif(0.05, 0.9999)
  
  pco2_m.pre ~ dgamma(pco2_m.pre.shp, pco2_m.pre.rate)
  pco2_m.pre.shp = 50
  pco2_m.pre.rate = 0.1
  
}