model {
  
  #Data model for pCO2 observations
  
  for(i in 1:length(pco2)){
    pco2[i] ~ dnorm(pco2_m[round(pco2.ai[i])], pco2.pre[i])
    
    pco2.ai[i] ~ dnorm(pco2.age[i], pco2.age.pre[i]) T (1, 700.99)
  }
  
  #Process model
  
  for(i in 2:al){
    pco2_m[i] = pco2_m[i-1] + pco2_m.eps[i]
    
    pco2_m.eps[i] ~ dnorm(pco2_m.eps[i-1] * pco2_m.eps.ac, pco2_m.pre) 
  }
  
  pco2_m.eps[1] ~ dnorm(0, pco2_m.pre)
  pco2_m[1] ~ dunif(4, 8)

  #Priors on model parameters  
  
  pco2_m.eps.ac ~ dunif(0.05, 0.9999)
  
  pco2_m.pre ~ dgamma(pco2_m.pre.shp, pco2_m.pre.rate)
  pco2_m.pre.shp = 5
  pco2_m.pre.rate = 0.1
  
}