model {
  
  #Data model for pCO2 observations
  
  for(i in 1:length(pco2)){
    pco2[i] ~ dnorm(pco2_m[pco2.aii[i]], pco2.pre.g[i])
    pco2.pre.g[i] = ifelse(group[i] == alkmeth || group[i] == 3, 
                           pco2.pre[i], pco2.pre[i] / 1e6)
  }
  
  #Age model
  for(i in 1:length(pco2)){
    pco2.aii[i] = max(min(round((70 - pco2.ai[i]) * 1 / ages.bin), al), 1)
  }
  pco2.ai ~ dmnorm(pco2.age, pco2.age.pre)
  
  #Dataset selector
  alkmeth ~ dcat(alkpp[])
  alkpp[1] = 0.5
  alkpp[2] = 0.5
  
  #Process model
  for(i in 2:al){
    pco2_m[i] = pco2_m[i-1] + pco2_m.eps[i]
    
    pco2_m.eps[i] ~ dnorm(pco2_m.eps[i-1] * pco2_m.eps.ac, pco2_m.pre) 
  }
  
  pco2_m.eps[1] ~ dnorm(0, pco2_m.pre)
  pco2_m[1] ~ dunif(6, 8)

  #Priors on model parameters  
  
  pco2_m.eps.ac ~ dunif(0.5, 0.9999)
  
  pco2_m.pre ~ dgamma(pco2_m.pre.shp, pco2_m.pre.rate)
  pco2_m.pre.shp = 1
  pco2_m.pre.rate = 0.01
  
}