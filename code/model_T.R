model {
  
  #Data model for temperature observations
  
  for(i in 1:length(tData)){
    tData[i] ~ dnorm(t_m[t.ind[i]], 1/4^2)
  }
  
  #Process model
  for(i in 2:al){
    t_m[i] = t_m[i-1] + t_m.eps[i]
    
    t_m.eps[i] ~ dnorm(t_m.eps[i-1] * t_m.eps.ac, t_m.pre) 
  }
  
  t_m.eps[1] ~ dnorm(0, t_m.pre)
  t_m[1] ~ dunif(5, 15)

  #Priors on model parameters  
  t_m.eps.ac ~ dunif(0.01, 0.99)
   
  t_m.pre ~ dgamma(t_m.pre.shp, t_m.pre.rate)
  t_m.pre.shp = 2
  t_m.pre.rate = 2
  
}