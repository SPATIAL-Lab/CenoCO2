#####
#Preliminaries
#####

#Load libraries
library(rjags)
library(R2jags)
library(openxlsx)

#Read proxy data
df = "data/220602_proxies.xlsx"
d = read.xlsx(df, sheet = "all data product")

#Data subset 
d = d[,c("CO2_ppm", "CO2_uncertainty_pos_ppm", "CO2_uncertainty_neg_ppm",
           "age_Ma", "Age_uncertainty_pos_Ma", "Age_uncertainty_neg_Ma",
           "locality")]
mod = data.frame(280, 5, 5, 0, 0.001, 0.001, "Keeling")
names(mod) = names(d)
d = rbind(d, mod)

#Set up ages vector
ages.bin = 0.5
ages = seq(70, 0, by = 0 - ages.bin) - ages.bin / 2
ages.len = length(ages)

#Parse data - co2 mean and uncertainty
pco2 = log(d$CO2_ppm)
##Max and min in log deviations
pco2.mm = log(d$CO2_ppm + d$CO2_uncertainty_pos_ppm/2) - pco2
pco2.mm = cbind(pco2.mm, pco2 - log(d$CO2_ppm - d$CO2_uncertainty__neg_ppm/2))
##Average 1sd in log units
pco2.sd = apply(pco2.mm, 1, mean)
pco2.pre = 1 / pco2.sd^2

#Parse data - ages and uncertainty
pco2.age = d$age_Ma
pco2.age.sd = apply(cbind(d$Age_uncertainty_pos_Ma, d$Age_uncertainty_neg_Ma), 1, mean)
pco2.loc = d$locality

#Create covariance matrix for ages
pco2.vcov = matrix(nrow = length(pco2), ncol = length(pco2))
diag(pco2.vcov) = pco2.age.sd^2
for(i in 1:(length(pco2)-1)){
  for(j in (i+1):length(pco2)){
    if(pco2.loc[i] == pco2.loc[j]){
      pco2.vcov[i, j] = pco2.vcov[j, i] = 
        0.90 * pco2.age.sd[i] * pco2.age.sd[j] 
    } else{
      pco2.vcov[i, j] = pco2.vcov[j, i] = 0
    }
  }
}
pco2.age.pre = solve(pco2.vcov)

##Data to pass to BUGS model
dat = list(pco2.age = pco2.age, pco2.age.pre = pco2.age.pre, 
           pco2 = pco2, pco2.pre = pco2.pre,
           al = ages.len, ages.bin = ages.bin)

##Parameters to save
parameters = c("pco2_m", "pco2_m.pre", "pco2_m.eps.ac")

##Run it
n.iter = 55000
n.burnin = 5000
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
#p = jags(model.file = "code/parametric_model_walk.R", parameters.to.save = parameters, 
#         data = dat, inits = NULL, n.chains=3, n.iter = n.iter, 
#         n.burnin = n.burnin, n.thin = n.thin)
p = do.call(jags.parallel, list(model.file = "code/model.R", parameters.to.save = parameters, 
                                      data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                      n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

save(p, file = "out/postCeno.rda")
