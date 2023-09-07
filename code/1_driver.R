# Preliminaries ----

## Load libraries
library(R2jags)
library(openxlsx)
source("code/8_Helpers.R")

# 500 kyr bins for main text ----

## Set up ages vector
ages.bin = 0.5
ages = agevec(70, ages.bin)
ages.len = length(ages)

## Prep data
d = prepit("230902_proxies.xlsx")

## Parse localities
locs = unique(d$pco2.loc)
lc1 = lc2 = numeric(length(locs))
lp = numeric(length(locs))

## Get locality indices (first-last) and age uncertainties
for(i in 1:length(locs)){
  if(i == 1){lc1[i] = 1}else{lc1[i] = lc2[i-1] + 1}
  lc2[i] = lc1[i] + sum(d$pco2.loc == locs[i]) - 1
  lp[i] = 1 / mean(d$pco2.age.sd[d$pco2.loc == locs[i]])^2
}
lc = cbind(lc1, lc2)

## Data to pass to BUGS model
dat = list(pco2.age = d$pco2.age, lc = lc, lp = lp, 
           pco2 = d$pco2, pco2.pre = d$pco2.pre,
           al = ages.len, ages.bin = ages.bin)

## Parameters to save
parameters = c("pco2_m", "pco2_m.pre", "pco2_m.eps.ac", "pco2.off", "pco2.ai")

## Run it
n.iter = 5e5
n.burnin = 5e4
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
p = do.call(jags.parallel, list(model.file = "code/models/model.R", parameters.to.save = parameters, 
                                      data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                      n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

save(p, file = "bigout/postCenoLERAM.rda")

# 1 Myr bins for SI ----

## Set up ages vector
ages.bin = 1
ages = agevec(70, ages.bin)
ages.len = length(ages)

## Prep data
d = prepit()

## Parse localities
locs = unique(d$pco2.loc)
lc1 = lc2 = numeric(length(locs))
lp = numeric(length(locs))

for(i in 1:length(locs)){
  if(i == 1){lc1[i] = 1}else{lc1[i] = lc2[i-1] + 1}
  lc2[i] = lc1[i] + sum(d$pco2.loc == locs[i]) - 1
  lp[i] = 1 / mean(d$pco2.age.sd[d$pco2.loc == locs[i]])^2
}
lc = cbind(lc1, lc2)

## Data to pass to BUGS model
dat = list(pco2.age = d$pco2.age, lc = lc, lp = lp, 
           pco2 = d$pco2, pco2.pre = d$pco2.pre,
           al = ages.len, ages.bin = ages.bin)

## Parameters to save
parameters = c("pco2_m", "pco2_m.pre", "pco2_m.eps.ac", "pco2.off", "pco2.ai")

## Run it
n.iter = 5e5
n.burnin = 5e4
n.thin = trunc((n.iter - n.burnin) / 2500)

pt = proc.time()
p = do.call(jags.parallel, list(model.file = "code/models/model.R", parameters.to.save = parameters, 
                                data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

save(p, file = "bigout/postCeno1Myr.rda")

# 100 kyr bins for SI ----

## Set up ages vector
ages.bin = 0.1
ages = agevec(70, ages.bin)
ages.len = length(ages)

## Prep data
d = prepit()

## Parse localities
locs = unique(d$pco2.loc)
lc1 = lc2 = numeric(length(locs))
lp = numeric(length(locs))

for(i in 1:length(locs)){
  if(i == 1){lc1[i] = 1}else{lc1[i] = lc2[i-1] + 1}
  lc2[i] = lc1[i] + sum(d$pco2.loc == locs[i]) - 1
  lp[i] = 1 / mean(d$pco2.age.sd[d$pco2.loc == locs[i]])^2
}
lc = cbind(lc1, lc2)

## Data to pass to BUGS model
dat = list(pco2.age = d$pco2.age, lc = lc, lp = lp, 
           pco2 = d$pco2, pco2.pre = d$pco2.pre,
           al = ages.len, ages.bin = ages.bin)

## Parameters to save
parameters = c("pco2_m", "pco2_m.pre", "pco2_m.eps.ac", "pco2.off", "pco2.ai")

## Run it
n.iter = 5e5
n.burnin = 5e4
n.thin = trunc((n.iter - n.burnin) / 2500)

pt = proc.time()
p = do.call(jags.parallel, list(model.file = "code/models/model.R", parameters.to.save = parameters, 
                                data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

save(p, file = "bigout/postCeno100kyr.rda")

# Marine only for SI ----

## Set up ages vector
ages.bin = 0.5
ages = agevec(70, ages.bin)
ages.len = length(ages)

## Prep data
d = prepit("230902_proxies.xlsx")
d = d[d$pco2.prox %in% c("boron isotopes", "phytoplankton"),]

## Parse localities
locs = unique(d$pco2.loc)
lc1 = lc2 = numeric(length(locs))
lp = numeric(length(locs))

## Get locality indices (first-last) and age uncertainties
for(i in 1:length(locs)){
  if(i == 1){lc1[i] = 1}else{lc1[i] = lc2[i-1] + 1}
  lc2[i] = lc1[i] + sum(d$pco2.loc == locs[i]) - 1
  lp[i] = 1 / mean(d$pco2.age.sd[d$pco2.loc == locs[i]])^2
}
lc = cbind(lc1, lc2)

## Data to pass to BUGS model
dat = list(pco2.age = d$pco2.age, lc = lc, lp = lp, 
           pco2 = d$pco2, pco2.pre = d$pco2.pre,
           al = ages.len, ages.bin = ages.bin)

## Parameters to save
parameters = c("pco2_m", "pco2_m.pre", "pco2_m.eps.ac", "pco2.off", "pco2.ai")

## Run it
n.iter = 5e5
n.burnin = 5e4
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
p = do.call(jags.parallel, list(model.file = "code/models/model.R", parameters.to.save = parameters, 
                                data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt
