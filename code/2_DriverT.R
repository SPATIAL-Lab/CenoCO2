# Prep ----

## Load libraries
library(rjags)
library(R2jags)
library(openxlsx)
source("code/Helpers.R")

## Read proxy data
df = "data/Westerhold.xlsx"
d = read.xlsx(df, sheet = "data")

## Data subset 
d = d[,c(1,3)]
names(d) = c("age", "temp")

# 500 kyr for paper ----

## Set up ages vector
ages.bin = 0.5
ages = seq(70, 0, by = 0 - ages.bin) - ages.bin / 2
ages.len = length(ages)

## Age index
d$ai = ceiling((70 - d$age) / ages.bin)

## Data to pass to BUGS model
dat = list("tData" = d$temp, "t.ind" = d$ai, "al" = ages.len)

## Parameters to save
parameters = c("t_m", "t_m.pre", "t_m.eps.ac")

## Run it
n.iter = 12000
n.burnin = 2000
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
p = do.call(jags.parallel, list(model.file = "code/models/model_T.R", parameters.to.save = parameters, 
                                      data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                      n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

save(p, file = "bigout/postTemp.rda")

# 100 kyr for stripes ----

## Set up ages vector
ages.bin = 0.1
ages = agevec(70, ages.bin)
ages.len = length(ages)

## Age index
d$ai = ceiling((70 - d$age) / ages.bin)

## Data to pass to BUGS model
dat = list("tData" = d$temp, "t.ind" = d$ai, "al" = ages.len)

## Parameters to save
parameters = c("t_m", "t_m.pre", "t_m.eps.ac")

## Run it
n.iter = 12000
n.burnin = 2000
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
p = do.call(jags.parallel, list(model.file = "code/models/model_T.R", parameters.to.save = parameters, 
                                data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

save(p, file = "bigout/postTemp100kyr.rda")
