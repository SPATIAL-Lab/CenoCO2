#####
#Preliminaries
#####

#Load libraries
library(rjags)
library(R2jags)
library(openxlsx)

#Read proxy data
df = "data/Westerhold.xlsx"
d = read.xlsx(df, sheet = "data")

#Data subset 
d = d[,c(1,3)]
names(d) = c("age", "temp")

#Set up ages vector
ages.bin = 0.5
ages = seq(70, 0, by = 0 - ages.bin) - ages.bin / 2
ages.len = length(ages)

#Age index
d$ai = ceiling((70 - d$age) / ages.bin)

##Data to pass to BUGS model
dat = list("tData" = d$temp, "t.ind" = d$ai, "al" = ages.len)

##Parameters to save
parameters = c("t_m", "t_m.pre", "t_m.eps.ac")

##Run it
n.iter = 12000
n.burnin = 2000
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
#p = jags(model.file = "code/model_T.R", parameters.to.save = parameters, 
#         data = dat, inits = NULL, n.chains=3, n.iter = n.iter, 
#         n.burnin = n.burnin, n.thin = n.thin)
p = do.call(jags.parallel, list(model.file = "code/model_T.R", parameters.to.save = parameters, 
                                      data = dat, inits = NULL, n.chains = 4, n.iter = n.iter, 
                                      n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

save(p, file = "out/postTemp.rda")
