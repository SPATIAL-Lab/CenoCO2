#####
#Preliminaries
#####

#Load libraries
library(rjags)
library(R2jags)
library(openxlsx)

#Read proxy data
d = read.xlsx("data/211129_proxies.xlsx", sheet = "all data product")

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
pco2.loc = d$Locality

#Creat covariance matrix for ages
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

#ages.s = (70 - d$age_Ma) * 10
#ages.sSD = d$age_uncert*10
#ages = unique(ages)
#ages = sort(ages, decreasing = TRUE)
#age.ind = match(d$age_Ma, ages)

##Data to pass to BUGS model
dat = list(pco2.age = pco2.age, pco2.age.pre = pco2.age.pre, al = ages.len,
           pco2 = pco2, pco2.pre = pco2.pre, ages.bin = ages.bin)

##Parameters to save
parameters = c("pco2_m", "pco2_m.pre", "pco2_m.eps.ac")

##Run it
n.iter = 202000
n.burnin = 2000
n.thin = trunc((n.iter - n.burnin) / 2500)
pt = proc.time()
#p = jags(model.file = "code/parametric_model_ages.R", parameters.to.save = parameters, 
#         data = dat, inits = NULL, n.chains=3, n.iter = n.iter, 
#         n.burnin = n.burnin, n.thin = n.thin)
p = do.call(jags.parallel, list(model.file = "code/parametric_model_ages.R", parameters.to.save = parameters, 
                                      data = dat, inits = NULL, n.chains=6, n.iter = n.iter, 
                                      n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

R2jags::traceplot(p, varname = c("pco2_m"))
R2jags::traceplot(p, varname = c("pco2_m.pre", "pco2_m.eps.ac"))
dev.off()

View(p$BUGSoutput$summary)

save(p, file = "out/post.rda")

sl = p$BUGSoutput$sims.list
su = p$BUGSoutput$summary
sims = nrow(sl$pco2_m)

plot(density(sl$pco2_m.eps.ac), xlim = c(0.5, 1), col = "red")
lines(density(runif(1e6, 0.5, 0.95)))

plot(density(sl$pco2_m.pre), col = "red")
lines(density(rgamma(1e6, shape = 1, rate = 0.01)))

png("out/jpi_simple.png", width = 8, height = 6, units = "in", res = 600)
plot(-10, 0, xlab="Age (Ma)", ylab = expression("pCO"[2]), 
     xlim=c(65,0), ylim=c(100,3000))
for(i in seq(1, sims, by = max(floor(sims / 500),1))){
  lines(ages, exp(sl$pco2_m[i,]), col = rgb(0,0,0, 0.02))
}

arrows(pco2.age, exp(pco2 + 2 * pco2.sd), 
       pco2.age, exp(pco2 - 2 * pco2.sd), 
       length = 0, angle = 90, code = 3, col = "light blue")
arrows(pco2.age - pco2.age.sd, exp(pco2), pco2.age + pco2.age.sd, 
       exp(pco2), length = 0, angle = 90, code = 3, col = "light blue")

lines(ages, exp(su[1:ages.len + 1, 5]), col="red", lw=2)
lines(ages, exp(su[1:ages.len + 1, 3]), col="red", lty=3, lw=2)
lines(ages, exp(su[1:ages.len + 1, 7]), col="red", lty=3, lw=2)

points(pco2.age, exp(pco2), pch=21, bg="white", cex=0.5)

dev.off()

pout = data.frame("Age" = ages, "Mean" = exp(su[1:ages.len + 1, 5]), "ptile_2.5" = exp(su[1:ages.len + 1, 3]), "ptile_97.5" = exp(su[1:ages.len + 1, 7]))
write.csv(pout, "out/pCO2_JPI.csv", row.names = FALSE)

#Make space
mod.p = double()

#Calculate the CDF and find quantile value of modern median in it
for(j in 1:length(ages)){
  cdf = ecdf(exp(sl$pco2_m[, j]))
  mod.p[j] = cdf(408.5)
}

mod.pp = pmin(mod.p, 1 - mod.p) * 2

png("out/modprog.png", width = 8, height = 6, units = "in", res = 600)
plot(ages, mod.pp, type = "l", ylab = "p(408.5)", xlab = "Age (Ma)", ylim=c(0, 0.2), xlim=c(30,0))
abline(0.05, 0, col = "red", lty = 3)
dev.off()
