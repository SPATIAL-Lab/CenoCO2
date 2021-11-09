#####
#Preliminaries
#####

#Load libraries
library(rjags)
library(R2jags)
library(openxlsx)

#Read proxy data
d = read.xlsx("data/211109_proxies.xlsx", sheet = "all data product")

#Set up ages vector
ages = seq(70, 0, by = -0.1)
ages.len = length(ages)

#Parse data - co2 mean and uncertainty
pco2 = log(d$CO2_ppm)
##Max and min in log deviations
pco2.mm = log(d$CO2_ppm + d$CO2_uncertainty_pos_ppm/2) - pco2
pco2.mm = cbind(pco2.mm, pco2 - log(d$CO2_ppm - d$CO2_uncertainty__neg_ppm/2))
##Average 1sd in log units
pco2.sd = apply(pco2.mm, 1, mean)

#Parse data - ages and uncertainty
pco2.age = d$age_Ma
pco2.age.sd = apply(d[,c("Age_uncertainty_pos_Ma", "Age_uncertainty_neg_Ma")], 1, mean)
pco2.age.loc = d$Locality

#ages.s = (70 - d$age_Ma) * 10
#ages.sSD = d$age_uncert*10
#ages = unique(ages)
#ages = sort(ages, decreasing = TRUE)
#age.ind = match(d$age_Ma, ages)

##Data to pass to BUGS model
dat = list(pco2.age = ages.s, pco2.age.pre = 1 / ages.sSD^2, al = ages.len,
           pco2 = d$ln_CO2mean, pco2.pre = 1 / (d$ln_2sig/2)^2)

##Parameters to save
parameters = c("pco2_m", "pco2_m.pre", "pco2_m.eps.ac")

##Run it
n.iter = 51000
n.burnin = 1000
n.thin = 5
pt = proc.time()
#p = jags(model.file = "parametric_model_ages.R", parameters.to.save = parameters, 
#         data = dat, inits = NULL, n.chains=3, n.iter = n.iter, 
#         n.burnin = n.burnin, n.thin = n.thin)
p = do.call(jags.parallel, list(model.file = "parametric_model_ages.R", parameters.to.save = parameters, 
                                      data = dat, inits = NULL, n.chains=3, n.iter = n.iter, 
                                      n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

R2jags::traceplot(p, varname = c("pco2_m"[5], "pco2_m[105]", "pco2_m[305]", "pco2_m[505]", "pco2_m[605]", 
                                 "pco2_m[705]", "pco2_m[905]", "pco2_m[1105]", "pco2_m[1205]", "pco2_m.pre", 
                                 "pco2_m.eps.ac"))
View(p$BUGSoutput$summary)

sl = p$BUGSoutput$sims.list
su = p$BUGSoutput$summary
sims = nrow(sl$pco2_m)

png("../jpi_simple.png", width = 8, height = 6, units = "in", res = 600)
plot(-10, 0, xlab="Age (Ma)", ylab ="pCO2", xlim=c(0,68), ylim=c(100,4000))
for(i in seq(1, sims, by = max(floor(sims / 500),1))){
  lines(ages, exp(sl$pco2_m[i,]), col = rgb(0,0,0, 0.02))
}

arrows(d$age_Ma, exp(d$ln_CO2mean+d$ln_2sig), d$age_Ma, exp(d$ln_CO2mean-d$ln_2sig), 
       length = 0, angle = 90, code = 3, col = "light blue")
arrows(d$age_Ma - d$age_uncert, exp(d$ln_CO2mean), d$age_Ma + d$age_uncert, 
       exp(d$ln_CO2mean), length = 0, angle = 90, code = 3, col = "light blue")

lines(ages, exp(su[1:ages.len + 1, 5]), col="red", lw=2)
lines(ages, exp(su[1:ages.len + 1, 3]), col="red", lty=3, lw=2)
lines(ages, exp(su[1:ages.len + 1, 7]), col="red", lty=3, lw=2)

points(d$age_Ma, exp(d$ln_CO2mean), pch=21, bg="white", cex=0.5)

dev.off()

pout = data.frame("Age" = ages, "Mean" = exp(su[1:ages.len + 1, 5]), "ptile_2.5" = exp(su[1:ages.len + 1, 3]), "ptile_97.5" = exp(su[1:ages.len + 1, 7]))
write.csv(pout, "../pCO2_JPI.csv", row.names = FALSE)

#Make space
mod.p = double()

#Calculate the CDF and find quantile value of modern median in it
for(j in 1:length(ages)){
  cdf = ecdf(exp(sl$pco2_m[, j]))
  mod.p[j] = cdf(408.5)
}

mod.pp = pmin(mod.p, 1 - mod.p) * 2

png("../modprog.png", width = 8, height = 6, units = "in", res = 600)
plot(ages, mod.pp, type = "l", ylab = "p(408.5)", xlab = "Age (Ma)", ylim=c(0, 0.2), xlim=c(0,30))
abline(0.05, 0, col = "red", lty = 3)
dev.off()
