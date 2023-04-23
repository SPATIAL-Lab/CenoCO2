# MCMC diagnostics ----
library(R2jags)
load("bigout/postCenoLERAM.rda")

# Traceplots
R2jags::traceplot(p, varname = "pco2_m")
R2jags::traceplot(p, varname = c("pco2_m.pre", "pco2_m.eps.ac"))
dev.off()

# Summary
View(p$BUGSoutput$summary)

# Process and save data ----
source("code/8_Helpers.R")
## 500 kyr CO2 ----
# Set up ages vector
ages.bin = 0.5
ages.max = 68
ages = agevec(ages.max, ages.bin)
ages.len = length(ages)
ages = ages[-length(ages)]

#Load data
cp = p$BUGSoutput$sims.list$pco2_m
cp = cp[,-(ncol(cp))]
cp = cp[,-(1:4)]

# Stats for timeseries plot
pts = apply(cp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
pts = t(pts)
pts = cbind(ages, pts)

# Save
write.csv(pts, "out/500kyrCO2.csv", row.names = FALSE)

## 500 kyr Temp ----
load("bigout/postTemp.rda")
tp = p$BUGSoutput$sims.list$t_m
tp = tp[,-(ncol(tp))]

# Trim posterior ts to get only 68 Ma to present
tp = tp[,-(1:4)]

# Stats for timeseries plot
tpts = apply(tp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
tpts = t(tpts)
tpts = cbind(ages, tpts)

# Save
write.csv(tpts, "out/500kyrTemp.csv", row.names = FALSE)

## Ring temperature dataset ----
tring = data.frame("Age_min" = c(48, 42, 33.9, 27.8, 20.3, 14.7, 7.2, 3),
                   "Age_max" = c(55, 46, 37.8, 33.9, 23.0, 17.0, 11.6, 3.3),
                   "T_025" = c(10.5, 8.2, 7.7, 5.6, 5.4, 6.4, 4.5, 3),
                   "T_5" = c(12.5, 10.5, 9.6, 7.5, 7.1, 8.2, 5.9, 3.9),
                   "T_975" = c(14.5, 12.8, 11.5, 9.4, 8.8, 10, 7.3, 4.8))
tring$Age_mean = apply(tring[,1:2], 1, mean)
tring$Bin_min = round((ages.max - tring$Age_min) / ages.bin)
tring$Bin_max = round((ages.max - tring$Age_max) / ages.bin)
tring$Bin_mean = round((ages.max - tring$Age_mean) / ages.bin)
tring$C_975 = tring$C_5 = tring$C_025 = rep(0)

for(i in 1:nrow(tring)){
  tring[i, 10:12] = (quantile(cp[,tring$Bin_max[i]:tring$Bin_min[i]], 
                              probs = c(0.025, 0.5, 0.975)) - log(280)) / log(2)
}

write.csv(tring, "out/RingTemp.csv")

## 1 Myr CO2 ----
load("bigout/postCeno1Myr.rda")
cp = p$BUGSoutput$sims.list$pco2_m
cp = cp[,-(ncol(cp))]

# Set up ages vector
ages.bin = 1
ages.max = 68
ages = agevec(ages.max, ages.bin)
ages.len = length(ages)
ages = ages[-length(ages)]

#trim posterior ts
cp = cp[,-(1:2)]

# Stats for timeseries plot
pts = apply(cp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
pts = t(pts)
pts = cbind(ages, pts)

# Save
write.csv(pts, "out/1MyrCO2.csv", row.names = FALSE)

## 100 kyr CO2 ----
# Set up ages vector
ages.bin = 0.1
ages.max = 68
ages = agevec(ages.max, ages.bin)
ages.len = length(ages)
ages = ages[-length(ages)]

# Load data
load("bigout/postCeno100kyr.rda")
cp = p$BUGSoutput$sims.list$pco2_m
cp = cp[,-(ncol(cp))]

#trim posterior ts
cp = cp[,-(1:20)]

# Stats for timeseries plot
pts = apply(cp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
pts = t(pts)
pts = cbind(ages, pts)

# Save
write.csv(pts, "out/100kyrCO2.csv", row.names = FALSE)
