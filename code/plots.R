load("out/post.rda")
source("code/tsdens.R")

sl = p$BUGSoutput$sims.list
su = p$BUGSoutput$summary
sims = nrow(sl$pco2_m)

#parameter posteriors
plot(density(sl$pco2_m.eps.ac), xlim = c(0.5, 1), col = "red")
lines(density(runif(1e6, 0.5, 0.95)))

plot(density(sl$pco2_m.pre), col = "red")
lines(density(rgamma(1e6, shape = 1, rate = 0.01)))

#timeseries plot
png("out/jpi_simple.png", width = 8, height = 6, units = "in", res = 600)
plot(-10, 0, xlab="Age (Ma)", ylab ="pCO2", xlim=c(65,0), ylim=c(100,3000))
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

points(pco2.age, exp(pco2), 
       pch=21, bg="white", cex=0.5)

dev.off()

#write output as CSV
pout = data.frame("Age" = ages, "Mean" = exp(su[1:ages.len + 1, 5]), "ptile_2.5" = exp(su[1:ages.len + 1, 3]), "ptile_97.5" = exp(su[1:ages.len + 1, 7]))
write.csv(pout, "out/pCO2_timeseries.csv", row.names = FALSE)

#Change vs modern
##Make space
mod.p = double()

##Calculate the CDF and find quantile value of modern median in it
for(j in 1:length(ages)){
  cdf = ecdf(exp(sl$pco2_m[, j]))
  mod.p[j] = cdf(408.5)
}

mod.pp = pmin(mod.p, 1 - mod.p) * 2

png("out/modprog.png", width = 8, height = 6, units = "in", res = 600)
plot(ages, mod.pp, type = "l", ylab = "p(408.5)", xlab = "Age (Ma)", ylim=c(0, 0.2), xlim=c(0,30))
abline(0.05, 0, col = "red", lty = 3)
dev.off()
