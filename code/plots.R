##climate sense
load("out/postCeno.rda")
cp = p$BUGSoutput$sims.list$pco2_m
cp = cp[,-(ncol(cp))]
load("out/postTemp.rda")
tp = p$BUGSoutput$sims.list$t_m
tp = tp[,-(ncol(tp))]
ages = ages[-length(ages)]
source("code/tsdens.R")

#timeseries plot
pts = apply(cp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
pts = t(pts)

tpts = apply(tp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
tpts = t(tpts)

#Temperature proxy data
tdat = read.xlsx("data/Westerhold.xlsx", sheet = "data")

#Data subset 
tdat = tdat[,c(1,3)]
names(tdat) = c("age", "temp")

png("out/CenozoicCO2.png", width = 8, height = 11, units = "in", res = 600)
#cairo_ps("out/CenozoicCO2.eps", width = 8, height = 6,
#         fallback_resolution = 600)
layout(matrix(c(1, 2), nrow = 2), heights = c(lcm(5*2.54), lcm(6*2.54)))
par(mai = c(0.1, 1.1, 0.1, 0.1))
plot(-10, 0, ylab = expression("CO"[2]*" (ppmv)"), 
     xlim=c(65,0), ylim=c(100,3000), axes = FALSE)

arrows(pco2.age, exp(pco2 + 2 * pco2.sd), 
       pco2.age, exp(pco2 - 2 * pco2.sd), 
       length = 0, angle = 90, code = 3, col = "light grey")
arrows(pco2.age - pco2.age.sd, exp(pco2), pco2.age + pco2.age.sd, 
       exp(pco2), length = 0, angle = 90, code = 3, col = "light grey")

points(pco2.age, exp(pco2), cex=0.5, col = "dark grey")

tsdens(cbind(ages, exp(pts)), "dodgerblue4")
axis(2)

par(mai = c(1.1, 1.1, 0.1, 0.1))
plot(-10, 0, xlab="Age (Ma)", ylab = "GMST (relative to preindustrial)", 
     xlim=c(65,0), ylim = range(tdat$temp), axes = FALSE)

points(tdat, cex=0.5, col = "dark grey")

tsdens(cbind(ages, tpts), "firebrick4")
axis(1)
axis(2)

dev.off()

#write output as CSV
pout = data.frame("Age" = ages, "Mean" = exp(su[1:ages.len + 1, 5]), "ptile_2.5" = exp(su[1:ages.len + 1, 3]), "ptile_97.5" = exp(su[1:ages.len + 1, 7]))
write.csv(pout, "out/pCO2_timeseries.csv", row.names = FALSE)

#Change vs modern
##Make space
mod.p = double()

##Calculate the CDF and find quantile value of modern median in it
for(j in 1:length(ages)){
  cdf = ecdf(exp(cp[, j]))
  mod.p[j] = cdf(418)
}

mod.pp = 1 - mod.p

png("out/modprog.png", width = 8, height = 6, units = "in", res = 600)
plot(ages, mod.pp, type = "l", ylab = "p(418)", xlab = "Age (Ma)", ylim=c(0, 0.2), xlim=c(0,30))
abline(0.05, 0, col = "red", lty = 3)
dev.off()

##climate sense
load("out/postCeno.rda")
cp = p$BUGSoutput$sims.list$pco2_m
cp = cp[,-(ncol(cp))]
load("out/postTemp.rda")
tp = p$BUGSoutput$sims.list$t_m
tp = tp[,-(ncol(tp))]

#stats
cps = (apply(cp, 2, quantile, probs = c(0.025, 0.5, 0.975)) - log(280)) / log(2)
tps = (apply(tp, 2, quantile, probs = c(0.025, 0.5, 0.975)))

#doublings vs T
png("out/CimSens.png", width = 8, height = 6, units = "in", res = 600)
cols = colorRampPalette(c("red", "blue"))
plot(cps[2,], tps[2,], xlim = range(cps), ylim = range(tps),
     xlab = expression("CO"[2]*" doublings (relative to preindustrial)"),
     ylab = "GMST (relative to preindustrial)")
abline(0, 8, lty = 3)
abline(0, 5, lty = 3)
lines(cps[2,], tps[2,])
arrows(cps[1,], tps[2,], cps[3,], tps[2,], length = 0, col = "dark grey")
arrows(cps[2,], tps[1,], cps[2,], tps[3,], length = 0, col = "dark grey")
points(cps[2,], tps[2,], pch = 21, bg = cols(length(cps[2,])))
dev.off()
       