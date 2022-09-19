##climate sense
load("out/postCenoLERAM.rda")
cp = p$BUGSoutput$sims.list$pco2_m
cp = cp[,-(ncol(cp))]
load("out/postTemp.rda")
tp = p$BUGSoutput$sims.list$t_m
tp = tp[,-(ncol(tp))]
ages = ages[-length(ages)]
source("code/tsdens.R")
library(openxlsx)

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
plot(-10, 0, ylab = expression("CO"[2]*" (ppmv)"), xlab = "", 
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
  mod.p[j] = cdf(417.58)
}

mod.pp = 1 - mod.p

##Number of curves exceeding
curve.p = double()
cp.ex = exp(cp) > 417.58
for(i in 1:ncol(cp.ex)){
  if(ncol(cp.ex) - i > 0){
    ct = apply(cp.ex[, i:ncol(cp.ex)], 1, any)
  }else{
    ct = any(cp.ex[, ncol(cp.ex)])
  }
  curve.p[i] = sum(ct) / nrow(cp.ex)
}

png("out/modprob.png", width = 9, height = 5, units = "in", 
    res = 600)
par(mai = c(1.1, 1.1, 0.1, 0.8))
plot(0, 0, xlim = c(16, 0), ylim = c(100, 600), type = "n",
     axes = FALSE, xlab = "Age (Ma)", ylab = "")
lines(c(30, 0), rep(417.58, 2), col = "red", lty = 2, lwd = 2)
for(i in 1:250){
  lines(ages, exp(cp[runif(1, 1, nrow(cp)),]), 
        col = rgb(0, 0, 0, 0.05), lwd=2)
}
axis(1)
axis(2, at = c(200, 300, 400, 500, 600))
mtext(expression("CO"[2]*" (ppmv)"), 2, line = 3, at = 400)

par(new = TRUE)

plot(ages, curve.p ^ (1/3), type = "l", axes = FALSE, lwd = 2,
     ylim=c(0, 2), xlim=c(16,0), xlab = "", ylab = "")
lines(c(30, 0), rep(0.05 ^ (1/3), 2), col = "red", lty = 3, lwd = 2)
axis(4, c(0.0001, 0.05, 0.5) ^ (1/3),
     c(expression("10"^{-4}), 0.05, 0.5), pos = 0)
mtext(expression("P(CO"[2]*" > 417.58)"), 4, line = 2,
      at = 0.05 ^ (1/3))
dev.off()

#only Cenozoic
cp.c = cp[,-(1:8)]
tp.c = tp[,-(1:8)]
ages.c = ages[-(1:8)]

#stats
cps = (apply(cp.c, 2, quantile, probs = c(0.025, 0.5, 0.975)) - log(280)) / log(2)
tps = (apply(tp.c, 2, quantile, probs = c(0.025, 0.5, 0.975)))

#colors
cols = colorRampPalette(c("blue", "red"))
cols = cols(6)
epochs = c(0, 2.58, 5.33, 23, 33.9, 56)
ci = findInterval(ages.c, epochs)

#doublings vs T
png("out/CimSens.png", width = 6, height = 7, units = "in", res = 600)
par(mai = c(2, 1, 0.2, 0.2))
plot(cps[2,], tps[2,], xlim = range(cps), ylim = range(tps),
     xlab = expression("CO"[2]*" doublings"),
     ylab = expression(Delta*" GMST"))
abline(0, 8, lty = 3)
abline(0, 5, lty = 3)
lines(cps[2,], tps[2,])
arrows(cps[1,], tps[2,], cps[3,], tps[2,], length = 0, col = "dark grey")
arrows(cps[2,], tps[1,], cps[2,], tps[3,], length = 0, col = "dark grey")
points(cps[2,], tps[2,], pch = 21, bg = cols[ci], cex = 1.25)
text(cps[2, 13], tps[2, 13], "60 Ma", pos = 1, offset = 0.7)
text(cps[2, 33], tps[2, 33], "50 Ma", pos = 3, offset = 0.7)
text(cps[2, 53], tps[2, 53] - 0.3, "40 Ma", pos = 4)
text(cps[2, 73], tps[2, 73], "30 Ma", pos = 1, offset = 0.7)
text(cps[2, 93], tps[2, 93], "20 Ma", pos = 3, offset = 1.2)
text(cps[2, 113], tps[2, 113], "10 Ma", pos = 2, offset = 1.3)
text(2.4, 11.4, "5 \u00B0C/doubling",
     srt = (sin(5 / (diff(par("usr")[3:4]) / diff(par("usr")[1:2]))))/pi*180)
text(0.8, 7, "8 \u00B0C/doubling", 
     srt = (sin(8 / (diff(par("usr")[3:4]) / diff(par("usr")[1:2]))))/pi*180)
legend("bottomright", legend = c("Pleistocene", "Pliocene", 
                                 "Miocene", "Oligocene", 
                                 "Eocene", "Paleocene"),
       pt.bg = cols, pch = 21, bty = "n")
axis(1, at = (log(c(250, 500, 1000, 1500)) - log(280)) / log(2), 
     labels = c("250", "500", "1000", "1500"), line = 5)
mtext(expression("CO"[2]*" (ppmv)"), 1, line = 8)
dev.off()
       
#version for talks
png("out/CenozoicCO2_slide.png", width = 9, height = 6, units = "in", res = 600)

par(mai = c(1.1, 1.1, 0.1, 0.9))
plot(-10, 0, ylab = "", xlab="Age (Ma)",  
     xlim=c(65,0), ylim=c(100,3500), axes = FALSE)

sc = rgb2hsv(col2rgb("dodgerblue"))
points(pco2.age, exp(pco2), cex=0.5, 
       col = hsv(sc[1], sc[2]/3, sc[3]))

tsdens(cbind(ages, exp(pts)), "dodgerblue4")
axis(2, seq(0, 2500, by = 500), c(0, 500, 1000, 1500, 2000, 2500))
axis(1)
mtext(expression("CO"[2]*" (ppmv)"), 2, line = 3, at = 1250)

par(new = TRUE)
sc = rgb2hsv(col2rgb("firebrick4"))
plot(tdat, xlim=c(65,0), ylim = c(-20, 22), axes = FALSE,
     cex = 0.5, col = hsv(sc[1], sc[2]/10, 0.8),
     xlab = "", ylab = "")

tsdens(cbind(ages, tpts), "firebrick4")
axis(4, seq(-5, 20, by=5), pos = -0.15)
mtext("GMST (relative to preindustrial)", 4, line = 2, at = 7.5)

dev.off()
