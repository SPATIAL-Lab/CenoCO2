source("code/PrepForPlots.R")

# Timeseries slide for talks ----
png("out/other_figs/CenozoicCO2_slide.png", width = 9, height = 6, units = "in", res = 600)

par(mai = c(1.1, 1.1, 0.1, 0.9))
plot(-10, 0, ylab = "", xlab="Age (Ma)",  
     xlim=c(65,0), ylim=c(100,3500), axes = FALSE)

sc = rgb2hsv(col2rgb("dodgerblue"))
points(dat$pco2.age, exp(dat$pco2), cex=0.5, 
       col = hsv(sc[1], sc[2]/3, sc[3]))

tsdens(cbind(cp[,1], exp(cp[,-1])), "dodgerblue4")
axis(2, seq(0, 2500, by = 500), c(0, 500, 1000, 1500, 2000, 2500))
axis(1)
mtext(expression("CO"[2]*" (ppm)"), 2, line = 3, at = 1250)

par(new = TRUE)
sc = rgb2hsv(col2rgb("firebrick4"))
plot(tdat, xlim=c(65,0), ylim = c(-20, 22), axes = FALSE,
     cex = 0.5, col = hsv(sc[1], sc[2]/10, 0.8),
     xlab = "", ylab = "")

tsdens(tp, "firebrick4")
axis(4, seq(-5, 20, by=5), pos = -0.15)
mtext("GMST (relative to preindustrial)", 4, line = 2, at = 7.5)

dev.off()

# Data density plots ----

## Stats for data
a1 = seq(0.5, 64.5, by = 0.1)
nd1 = np1 = a1
for(i in 1:length(a1)){
  ds = dat[dat$pco2.age < a1[i] + 0.5 & dat$pco2.age >= a1[i] - 0.5,]
  nd1[i] = nrow(ds)
  np1[i] = length(unique(ds$pco2.prox))
}

a5 = seq(2.5, 64.5, by = 0.1)
nd5 = np5 = a5
for(i in 1:length(a5)){
  ds = dat[dat$pco2.age < a5[i] + 2.5 & dat$pco2.age >= a5[i] - 2.5,]
  nd5[i] = nrow(ds)
  np5[i] = length(unique(ds$pco2.prox)) 
}

png("out/other_figs/DataDens1.png", 9, 4.25, units = "in", res = 600)
par(mar = c(5,5,1,1))
plot(a1, log(nd1), xlim = c(67, 0), type = "l", axes = FALSE,
     lwd = 3, xlab = "Age (Ma)", ylab = "# data per Myr",
     col = rgb(237, 125, 49, maxColorValue = 255), cex.lab = 1.5)
axis(1, cex.axis = 1.5)
axis(2, at = log(c(1, 10, 100)), labels = c(1, 10, 100), 
     cex.axis = 1.5)
box()
lines(a5, log(nd5/5), lwd = 3, lty = 3, 
      col = rgb(237, 125, 49, maxColorValue = 255))
legend(67, log(300), c("1 Myr bin", "5 Myr bin"), lty = c(1, 3), 
       lwd = c(2,2), bty = "n", cex = 1.5,
       col = rgb(237, 125, 49, maxColorValue = 255))
dev.off()

png("out/other_figs/DataDens2.png", 9, 4.25, units = "in", res = 600)
par(mar = c(5,5,1,1))
plot(a1, np1, xlim = c(67, 0), type = "l", cex.axis = 1.5,
     ylim = c(0, 5), lwd = 3, xlab = "Age (Ma)", ylab = "# proxies",
     col = rgb(112, 173, 71, maxColorValue = 255), cex.lab = 1.5)
lines(a5, np5, lwd = 3, lty = 3, 
      col = rgb(112, 173, 71, maxColorValue = 255))
legend(67, 5.4, c("1 Myr bin", "5 Myr bin"), lty = c(1, 3), 
       lwd = c(2,2), bty = "n", cex = 1.5,
       col = rgb(112, 173, 71, maxColorValue = 255))
dev.off()

# Illustrate age model uncertainty ----

load("bigout/postCenoLERAM.rda")
sl = p$BUGSoutput$sims.list
sl$pco2_m = sl$pco2_m[, -ncol(sl$pco2_m)]
# Problem point at 42.5 Ma
le.ind = dat$pco2.age == 42.5

png("out/other_figs/Age_uncert.png", width = 9, height = 6, units = "in", res = 600)
par(mai = c(1.1, 1.1, 0.2, 1.1))
plot(dat$pco2.age[le.ind], dat$pco2[le.ind], xlim = c(55, 25), ylim = c(5.5, 9.5),
     pch = 21, bg = "dark gray", axes = FALSE, xlab = ("Age(Ma)"),
     ylab = "")
arrows((dat$pco2.age - dat$pco2.age.sd)[le.ind], dat$pco2[le.ind],
       (dat$pco2.age + dat$pco2.age.sd)[le.ind], dat$pco2[le.ind], length = 0)
arrows(dat$pco2.age[le.ind], (dat$pco2 - dat$pco2.sd)[le.ind],
       dat$pco2.age[le.ind], (dat$pco2 + dat$pco2.sd)[le.ind], length = 0)
points(dat$pco2.age[le.ind], dat$pco2[le.ind], pch = 21, bg = "dark gray", cex = 2)
axis(1)
axis(2, at = c(log(250), log(500), log(1000), log(2000)), labels = c(250, 500, 1000, 2000))
mtext(expression("CO"[2]*" (ppm)"), 2, line = 3, at = 6.5)
tsdens(cp, "dodgerblue4")
set.seed(101)
ri = sample(1:nrow(sl$pco2.ai), 500)
for(i in 1:length(ri)){
  points(sl$pco2.ai[ri[i], le.ind], 
         sl$pco2.off[ri[i], le.ind],
         col = rgb(1, 0, 0, 0.5))
}
par(new = TRUE)
plot(density(rnorm(1e6, dat$pco2.age[le.ind], dat$pco2.age.sd[le.ind])),
     ylim = c(-0.18, 0.14), xlim = c(55, 25), axes = FALSE, xlab = "",
     ylab = "", main = "", zero.line = FALSE)
lines(density(sl$pco2.ai[, le.ind]), col = "red")
axis(4, at = c(0, 0.05, 0.1, 0.15))
mtext("P(Age)", 4, line = 3, at = 0.075)
dev.off()
