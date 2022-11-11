load("out/postTemp.rda")
tp = p$BUGSoutput$sims.list$t_m
tp = tp[,-(ncol(tp))]
load("out/postCenoLERAM.rda")
cp = p$BUGSoutput$sims.list$pco2_m
cp = cp[,-(ncol(cp))]
source("code/helpers.R")
library(openxlsx)

# prep data
dat = prepit()

# Set up ages vector
ages.bin = 0.5
ages = agevec(70, ages.bin)
ages.len = length(ages)
ages = ages[-length(ages)]

# timescale and colors
cols = rev(rgb(matrix(c(249, 169, 112, 252, 188, 134, 254, 219, 171,
                        255, 242, 0, 255, 249, 174, 254, 242, 227),
                      ncol = 3, byrow = TRUE), maxColorValue = 255))
epochs = c(0, 2.58, 5.33, 23, 33.9, 56)

# timeseries plot
pts = apply(cp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
pts = t(pts)

tpts = apply(tp, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
tpts = t(tpts)

# Temperature proxy data
tdat = read.xlsx("data/Westerhold.xlsx", sheet = "data")

# Data subset 
tdat = tdat[,c(1,3)]
names(tdat) = c("age", "temp")

# only Cenozoic
cp.c = cp[,-(1:8)]
tp.c = tp[,-(1:8)]
ages.c = ages[-(1:8)]

# stats
cps = (apply(cp.c, 2, quantile, probs = c(0.025, 0.5, 0.975)) - log(280)) / log(2)
tps = (apply(tp.c, 2, quantile, probs = c(0.025, 0.5, 0.975)))

# Ring dataset
tring = data.frame("Age_min" = c(48, 42, 33.9, 27.8, 20.3, 14.7, 7.2, 3),
                   "Age_max" = c(55, 46, 37.8, 33.9, 23.0, 17.0, 11.6, 3.3),
                   "T_025" = c(10.5, 8.2, 7.7, 5.6, 5.4, 6.4, 4.5, 3),
                   "T_5" = c(12.5, 10.5, 9.6, 7.5, 7.1, 8.2, 5.9, 3.9),
                   "T_975" = c(14.5, 12.8, 11.5, 9.4, 8.8, 10, 7.3, 4.8))
tring$Age_mean = apply(tring[,1:2], 1, mean)
tring$Bin_min = pmax(pmin(round((66 - tring$Age_min) / ages.bin), length(ages.c)), 1)
tring$Bin_max = pmax(pmin(round((66 - tring$Age_max) / ages.bin), length(ages.c)), 1)
tring$Bin_mean = pmax(pmin(round((66 - tring$Age_mean) / ages.bin), length(ages.c)), 1)
tring$C_975 = tring$C_5 = tring$C_025 = rep(0)

for(i in 1:nrow(tring)){
  tring[i, 10:12] = (quantile(cp.c[,tring$Bin_max[i]:tring$Bin_min[i]], 
                              probs = c(0.025, 0.5, 0.975)) - log(280)) / log(2)
}

# assign points to Epoch
ci = findInterval(ages.c, epochs)
tringi = findInterval(tring$Age_mean, epochs)

rp = function(){
  par(mai = c(0.1, 1.1, 1.1, 0.9))
  plot(-10, 0, ylab = "", xlab="Age (Ma)",  
       xlim=c(70,0), ylim=c(3.5,8.3), axes = FALSE)
  
  sc = rgb2hsv(col2rgb("dodgerblue2"))
  points(dat$pco2.age, dat$pco2, cex=0.5, 
         col = hsv(sc[1], sc[2]/3, sc[3]))
  
  tsdens(cbind(ages, pts), "dodgerblue4")
  axis(2, c(log(100), log(250), log(500), log(1000), log(2000)),
       c(100, 250, 500, 1000, 2000))
  axis(3, seq(70, 0, by = -10))
  mtext(expression("CO"[2]*" (ppm)"), 2, line = 3, at = 6.2)
  mtext("Age (Ma)", 3, line = 3)
  
  ptop = par("usr")[4]
  enames = c("Qu", "Pl", "Miocene", "Oligocene", "Eocene", "Paleocene")
  for(i in 1:(length(epochs))){
    polygon(c(rep(c(epochs, 66)[i], 2), rep(c(epochs, 66)[i+1], 2)),
            c(ptop, rep(ptop - 0.3, 2), ptop), col = cols[i])
    text(mean(c(epochs, 66)[i:(i+1)]), ptop - 0.15, enames[i])
  }
  
  par(new = TRUE)
  sc = rgb2hsv(col2rgb("grey50"))
  plot(tdat, xlim=c(70,0), ylim = c(-5, 34), axes = FALSE,
       cex = 0.5, col = hsv(sc[1], sc[2]/10, 0.8),
       xlab = "", ylab = "")
  
  tsdens(cbind(ages, tpts), "black")
  axis(4, seq(-5, 20, by=5), pos = -0.15)
  mtext("GMST (K, relative to preindustrial)", 4, line = 2, at = 7.5)
  
  rcol = col2rgb("grey40", TRUE)
  rcol[4] = 80
  for(i in 1:nrow(tring)){
    lines(c(tring$Age_min[i], tring$Age_max[i]), rep(tring$T_5[i], 2), 
          lw = 2, lend = 1)    
    polygon(c(rep(tring$Age_min[i], 2), rep(tring$Age_max[i], 2)), 
            c(tring$T_025[i], rep(tring$T_975[i], 2), tring$T_025[i]),
            col = rgb(rcol[1], rcol[2], rcol[3], rcol[4], maxColorValue = 255), 
            border = "grey60")
  }
  
}

png("out/CenozoicCO2.png", width = 9, height = 5.5, units = "in", res = 600)
#cairo_ps("out/CenozoicCO2.eps", width = 8, height = 6,
#         fallback_resolution = 600)
rp()
dev.off()

# ESS plot

# doublings vs T
png("out/CimSens.png", width = 6, height = 7, units = "in", res = 600)
par(mai = c(2, 1, 0.2, 0.2))
plot(cps[2,], tps[2,], xlim = range(cps), ylim = range(tps),
     xlab = expression("CO"[2]*" doublings"),
     ylab = expression(Delta*" GMST (\u00B0 C)"))
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
       col = "grey80")
for(i in seq(-25, 15, by = 5)){
  abline(i, 8, lty = 3, col = "white")
}
for(i in seq(-25, 15, by = 5)){
  abline(i, 5, lty = 2, col = "white")
}
arrows(cps[1,], tps[2,], cps[3,], tps[2,], length = 0, 
       lwd = 0.75, col = "grey40")
arrows(cps[2,], tps[1,], cps[2,], tps[3,], length = 0, 
       lwd = 0.75, col = "grey40")
lines(cps[2,], tps[2,], lwd = 2, col = "grey40")
points(cps[2,], tps[2,], pch = 21, bg = cols[ci], 
       col = "grey40", cex = 1.25)

arrows(tring$C_025, tring$T_5, tring$C_975, tring$T_5, length = 0,
       lwd = 0.75)
arrows(tring$C_5, tring$T_025, tring$C_5, tring$T_975, length = 0,
       lwd = 0.75)
lines(tring$C_5, tring$T_5, lwd = 2)
points(tring$C_5, tring$T_5, pch = 22, bg = cols[tringi], cex = 2.5)

text(tring$C_5, tring$T_5, round(tring$Age_mean), cex = 0.75)

text(cps[2, 13], tps[2, 13] - 0.3, "60", pos = 2, 
     col = "grey40")
text(cps[2, 33], tps[2, 33], "50", pos = 3, offset = 0.7, 
     col = "grey40")
text(cps[2, 53], tps[2, 53] - 0.3, "40", pos = 4, offset = 1,
     col = "grey40")
text(cps[2, 73], tps[2, 73], "30", pos = 1, offset = 1, 
     col = "grey40")
text(cps[2, 93] + 0.1, tps[2, 93], "20", pos = 3, offset = 1.2, 
     col = "grey40")
text(cps[2, 113], tps[2, 113], "10", pos = 2, offset = 2.7, 
     col = "grey40")

text(1.52, -2, "5 \u00B0C/doubling", col = "white",
     srt = (sin(4.8 / (diff(par("usr")[3:4]) / diff(par("usr")[1:2]))))/pi*180)
text(0.95, -1.8, "8 \u00B0C/doubling", col = "white",
     srt = (sin(8 / (diff(par("usr")[3:4]) / diff(par("usr")[1:2]))))/pi*180)
legend("bottomright", legend = c("Quaternary", "Pliocene", 
                                 "Miocene", "Oligocene", 
                                 "Eocene", "Paleocene"),
       pt.bg = cols, pch = 21, bg = "grey80")
axis(1, at = (log(c(250, 500, 1000, 1500)) - log(280)) / log(2), 
     labels = c("250", "500", "1000", "1500"), line = 5)
mtext(expression("CO"[2]*" (ppmv)"), 1, line = 8)
box()
dev.off()


# Change vs modern
## Make space
mod.p = double()

## Calculate the CDF and find quantile value of modern median in it
modCO2 = 417.93
for(j in 1:length(ages)){
  cdf = ecdf(exp(cp[, j]))
  mod.p[j] = cdf(modCO2)
}

mod.pp = 1 - mod.p

## Number of curves exceeding
curve.p = double()
cp.ex = exp(cp) > modCO2
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
lines(c(30, 0), rep(modCO2, 2), col = "red", lty = 2, lwd = 2)
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
mtext("P(Exceedance)", 4, line = 2, at = 0.05 ^ (1/3))
dev.off()


# version for talks
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

# Stats for data

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

png("out/ts1.png", 9, 4.25, units = "in", res = 600)
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


png("out/ts2.png", 9, 4.25, units = "in", res = 600)
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

# Show age model uncertainty
sl = p$BUGSoutput$sims.list
sl$pco2_m = sl$pco2_m[, -ncol(sl$pco2_m)]
# Problem point at 42.5 Ma
le.ind = dat$pco2.age == 42.5

png("out/ProblemPoint.png", width = 9, height = 6, units = "in", res = 600)
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
tsdens(cbind(ages, pts), "dodgerblue4")
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

