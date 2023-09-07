# Warming stripes ----
## Load data
source("code/PrepForPlots.R")

## Colors
wscols = c("#08306b", "#08519c", "#2171b5", "#4292c6", "#6baed6",
           "#9ecae1", "#c6dbef", "#deebf7", "#fee0d2", "#fcbba1",
           "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15",
           "#67000d")
cp.trunc = cp[7:nrow(cp),]

tcol = floor((log(tp$X50.+10) - min(log(tp$X50.+10))) / diff(range(log(tp$X50.+10))) * 
               (length(wscols) - 1e-6) + 1)

## Print figure
png("out/main_figs/PrintFig.png", width = 6, height = 2.9, units = "in", res = 600)
par(mai = c(0.5, 0, 0, 0.6))
plot(0, 0, type = "n", xlim = c(65, 0), ylim = c(0, 1), axes = FALSE,
     xlab = "", ylab = "")
### Warming stripes
for(i in 7:nrow(tp)){
  polygon(c(rep(tp$ages[i] - 0.25, 2), rep(tp$ages[i] + 0.25, 2)),
          c(0, 1, 1, 0), col = wscols[tcol[i]], border = NA)
}
### CO2
lines(cp.trunc$ages, (cp.trunc$X50. - min(cp.trunc$X50.)) / diff(range(cp.trunc$X50.))
      * 0.9 + 0.05, lw = 3, col = "white")
lines(cp.trunc$ages, (cp.trunc$X50. - min(cp.trunc$X50.)) / diff(range(cp.trunc$X50.))
      * 0.9 + 0.05, lw = 2)
### Axes
polygon(c(65, 65, 0, 0),
        c(1, 0, 0, 1))
axis(1, c(51, 33.9, 16, 2.6), labels = FALSE, pos = 0)
mtext(c(51, 33.9, 16, 2.6), 1, at = c(51, 33.9, 16, 2.6), cex = 0.7)
mtext("Millions of years before present", 1, 1, cex = 0.9)
ticks = c((log(270) - min(cp.trunc$X50.)) / diff(range(cp.trunc$X50.)) * 0.9 + 0.05,
          (log(480) - min(cp.trunc$X50.)) / diff(range(cp.trunc$X50.)) * 0.9 + 0.05,
          (log(720) - min(cp.trunc$X50.)) / diff(range(cp.trunc$X50.)) * 0.9 + 0.05,
          (log(1600) - min(cp.trunc$X50.)) / diff(range(cp.trunc$X50.)) * 0.9 + 0.05)
axis(4, ticks, labels = FALSE, pos = 0)
mtext(c(270, 480, 720, 1600), 4, -0.4, at = ticks, cex = 0.7)
mtext(expression("Atmospheric CO"[2]*" (ppm)"), 4, 1, cex = 0.9)
### Legend
legx.min = 61
legx.max = legx.min - 20
legy.min = 0.18
legy.max = legy.min + 0.07
n = length(wscols)
xint = (legx.max - legx.min) / n
polygon(c(legx.min, legx.min, legx.max, legx.max),
        c(legy.min, legy.max, legy.max, legy.min), lwd = 4, 
        border = "white", lend = 2, ljoin = 1)
for(i in seq(n)){
  polygon(c(rep(legx.min + xint * (i - 1), 2), rep(legx.min + xint * i, 2)),
          c(legy.min, legy.max, legy.max, legy.min), border = NA,
          col = wscols[i])
}
tticks = seq(legx.min, legx.max, length = 4)
tvals = min(log(tp$X50.+10)) + diff(range(log(tp$X50.+10))) * seq(0, 1, length = 4)
tvals = exp(tvals) - 10
for(i in 1:4){
  lines(rep(tticks[i], 2), c(legy.min, legy.min - 0.025), col = "white", 
        lw = 2, lend = 1)
  text(x = tticks[i], y = legy.min - 0.055, round(tvals[i], 1), col = "white", 
       cex = 0.7)
}
text(mean(tticks), legy.min - 0.11, "Temperature (\u00B0C)", cex = 0.7, 
     col = "white")
dev.off()

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
