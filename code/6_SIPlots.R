# Change vs modern ----

# Load data
load("bigout/postCenoLERAM.rda")
cp = p$BUGSoutput$sims.list$pco2_m
cp = cp[,-(ncol(cp))]
cp = cp[,-(1:4)]

# Age vector
source("code/8_Helpers.R")
ages.bin = 0.5
ages.max = 68
ages = agevec(ages.max, ages.bin)
ages.len = length(ages)
ages = ages[-length(ages)]

# Make space
mod.p = double()

# Approach 1, using modern value
modCO2 = 418

# Number of curves exceeding
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

#Plot it
png("out/FigS10_2.png", width = 9, height = 5, units = "in", 
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
mtext(expression("CO"[2]*" (ppm)"), 2, line = 3, at = 400)

par(new = TRUE)

plot(ages, curve.p ^ (1/3), type = "l", axes = FALSE, lwd = 2,
     ylim=c(0, 2), xlim=c(16,0), xlab = "", ylab = "")
lines(c(30, 0), rep(0.05 ^ (1/3), 2), col = "red", lty = 3, lwd = 2)
axis(4, c(0.0001, 0.05, 0.5) ^ (1/3),
     c(expression("10"^{-4}), 0.05, 0.5), pos = 0)
mtext("P(Exceedance)", 4, line = 2, at = 0.05 ^ (1/3))
dev.off()

# Approach 2, using change from preindustrial
cp.anom = cp
for(i in 1:nrow(cp)){
  cp.anom[i,] = exp(cp[i,]) - exp(cp[i, ncol(cp)])
}

modanom = 418 - 280

## Number of curves exceeding
curve.p = double()
cp.ex = cp.anom > modanom
for(i in 1:ncol(cp.ex)){
  if(ncol(cp.ex) - i > 0){
    ct = apply(cp.ex[, i:ncol(cp.ex)], 1, any)
  }else{
    ct = any(cp.ex[, ncol(cp.ex)])
  }
  curve.p[i] = sum(ct) / nrow(cp.ex)
}

#Plot it
png("out/FigS10_3.png", width = 9, height = 5, units = "in", 
    res = 600)
par(mai = c(1.1, 1.1, 0.1, 0.8))
plot(0, 0, xlim = c(16, 0), ylim = c(-150, 350), type = "n",
     axes = FALSE, xlab = "Age (Ma)", ylab = "")
lines(c(30, 0), rep(modanom, 2), col = "red", lty = 2, lwd = 2)
for(i in 1:250){
  lines(ages, cp.anom[runif(1, 1, nrow(cp.anom)),], 
        col = rgb(0, 0, 0, 0.05), lwd=2)
}
axis(1)
axis(2, at = c(0, 100, 200, 300))
mtext(expression("CO"[2]*" change (ppm)"), 2, line = 3, at = 150)

par(new = TRUE)

plot(ages, curve.p ^ (1/3), type = "l", axes = FALSE, lwd = 2,
     ylim=c(0, 2), xlim=c(16,0), xlab = "", ylab = "")
lines(c(30, 0), rep(0.05 ^ (1/3), 2), col = "red", lty = 3, lwd = 2)
axis(4, c(0.0001, 0.05, 0.5) ^ (1/3),
     c(expression("10"^{-4}), 0.05, 0.5), pos = 0)
mtext("P(Exceedance)", 4, line = 2, at = 0.05 ^ (1/3))
dev.off()

# Alt resolutions ----

source("code/4_PrepForPlots.R")

# 1 Myr
png("out/FigS10_1.png", width = 9, height = 6, units = "in", res = 600)

par(mai = c(0.1, 1.1, 1.1, 1.1))
plot(-10, 0, ylab = "", xlab="Age (Ma)",  
     xlim=c(67,0), ylim=c(2.5,8.3), axes = FALSE)

sc = rgb2hsv(col2rgb("dodgerblue2"))
points(dat$pco2.age, dat$pco2, cex=0.5, 
       col = hsv(sc[1], sc[2]/3, sc[3]))

tsdens(cp1M, "dodgerblue4")
axis(2, c(log(100), log(250), log(500), log(1000), log(2000)),
     c(100, 250, 500, 1000, 2000))
axis(3, c(66,0), lwd.ticks = 0, labels = FALSE)
axis(3, seq(60, 0, by = -10))
mtext(expression("CO"[2]*" (ppm)"), 2, line = 3, at = 6.2)
mtext("Age (Ma)", 3, line = 3)

ptop = par("usr")[4]
enames = c("Ple", "Pli", "Miocene", "Oligocene", "Eocene", "Paleocene")
for(i in 1:(length(epochs))){
  polygon(c(rep(c(epochs, 66)[i], 2), rep(c(epochs, 66)[i+1], 2)),
          c(ptop, rep(ptop - 0.3, 2), ptop), col = cols[i])
  text(mean(c(epochs, 66)[i:(i+1)]), ptop - 0.15, enames[i])
}

# 100 kyr
par(new = TRUE)
plot(-10, 0, ylab = "", xlab="Age (Ma)",  
     xlim=c(67,0), ylim=c(4.5,10.3), axes = FALSE)

sc = rgb2hsv(col2rgb("dodgerblue2"))
points(dat$pco2.age, dat$pco2, cex=0.5, 
       col = hsv(sc[1], sc[2]/3, sc[3]))

tsdens(cp100k, "dodgerblue4")
axis(4, c(log(100), log(250), log(500), log(1000), log(2000)),
     c(100, 250, 500, 1000, 2000))
mtext(expression("CO"[2]*" (ppm)"), 4, line = 3, at = 6.2)

dev.off()
