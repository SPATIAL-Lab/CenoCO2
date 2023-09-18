source("code/PrepForPlots.R")

# Fig 2 ----

## Plotting function
rp = function(){
  par(mai = c(0.1, 1.1, 1.1, 0.9))
  plot(-10, 0, ylab = "", xlab="Age (Ma)",  
       xlim=c(70,0), ylim=c(3.5,8.3), axes = FALSE)
  
  sc = rgb2hsv(col2rgb("dodgerblue2"))
  points(dat$pco2.age, dat$pco2, cex=0.5, 
         col = hsv(sc[1], sc[2]/3, sc[3]))
  
  tsdens.s(cp, "dodgerblue4")
  axis(2, c(log(100), log(250), log(500), log(1000), log(2000)),
       c(100, 250, 500, 1000, 2000))
  axis(3, seq(70, 0, by = -10))
  mtext(expression("CO"[2]*" (ppm)"), 2, line = 3, at = 6.2)
  mtext("Age (Ma)", 3, line = 3)
  
  ptop = par("usr")[4]
  enames = c("Ple", "Pli", "Miocene", "Oligocene", "Eocene", "Paleocene")
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
  
  tsdens.s(tp, "black")
  axis(4, seq(-5, 20, by=5), pos = -0.15)
  mtext("GMST (K, relative to preindustrial)", 4, line = 2, at = 7.5)
  
  rcol = col2rgb("grey40", TRUE)
  rcol[4] = 255
  for(i in 1:nrow(tring)){
    lines(c(tring$Age_min[i], tring$Age_max[i]), rep(tring$T_5[i], 2), 
          lw = 2, lend = 1)    
    polygon(c(rep(tring$Age_min[i], 2), rep(tring$Age_max[i], 2)), 
            c(tring$T_025[i], rep(tring$T_975[i], 2), tring$T_025[i]),
            col = rgb(rcol[1], rcol[2], rcol[3], rcol[4], maxColorValue = 255), 
            border = "grey60")
  }
  
}

## png("out/CenozoicCO2.png", width = 9, height = 5.5, units = "in", res = 600) ##
setEPS()
postscript("out/main_figs/Fig2.eps")
rp()
dev.off()

# Fig 3 ----

library(RColorBrewer)
cols = brewer.pal(6, "YlOrRd")

## Assign points to Epoch
ci = findInterval(cp.c$ages, epochs)
tringi = findInterval(tring$Age_mean, epochs)

## Plot it
png("out/main_figs/Fig3.png", width = 6, height = 7, units = "in", res = 600)
par(mai = c(2, 1, 0.2, 0.2))
plot(cp.c[,3], tp.c[,2], xlim = range(cp.c[,-1]), ylim = range(tp.c[,-1]),
     xlab = expression("CO"[2]*" doublings (relative to preindustrial)"),
     ylab = expression(Delta*" GMST (K relative to preindustrial)"))
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
       col = "grey80")

abline(0, 8, lty = 3, col = "white")
abline(10, 8, lty = 3, col = "white")
abline(-10, 8, lty = 3, col = "white")
abline(0, 5, lty = 2, col = "white")
abline(10, 5, lty = 2, col = "white")
abline(-10, 5, lty = 2, col = "white")

arrows(cp.c[,2], tp.c[,3], cp.c[,4], tp.c[,3], length = 0, 
       lwd = 0.75, col = "grey40")
arrows(cp.c[,3], tp.c[,2], cp.c[,3], tp.c[,4], length = 0, 
       lwd = 0.75, col = "grey40")
lines(cp.c[,3], tp.c[,3], lwd = 2, col = "grey40")
points(cp.c[,3], tp.c[,3], pch = 21, bg = cols[ci], 
       col = "grey40", cex = 1.25)

arrows(tring$C_025[-8], tring$T_5[-8], tring$C_975[-8], tring$T_5[-8], length = 0,
       lwd = 0.75)
arrows(tring$C_5[-8], tring$T_025[-8], tring$C_5[-8], tring$T_975[-8], length = 0,
       lwd = 0.75)
lines(tring$C_5[-8], tring$T_5[-8], lwd = 2)
points(tring$C_5[-8], tring$T_5[-8], pch = 22, bg = cols[tringi], cex = 2.5)

text(tring$C_5[-8], tring$T_5[-8], round(tring$Age_mean[-8]), cex = 0.75)

text(cp.c[13, 3], tp.c[13, 3] - 0.3, "60", pos = 2, 
     col = "grey40")
text(cp.c[33, 3], tp.c[33, 3], "50", pos = 3, offset = 0.7, 
     col = "grey40")
text(cp.c[53, 3], tp.c[53, 3] - 0.3, "40", pos = 4, offset = 1,
     col = "grey40")
text(cp.c[73, 3], tp.c[73, 3], "30", pos = 1, offset = 1, 
     col = "grey40")
text(cp.c[93, 3] + 0.1, tp.c[93, 3], "20", pos = 3, offset = 1.2, 
     col = "grey40")
text(cp.c[113, 3], tp.c[113, 3], "10", pos = 2, offset = 2.7, 
     col = "grey40")

text(1.52, -2, "5 \u00B0C/doubling", col = "white",
     srt = (sin(4.8 / (diff(par("usr")[3:4]) / diff(par("usr")[1:2]))))/pi*180)
text(0.95, -1.8, "8 \u00B0C/doubling", col = "white",
     srt = (sin(8 / (diff(par("usr")[3:4]) / diff(par("usr")[1:2]))))/pi*180)
legend("bottomright", legend = c("Pleistocene", "Pliocene", 
                                 "Miocene", "Oligocene", 
                                 "Eocene", "Paleocene"),
       pt.bg = cols, pch = 21, bg = "grey80")
axis(1, at = (log(c(250, 500, 1000, 1500)) - log(280)) / log(2), 
     labels = c("250", "500", "1000", "1500"), line = 5)
mtext(expression("CO"[2]*" (ppm)"), 1, line = 8)
box()
dev.off()

# Print figure ----

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
