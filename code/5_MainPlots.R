# Fig 2 ----
source("code/4_PrepForPlots.R")

# Plotting function
rp = function(){
  par(mai = c(0.1, 1.1, 1.1, 0.9))
  plot(-10, 0, ylab = "", xlab="Age (Ma)",  
       xlim=c(70,0), ylim=c(3.5,8.3), axes = FALSE)
  
  sc = rgb2hsv(col2rgb("dodgerblue2"))
  points(dat$pco2.age, dat$pco2, cex=0.5, 
         col = hsv(sc[1], sc[2]/3, sc[3]))
  
  tsdens(cp, "dodgerblue4")
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
  
  tsdens(tp, "black")
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

#png("out/CenozoicCO2.png", width = 9, height = 5.5, units = "in", res = 600)
setEPS()
postscript("out/Fig2.eps")
rp()
dev.off()

# Fig 3 ----

library(RColorBrewer)
cols = brewer.pal(6, "YlOrRd")

# Assign points to Epoch
ci = findInterval(cp.c$ages, epochs)
tringi = findInterval(tring$Age_mean, epochs)

# Plot it
png("out/Fig3.png", width = 6, height = 7, units = "in", res = 600)
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

arrows(cp.c[,2], tp.c[,3], cp.c[,4], tp.c[,2], length = 0, 
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
