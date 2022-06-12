##My local working directories
setwd("C:/Users/gjbowen/Dropbox/HypoMirror/Honish_pCO2/code/")

library(openxlsx)
library(raster)

##Read proxy data
d = read.xlsx("../191126_proxies.xlsx", sheet = 2)

csynth = matrix(nrow = nrow(d), ncol = 1000)
tsynth = matrix(nrow = nrow(d), ncol = 1000)
for(i in 1:nrow(d)){
  csynth[i,] = rnorm(1000, d$ln_CO2mean[i], d$ln_2sig[i]/2)
  tsynth[i,] = rnorm(1000, d$age_Ma[i], d$age_uncert[i])
}
csynth = exp(csynth)

dt = 0.5
dc = 50
t = seq(0 + dt/2, 67 - dt/2, by = dt)
c = seq(0+ dc/2, 5000 - dc/2, by = dc)

dd = matrix(rep(0), nrow = length(c), ncol = length(t))

for(i in 1:length(t)){
  for(j in 1:length(c)){
    dd[j, i] = sum(csynth >= (c[j] - dc/2) & csynth < (c[j] + dc/2) &
                     tsynth >= (t[i] - dt/2) & tsynth < (t[i] + dt/2))
  }
}

r = raster(dd, xmn = 0, xmx = 67, ymn = 0, ymx = 5000)
r = flip(r, direction = "y")

png("heatMap.png", width = 8, height = 5, units = "in", res = 600)
image(log(r), xlim = c(67,0), xlab = "Age (Ma)", ylab = "pCO2")
box()
dev.off()
