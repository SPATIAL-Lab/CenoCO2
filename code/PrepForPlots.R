# Load packages/code
library(openxlsx)
source("code/Helpers.R")

# Prep proxy dataset
dat = prepit("230902_proxies.xlsx")

# Read curves
cp = read.csv("out/500kyrCO2.csv")
cp1M = read.csv("out/1MyrCO2.csv")
cp100k = read.csv("out/100kyrCO2.csv")
cpm = read.csv("out/500kyrCO2MarOnly.csv")
tp = read.csv("out/500kyrTemp.csv")

# Trim to only Cenozoic and modify for ESS plot
cp.c = cp[-(1:4), -c(3,5)]
cp.c[,-1] = (cp.c[,-1] - log(280)) / log(2)
tp.c = tp[-(1:4), -c(3,5)]

# Timescale and colors
cols = rev(rgb(matrix(c(249, 169, 112, 252, 188, 134, 254, 219, 171,
                        255, 242, 0, 255, 249, 174, 254, 242, 227),
                      ncol = 3, byrow = TRUE), maxColorValue = 255))
epochs = c(0, 2.58, 5.33, 23, 33.9, 56)

# Marine temperature proxy data
tdat = read.xlsx("data/Westerhold.xlsx", sheet = "data")
tdat = tdat[,c(1,3)]
names(tdat) = c("age", "temp")

# Ring dataset
tring = read.csv("out/RingTemp.csv")

