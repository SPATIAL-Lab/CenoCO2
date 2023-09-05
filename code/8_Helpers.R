
# Adds time series to plot w 2 prob density envelopes
tsdens = function(d, base = "black"){
  #Check dimensions of d
  if(ncol(d) != 6){stop("d cols should be should be time, 5%, 25%, 50%, 75%, 95% CI")}
  
  base.rgb = col2rgb(base)
  cols = c(rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 0.25), 
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 0.25),
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 1))
  
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 2], rev(d[, 6])), 
          col = cols[1], border = NA)
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 3], rev(d[, 5])), 
          col = cols[2], border = NA)
  lines(d[, 1], d[, 4], col = cols[3], lwd = 2)
}

# Adds time series to plot w 2 prob density envelopes - solid color
tsdens.s = function(d, base = "black"){
  #Check dimensions of d
  if(ncol(d) != 6){stop("d cols should be should be time, 5%, 25%, 50%, 75%, 95% CI")}
  
  base.rgb = col2rgb(base)
  cols = c(rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255), 
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255),
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255))
  
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 2], rev(d[, 6])), 
          col = cols[1], border = NA)
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 3], rev(d[, 5])), 
          col = cols[2], border = NA)
  lines(d[, 1], d[, 4], col = cols[3], lwd = 2)
}

# Generates age vector
agevec = function(start, ages.bin){
  return(ages = seq(start, 0, by = 0 - ages.bin) - ages.bin / 2)
}

# Prep the data and timeseries
prepit = function(df){
  #Read proxy data
  df = file.path("data", df)
  d = read.xlsx(df, sheet = "all data product")
  
  #Data subset 
  d = d[,c("CO2_ppm", "CO2_uncertainty_pos_ppm", "CO2_uncertainty_neg_ppm",
           "age_Ma", "Age_uncertainty_pos_Ma", "Age_uncertainty_neg_Ma",
           "proxy", "locality")]
  mod = data.frame(280, 5, 5, 0, 0.001, 0.001, "Instrumental", "Keeling")
  names(mod) = names(d)
  d = rbind(d, mod)

  #Parse data - co2 mean and uncertainty
  pco2 = log(d$CO2_ppm)
  ##Max and min in log deviations
  pco2.mm = log(d$CO2_ppm + d$CO2_uncertainty_pos_ppm/2) - pco2
  pco2.mm = cbind(pco2.mm, pco2 - log(d$CO2_ppm - d$CO2_uncertainty__neg_ppm/2))
  ##Average 1sd in log units
  pco2.sd = apply(pco2.mm, 1, mean)
  pco2.pre = 1 / pco2.sd^2
  
  #Parse data - ages and uncertainty
  pco2.age = d$age_Ma
  pco2.age.sd = apply(cbind(d$Age_uncertainty_pos_Ma, d$Age_uncertainty_neg_Ma), 1, mean)
  pco2.loc = d$locality
  pco2.prox = d$proxy
  
  dat = data.frame(pco2, pco2.sd, pco2.pre, pco2.age, 
                   pco2.age.sd, pco2.prox, pco2.loc)
  dat = dat[order(dat$pco2.loc),]

  return(dat)
}