
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
