load("out/post.rda")

R2jags::traceplot(p, varname = c("pco2_m"))
R2jags::traceplot(p, varname = c("pco2_m.pre", "pco2_m.eps.ac"))
dev.off()

View(p$BUGSoutput$summary)

