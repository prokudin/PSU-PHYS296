"This script file will plot pT as a function of ALL with error stat_u and sys_u"

# Define the function ----
Rplot <- function(pT, value, stat_u, sys_u){
  # Pulling Error
  error = sqrt(stat_u^2+sys_u^2)
  
  # Ploting
  plot(pT, value,
       ylim=range(c(-1*delta/2, delta/2)),
       pch=19, xlab="pT", ylab="ALL",
       main="Scatter plot of pT as f(ALL)"
  )
  # add error bars
  arrows(pT, -1*delta/2, pT, delta/2, length=0.05, angle=90, code=3)
}

# Testing ----

library(xlsx)
data = 