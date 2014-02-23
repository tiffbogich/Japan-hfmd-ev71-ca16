require(plom)
source('func.r')

setwd('~/Dropbox/hfmd_3/results')

models <- dir(path = ".", pattern = "^hfmd_*")

for(m in models){
  print(m)
  root <- file.path(m, "model")
  s <- get.settings(root=root)

#  plot.endpoints(s, res='lhs_simplex', perc=95/100, box=TRUE)
  
  for(p in c(75, 85, 95)){
    png(file=paste('prs_', m, '_', p, '.png', sep=''), width=1000, height=1000)
    plot.prs(s, res='lhs_simplex', perc=p/100)
    dev.off()

    png(file=paste('endpoints_', m, '_', p, '.png', sep=''), width=1000, height=1000)
    plot.endpoints(s, res='lhs_simplex', perc=p/100, box=TRUE)
    dev.off()
  }
  
  png(file=paste('traj_', m, '.png', sep=''), width=700)
  plot.X(s)
  title(m)
  dev.off()
}
