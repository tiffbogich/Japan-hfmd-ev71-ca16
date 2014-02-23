require(plom)

setwd('..')

models <- dir(path = ".", pattern = "^hfmd_*")
m <- models[1]

for(m in models){
  print(m)
  root <- file.path(m, "model")
  s <- get.settings(root=root, theta="theta_lhs_simplex.json")
  
  for(p in c(75, 85, 95)){
    png(file=file.path('graphs', paste('prs_', m, '_', p, '.png', sep='')), width=1000, height=1000)
    plot.prs(s, res='lhs_simplex', perc=p/100, par_sv=TRUE, par_obs=TRUE)
    dev.off()

    png(file=file.path('graphs', paste('endpoints_', m, '_', p, '.png', sep='')), width=1000, height=1000)
    plot.endpoints(s, res='lhs_simplex', perc=p/100, box=TRUE, par_sv=TRUE, par_obs=TRUE)
    dev.off()
  }
  
  png(file=file.path('graphs', paste('traj_lhs_simplex_', m, '.png', sep='')), width=700)
  plot.X(s)
  title(m)
  dev.off()
  
  png(file=file.path('graphs', paste('sv_lhs_simplex_', m, '.png', sep='')), width=700)
  plot.SV(s)
  title(m)
  dev.off()
  
}
