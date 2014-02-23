require(plom)

setwd('..')

models <- dir(path = ".", pattern = "^hfmd_*")
m <- models[1]

for(m in models){
  print(m)
  root <- file.path(m, "model")
  s <- get.settings(root=root, theta="mle.json")
    
  png(file=file.path('graphs', paste('traj_ksimplex_', m, '.png', sep='')), width=700)
  plot.X(s)
  title(m)
  dev.off()
  
  png(file=file.path('graphs', paste('sv_ksimplex_', m, '.png', sep='')), width=700)
  plot.SV(s)
  title(m)
  dev.off()
  
}
