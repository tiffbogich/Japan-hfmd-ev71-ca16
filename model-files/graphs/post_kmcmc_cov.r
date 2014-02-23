require(plom)

setwd('..')

models <- dir(path = ".", pattern = "^hfmd_*")

for(m in models){
  print(m)
  root <- file.path(m, "model")
  s <- get.settings(root=root, theta='theta_kmcmc_cov.json')

  for(p in 0:4){
    x = read.csv(file.path(root, 'results', 'kmcmc_cov', paste('acc_', p, '.csv', sep='')), header=TRUE)
    if(nrow(x)>0){
  
    png(file=file.path('graphs', 'kmcmc_cov', paste('post_', m, '_', p, '.png', sep='')), width=1000, height=1000)
    plot.posteriors(s, res='results/kmcmc_cov', id=p)
    dev.off()

    png(file=file.path('graphs', 'kmcmc_cov', paste('trace_', m, '_', p, '.png', sep='')), width=1000, height=1000)
    plot.trace(s, res='results/kmcmc_cov', id=p)
    dev.off()
    }
  }

}
