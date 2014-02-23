require(plom)

setwd('..')

models <- dir(path = ".", pattern = "^hfmd_*")

for(m in models){
  print(m)
  root <- file.path(m, "model")
  s <- get.settings(root=root, theta='theta_pmcmc_ode_cov.json')

  for(p in 0:4){
    x = read.csv(file.path(root, 'results', 'pmcmc_ode_cov', paste('acc_', p, '.csv', sep='')), header=TRUE)
    if(nrow(x)>0){
  
    png(file=file.path('graphs', 'pmcmc_ode_cov', paste('post_', m, '_', p, '.png', sep='')), width=1000, height=1000)
    plot.posteriors(s, res='results/pmcmc_ode_cov', id=p)
    dev.off()

    png(file=file.path('graphs', 'pmcmc_ode_cov', paste('trace_', m, '_', p, '.png', sep='')), width=1000, height=1000)
    plot.trace(s, res='results/pmcmc_ode_cov', id=p)
    dev.off()
    }
  }

}
