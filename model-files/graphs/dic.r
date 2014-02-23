setwd('..')

models <- dir(path = ".", pattern = "^hfmd_*")

dic <- list()

for(m in models){
  print(m)

  dic[[m]] <- NULL

  for(p in 0:2){
    res <- read.csv(file.path(m, "model", 'results', 'kmcmc_cov', paste('trace_', p, '.csv', sep='')))   

    llike <- res[,ncol(res)]
    dev<- -2*llike    
    dic[[m]] <- c(dic[[m]], 2*mean(dev)-min(dev))    
  }

}

