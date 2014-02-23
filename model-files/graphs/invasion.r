models <- c('siri_hbrs', 'siqri_hbrs', 'siqri_hbrs_b')

for(m in models){
  design <- read.csv(file.path(m, 'model', 'results', 'invasion', 'z', 'all', 'design.csv'), header=TRUE)

  sigma <- seq(0, 1, length=nrow(design))
  fit <- matrix(NA, nrow(design), length(sigma))
  pop <- 126708453
  z <- design$z.all
  
  for(i in 1:nrow(design)){
    for(j in 1:length(sigma)){
      X <- read.csv(file.path(m, 'model', 'results', 'invasion', 'z', 'all', paste('X_', i-1, '.csv', sep='')), header=TRUE)
      X <- apply(X, 2, mean)
        
      
      r0EV <- design$r0_1.all[i]
      v <- (1/(design$v.all[i]/7.0))
      mu_d <- 0.000181634909051118
      fit[i, j] <- r0EV * v / pop * (X['SS.Japan__all'] + sigma[j]*X['SR.Japan__all']) -v - mu_d
      fit[i, j] <- as.numeric(fit[i, j] > 0)      
    }
  }

  png(file=paste('invasion_', m, '.png', sep=''))
  image(z, sigma,  fit, col=c('white', 'black'))
  abline(0,1, col='red')
  title(m)  
  dev.off()
}


require(plom)
s=get.settings(root='siqri_hbrs/model')
plot.SV(s, res='results/invasion/z/all', id=1)
plot.X(s, res='results/invasion/z/all', o=F, id=19)
