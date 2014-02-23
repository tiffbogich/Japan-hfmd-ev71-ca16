require(plom)

setwd('~/Dropbox/hfmd_3/vacc')

models <- dir(path = ".", pattern = "*_v")

for(m in models){

  print(m)
  root <- file.path(m, 'model', 'results', 'vaccination', 'p', 'all')

  design <- read.csv(file.path(root,'design.csv'), header=TRUE)  
  p <- design$p.all
  H <- length(p)

  
  inf.ev <- rep(NA, H)
  inf.ca <- rep(NA, H)
  
  for(i in 1:H){
    path <- file.path(root, paste('X_', i-1, '.csv', sep=''))
    if(file.exists(path)){
      X <- read.csv(path, header=TRUE)
      x <- apply(X,2,mean)
    
      inf.ev[i] <- x['obs_mean.EV71_JAP__IASR__inc']
      inf.ca[i] <- x['obs_mean.CA16_JAP__IASR__inc']
    } else {
      inf.ev[i] <- NA
      inf.ca[i] <- NA
    }

  }

  pdf(file=paste('vacc_', m, '.pdf', sep=''))

  layout(matrix(1:2, 1, 2))
  plot(p, inf.ev, ylab="average observed weekly incidence EV71")
  title(m)
  plot(p, inf.ca, ylab="average observed weekly incidence CA16")  
  title(m)
  dev.off()
  
}


##design <- read.csv('design.csv', header=TRUE)
##
##p <- design$p.all
##H <- length(p)
##
##inf.ev <- rep(NA, H)
##reinf.ev <- rep(NA, H)
##inf.ca <- rep(NA, H)
##reinf.ca <- rep(NA, H)
##
##i <- H-1
##X <- read.csv(paste('X_', i-1, '.csv', sep=''), header=TRUE)
##plot(X$obs_mean.EV16_inf__theo__inc, type='l')
##
##for(i in 1:H){
##  X <- read.csv(paste('X_', i-1, '.csv', sep=''), header=TRUE)
##  x <- apply(X,2,mean)
##
##  inf.ev[i] <- x['obs_mean.EV71_inf__theo__inc']
##  reinf.ev[i] <- x['obs_mean.EV71_reinf__theo__inc']
##  inf.ca[i] <- x['obs_mean.CA16_inf__theo__inc']
##  reinf.ca[i] <- x['obs_mean.CA16_reinf__theo__inc'] 
##}
##
##
##pdf(file='vaccination.pdf')
##layout(matrix(1:2, 1, 2))
##plot(p, inf.ev, ylim=c(min(inf.ev, reinf.ev), max(inf.ev, reinf.ev)), ylab="average observed weekly incidence EV71")
##points(p, reinf.ev, col='red')
##legend('topright', c('infection', 're-infection'), lty=1, col=c('black', 'red'))
##
##plot(p, inf.ca, ylim=c(min(inf.ca, reinf.ca), max(inf.ca, reinf.ca)), ylab="average observed weekly incidence CA16")
##points(p, reinf.ca, col='red')
##
##dev.off()
