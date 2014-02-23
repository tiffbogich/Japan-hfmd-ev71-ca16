###########################
##LHS
###########################


#' Plot the results of a Latin Hypercube Sampling design (LHS)
#'
#' @param s a settings object (obtained with \code{\link{get.settings}})
#' @param res the directory name (not a path) where the results of the lhs design can be found
#' @param navg integer Instead of plotting the log likelihood value of the last iteration, an average of the back last points will be used
#' @param par_sv logical Should the initial conditions be plotted
#' @param par_sv logical Should the parameter of the porcess model be plotted
#' @param par_sv logical Should the parameter of the observation process be plotted
#' @param perc proportion keep only points whose likelihood is above the specified percentile
#' @export
plot.prs <- function(s, res='lhs', navg=NA, perc=0.0, par_sv=FALSE, par_proc=TRUE, par_obs=FALSE){
  
  path.design <- file.path(s$root,'results', res, 'design.csv')

  if (file.exists(path.design)) {

    s <- add_theta(s, file.path(s$root, 'results', res, 'theta.json'))
    des <- read.csv(path.design, header=TRUE)
    resbest <- NULL
    for(i in 1:nrow(des)) {
      mypath <- file.path(s$root, 'results', res, paste("best_", i-1, ".csv", sep=''))

      best <- tryCatch(read.csv(mypath, header=TRUE),
                       error=function(e){
                         nabest <- t(as.matrix(rep(NA, ncol(des)+1)))
                         colnames(nabest) <- c(names(des), 'log_like')
                         nabest
                       })

      if(!is.na(navg) & (navg <= nrow(best))){
          best <- apply(best[(nrow(best)-navg):nrow(best), 2:ncol(best)], 2, mean, na.rm=TRUE)
      } else {
          best <- best[nrow(best), 2:ncol(best)]
      }
      resbest <- rbind(resbest, best)      
    }

    ##trash what we don't want par_sv par_proc or par_obs AND sd_transf > 0.0
    keep <- NULL
    offset <- 1
    partype <- list('par_sv' = par_sv, 'par_proc' = par_proc, 'par_obs' = par_obs)
    for(ptype in names(partype)){
      for(par in s$settings$orders[[ ptype ]]){      
        for(group in s$settings$partition[[ s$settings$parameters[[par]][['partition_id']] ]][['group']]){
          if(partype[[ptype]] & (s$settings$parameters[[par]][['sd_transf']][[ group$id ]] > 0.0) ){
            keep <- c(keep, offset)
          }
          offset <- offset+1
        }
      }
    }

    like <- c(rep(NA, nrow(des)), as.numeric(resbest[,ncol(resbest)]))  
    resbest <- resbest[,-ncol(resbest)]
    full <- as.data.frame(rbind(as.matrix(des[,2:ncol(des)]),as.matrix(resbest)))
    full <- full[,keep]

    full$iota_1.all[full$iota_1.all < -3] <- -4
    full$iota_2.all[full$iota_2.all < -3] <- -4

    
    ##trash quantile specified by auto.like 
    qt <- quantile(like, perc, na.rm=TRUE)
    keep <- (like>=qt | is.na(like))
    full <- full[keep,]
    like <- like[keep]
    
    ##color (below media = blue, above=red)
    cols <- rep(NA, length(like))
    keep <- ((like>=median(like, na.rm=TRUE)) & !is.na(like))
    alpha <- rescale(like[keep], min(like[keep], na.rm=TRUE), max(like[keep], na.rm=TRUE), 0, 1)
    cols[keep] <- rgb(1,0,0, alpha=alpha)

    keep <- ((like<median(like, na.rm=TRUE)) & !is.na(like))
    alpha <- rescale(like[keep], min(like[keep], na.rm=TRUE), max(like[keep], na.rm=TRUE), 0, 0.4)
    cols[keep] <- rgb(0,0,1, alpha=alpha)
    
    pairs(full, pch=21, bg=cols, col=rgb(0.1,0.1,0.1,alpha=0.05))
    
  } else {
    print (paste("could not find the design:", path.design))
  }
}



#' Plot the results of a Latin Hypercube Sampling design (LHS)
#'
#' @param s a settings object (obtained with \code{\link{get.settings}})
#' @param res the directory name (not a path) where the results of the lhs design can be found
#' @param navg integer Instead of plotting the log likelihood value of the last iteration, an average of the back last points will be used
#' @param par_sv logical Should the initial conditions be plotted
#' @param par_sv logical Should the parameter of the porcess model be plotted
#' @param par_sv logical Should the parameter of the observation process be plotted
#' @param perc proportion keep only points whose likelihood is above the specified percentile
#' @param box logical keep only points whose value are within the min and max values of the design (extended to encompass the MLE)
#' @export
plot.endpoints <- function(s, res='lhs', navg=NA, perc=0.0, box=FALSE, par_sv=FALSE, par_proc=TRUE, par_obs=FALSE){
  
  path.design <- file.path(s$root,'results', res, 'design.csv')
  
  if (file.exists(path.design)) {
    
    s <- add_theta(s, file.path(s$root, 'results', res, 'theta.json'))
    des <- read.csv(path.design, header=TRUE)
    resbest <- NULL
    for(i in 1:nrow(des)) {
      mypath <- file.path(s$root, 'results', res, paste("best_", i-1, ".csv", sep=''))

      best <- tryCatch(read.csv(mypath, header=TRUE),
                       error=function(e){
                         nabest <- t(as.matrix(rep(NA, ncol(des)+1)))
                         colnames(nabest) <- c(names(des), 'log_like')
                         nabest
                       })

      if(!is.na(navg) & (navg <= nrow(best))){
        best <- apply(best[(nrow(best)-navg):nrow(best), 2:ncol(best)], 2, mean, na.rm=TRUE)
      } else {
        best <- best[nrow(best), 2:ncol(best)]
      }
      resbest <- rbind(resbest, best)      
    }
    
    ##trash what we don't want par_sv par_proc or par_obs AND sd_transf > 0.0
    keep <- NULL
    offset <- 1
    partype <- list('par_sv' = par_sv, 'par_proc' = par_proc, 'par_obs' = par_obs)
    for(ptype in names(partype)){
      for(par in s$settings$orders[[ ptype ]]){      
        for(group in s$settings$partition[[ s$settings$parameters[[par]][['partition_id']] ]][['group']]){
          if(partype[[ptype]] & (s$settings$parameters[[par]][['sd_transf']][[ group$id ]] > 0.0) ){
            keep <- c(keep, offset)
          }
          offset <- offset+1
        }
      }
    }
    
    like <- resbest[,ncol(resbest)]
    resbest <- resbest[,keep]
    des <- des[,keep+1]

    resbest$iota_1.all[resbest$iota_1.all < -3] <- -4
    resbest$iota_2.all[resbest$iota_2.all < -3] <- -4

    
    des.min <- apply(des, 2, min, na.rm=TRUE)
    des.max <- apply(des, 2, max, na.rm=TRUE)
    
    ind.mle <- which(like == max(like, na.rm=TRUE))
   
    myylim <- c(quantile(like, perc, na.rm=TRUE), max(like, na.rm=TRUE))
    
    par(mar=c(3,3,0.5,0.5))
    n <- length(keep)
    layout(matrix(1:(ceiling(sqrt(n))*round(sqrt(n))), round(sqrt(n)), ceiling(sqrt(n))))

    for(offset in 1:length(keep)){

      x <- resbest[, offset]
      xkeep <- x[like>=myylim[1]]

      if(box){
        myxlim <- c(min(des.min[offset], resbest[ind.mle, offset], na.rm=TRUE), max(des.max[offset], resbest[ind.mle, offset], na.rm=TRUE))
      }else{
        myxlim <- c(min(des.min[offset], xkeep[xkeep > -Inf], na.rm=TRUE), max(des.max[offset], xkeep[xkeep < Inf], na.rm=TRUE))
      }
      
      plot(resbest[,offset], like,
           xlim=myxlim,
           ylim=myylim,
           xlab=names(resbest)[offset],
           ylab="log likelihood",
           col='grey',
           pch='.',
           mgp=c(2,1,0))
      points(resbest[ind.mle, offset], like[ind.mle], col='red')

      abline(v=des.min[offset], lty=2, col='grey')
      abline(v=des.max[offset], lty=2, col='grey')
      text(resbest[, offset], like, 1:nrow(resbest)-1, cex=0.5)                 
    }
    
    for(m in ind.mle){
      print(paste("best index:", m-1, ": ", like[m],  sep=" "), quote=FALSE)
      print(resbest[m,], row.names=FALSE)
    }
    
  } else {
    print (paste("could not find the design:", path.design))
  }  
}



