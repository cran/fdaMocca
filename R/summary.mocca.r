

summary.mocca <- function (object,...) 
## summary method for mocca object...
{
  K <- object$initials$K ## number of clusters
  ## table for model parameters...
  p <- numOfEst_pars(object) ## total number of the estimated parameters
  tab_numOfCurves_cluster <- t(numOfCurves_cluster(object)) ## number of years/curves per cluster

  ## probs of cluster membership as a table...
  t.prob <- matrix(object$pars$probs,1,K)
  dimnames(t.prob) <- list("probs",paste("cluster:",1:K,sep=""))
  aic.bic.entropy <- criteria.mocca(object) 
  dimnames(aic.bic.entropy)[[1]] <- ""
  crita <- cbind(object$loglik,aic.bic.entropy)
  colnames(crita)[1] <- "Loglik"
     
  ret <- list(sig2=object$sig2, K=K, r=object$initials$r,N=object$initials$N, 
         p=p, tab_numOfCurves_cluster=tab_numOfCurves_cluster, covariates_est=object$pars$vk, t.probs=t.prob,
         crita=crita,moc=object$initials$moc)

  class(ret)<-"summary.mocca"
  ret
} ## end summary.mocca


print.summary.mocca <- function(x, digits = max(3, getOption("digits") - 3), ...)
## print method for mocca summary method...
{  cat("\n")  
   if (x$moc){ ## model with covariates
     cat("Model: with scalar covariates; N = ", x$N, sep = "")
   } else cat("Model: without scalar covariates; N = ", x$N, sep = "")
   cat("\n")
   cat("\n Number of clusters: ", x$K ,sep = "")
   cat("\n Number of estimated parameters: ", x$p ,sep = "")
   cat("\n")
   cat("\n Probabilities of cluster membership:", sep="")
   cat("\n")
   print(x$t.prob)
   cat("\n Number of objects in each cluster:", sep="")
   cat("\n")
   print(x$tab_numOfCurves_cluster)
   #cat("\n")
   cat("\n Residual variance est. of functional data = ", round(x$sig2[1], digits=3), sep="")

   if (x$moc){ ## model with covariates
     cat("\n")
     cat("\n ------------------------------ ")
     cat("\n Number of covariates: ", x$r ,sep = "")
     cat("\n")
     cat("\n Mean value estimates for scalar covariates:\n", sep = "")
   #  cat("\n")
     print(x$covariates_est)
     cat("\n Residual variance est. of scalar covariates = ", round(x$sig2[2], digits=3), sep="")
     cat("\n ------------------------------ ")
     cat("\n")
   }
   
   cat("\n")
   print(x$crita)
   
   cat("\n")
   invisible(x)
} ## end print.summary.mocca




