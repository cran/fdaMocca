###############################
## printing mocca() results ...
###############################

print.mocca <- function (x,...) 
## default print function for mocca objects...
{
   # cat("\n")
   if (x$initials$moc){ ## model with covariates
     cat("Model: with scalar covariates\n", sep = "")
     cat("\n Number of covariates: ", x$initials$r ,sep = "")
   } else cat("Model: without scalar covariates\n", sep = "")
     
  # cat("\n")
   cat("\n Number of clusters: ", x$initials$K ,sep = "")
  cat("\n")
   cat("\n Residual variance est. for functional data = ", round(x$sig2[1], digits=3), sep="")
   if (x$initials$moc){ ## model with covariates
     cat("\n Residual variance est. for scalar covariates = ", round(x$sig2[2], digits=3), sep="")
   }
   cat("\n")
   cat("\n Loglik = ", round(x$loglik, digits=3), ";  N = ", x$initials$N, "\n",sep="")
   
  cat("\n")
  invisible(x)
}


