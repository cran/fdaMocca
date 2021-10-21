
## function to calculate number of the estimated parameters to be used with loglik, entropy, AIC, BIC ...
numOfEst_pars <- function(x){ ## lista=svar6cl4cov)
 ## function to sum up the number of estimated parameters in the model...
 l0 <- length(x$pars$lambda.zero)
 La <- length(x$pars$Lambda)
 al <- length(x$pars$alpha)
 pi <- length(x$pars$probs)
 si <- 1 ## length(x$pars$sig2)
 if(x$initials$moc) { ## model with covariates
    d <- dim(x$pars$Delta)
    lenG <- d[1]*(d[2]+1)*d[3]/2
    vk <- length(x$pars$vk)
    si <- si +1
 } else { ## without covariates
     d <- dim(x$pars$Gamma)
     lenG <- d[1]*(d[2]+1)*d[3]/2 ## dim(x$pars$Gamma)[1]*(1+dim(x$pars$Gamma)[2])/2
     vk <- 0
   }
 num <- l0+La+al+pi+si+lenG+vk
 num
}

## function to calculate the model entropy, AIC and BIC...
entropy.function <- function(vek) {
     good <- vek >0
     vek[!good] <- 1e-12 ## set the very small value if probability is zero
    -sum(vek*log(vek))
}

criteria.mocca <- function(x){
 ## function to calculate entropy (the Shannon entropy), AIC, and BIC criteria of the fitted mocca model...
  npar <- numOfEst_pars(x)
  if(x$initials$moc) 
        probabilities <- x$vars$tag.piigivej
    else probabilities <- x$vars$piigivej
  n <- length(unique(x$data$curve))
  logl <- x$loglik
  AIC <- 2*npar - 2*logl
  BIC <- npar*log(n) - 2*logl
  entropy <- entropy.function(probabilities)
  out <- as.matrix(c(AIC=AIC,BIC=BIC,entropy=entropy/n),ncol=1)
  dimnames(out)[[2]] <- deparse(substitute(x))
  t(out)
}

## function to calculate number of subjects/curves in each cluster...
numOfCurves_cluster <- function(x){
 ## hard coded cluster belonging. The year goes to the cluster with the highest probability...
 cluster.prob <- x$vars$tag.piigivej
 if(is.null(cluster.prob)) cluster.prob <- x$vars$piigivej
 num.of.clusters <- dim(cluster.prob)[2]
 cluster <- apply(cluster.prob,1,order)[num.of.clusters,]
 ttab <- t(table(cluster))
 ttab <- unname(t(ttab))
 dimnames(ttab)[[1]]<-paste("cluster:",1:nrow(ttab),sep="")
 dimnames(ttab)[[2]]<-"#curves"
 ttab
}


###############################################
## logLik for a fitted mocca as S3 method ...
##############################################
logLik.mocca <- function(object,...){
  ## extracting loglik from the output of the fitted mocca, rather than re-calculating it as in the version v0 below...
   ll <- object$loglik
   p <- numOfEst_pars(object)
   attr(ll, "df") <- p 
   class(ll) <- "logLik"
   ll
} ## logLik.mocca


logLik_v0.mocca <- function(object,...){
  data <- object$data
  parameters <- object$pars
  design <- object$design
  N<-length(unique(data$curve))
  Lambda.alpha <- as.vector(parameters$lambda.zero) + parameters$Lambda %*% t(parameters$alpha)

  nc <- parallel::detectCores()
  cl <- parallel:: makeCluster(nc-1)
  doParallel:: registerDoParallel(cl)
 
  ## First, if we have a single Gamma-matrix...
  if(is.null(parameters$Gamma) & is.null(parameters$Delta)) {
    K<-length(parameters$probs)
    parameters$Gamma<-array(0,c(K,length(parameters$lambda.zero),length(parameters$lambda.zero)))
    for(i in 1:K) parameters$Gamma[i,,]<-parameters$Gamma
  }

  ## Next, if we have separate Gamma-matrices and covariates...
  j <- NULL ## binding the variable locally to the function, so the 'R CMD check' has nothing to complain about 
  if(!is.null(design$tag.S)) {
    ## print('hello')
    K <- length(parameters$probs)
    tag.S<-design$tag.S
    antcov<-dim(data$covariates)[2]
    tag.GammaK <- parameters$Delta
    tag.Lambda.alpha <- rbind(Lambda.alpha,t(parameters$vk))
    summa <- foreach::foreach(j = 1:N, .combine = cbind, .packages = "mvtnorm") %dopar% {
         tag.Sj <- tag.S[data$tag.curve == j,  ]
         nj <- sum(data$tag.curve == j)
         varm <- diag(parameters$sig2, nj)
         nj<-nj+1
         for(sista in 1:antcov){
               nj<-nj-1
               varm[nj,nj]<-parameters$sig2x ## matrix R in the paper
         }
         centx <- data$tag.x[data$tag.curve == j] - tag.Sj %*% tag.Lambda.alpha
         d <- NULL
         for(k in 1:K){
             covx <- tag.Sj %*% tag.GammaK[k,,] %*% t(tag.Sj) + varm ## covariance matrix of the model with covariates (ui|zi)
             d <- c(d,dmvnorm(centx[,k],rep(0,length(centx[,k])),covx) * parameters$probs[k])
         }
         d
      } ## end foreach
  } else{
     K <- length(parameters$probs)
     S<-design$S
     GammaK <- parameters$Gamma
    ## j <- NULL ## binding the variable locally to the function, so the 'R CMD check' has nothing to complain about 
     summa <- foreach:: foreach(j = 1:N, .combine = cbind, .packages = "mvtnorm") %dopar% {
                Sj <- S[data$curve == j,  ]
                nj <- sum(data$curve == j)
                varm <- diag(parameters$sig2, nj)
                centx <- data$x[data$curve == j] - Sj %*% Lambda.alpha
                d <- NULL
                for(k in 1:K){
                       covx <- Sj %*% GammaK[k,,] %*% t(Sj) + varm ## covariance matrix of the model without covariates
                       d <- c(d,dmvnorm(centx[,k],rep(0,length(centx[,k])),covx) * parameters$probs[k])
                }
                d
             } ## end foreach
   } ## end else

   if (!is.null(cl)) parallel::stopCluster(cl) 

   summa<-apply(summa,2,sum)
   summa<- log(summa)  ## log(summa+1e-5)
   summa<-sum(summa)
   p <- numOfEst_pars(object)
   attr(summa, "df") <- p 
   class(summa) <- "logLik"
   summa
} ## logLik_v0.mocca





