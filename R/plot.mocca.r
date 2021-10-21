## plot functions ...

##########################################################################
## function to illustrate the covariance structure within each cluster....
##########################################################################
covariance.illustration <- function(x,transform=FALSE,covariance=TRUE,covariates =FALSE){ ## splines=TRUE
## routine to illustrate the covariance structure within each cluster...
## covariance' - whether covariance or correlation matrices should be plotted,
## covariates' - whether covariates should be included
## transform' - whether svd back-transformation of the spline model matrix should be applied
  parameters <- x$pars
  vars <- x$vars
  Design <- x$design
  K <- x$initials$K
  d.splines<- x$initials$q ## number of spline coefficients
  moc <- x$initials$moc ## whether model with covariates (TRUE)
 ## r <- ceiling(sqrt(K))
  
  if (moc){ ## model with covariates
     d.cov<- x$initials$r # number of covariates
     no.bspl<-x$initials$Q # number of spline coeffs and covariates
  } else{
      d.cov <- 0
      no.bspl<- d.splines # number of spline coeffs
    }

  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))  

  ## if(K==5 | K==6) par(mfrow=c(2,3)) else par(mfrow=c(r,r))

  if (!moc){ ## model without covariates
     if(covariates) ## covariates should be included
        warning("model without covariates, so covariance matrices include only spline coefficients")
     covmat<-parameters$Gamma
  } else{ ## model with covariates
          if(!covariates) ## covariates should not be included
             covmat<-parameters$Delta[,1:d.splines,1:d.splines] 
          else covmat<-parameters$Delta
    }
  
 
  if(transform) {
    svdFullS<-svd(Design$FullS.bmat)
    #Transf<-diag(svdFullS$d)%*%t(svdFullS$v)
    BackTransf<-svdFullS$v%*%diag(1/svdFullS$d)
    #round(BackTransf%*%parameters1$lambda.zero,2)
    #round(BackTransf%*%parameters1$Lambda,2)
    #round(parameters1$alpha,2)
    #round(parameters1$sigma,2)
    #round(parameters1$pi,2)
    if(!covariates) for(k in 1:K) covmat[k,,]<-BackTransf%*%t(covmat[k,,])%*%t(BackTransf)
     else{
            for(k in 1:K) {
             covmat[k,1:d.splines,1:d.splines]<-BackTransf%*%covmat[k,1:d.splines,1:d.splines]%*%t(BackTransf)
             temp <- BackTransf%*%t(covmat[k,c(d.splines+1):no.bspl,1:d.splines])
             covmat[k,c(d.splines+1):no.bspl,1:d.splines]<-t(temp)
             covmat[k,1:d.splines,c(d.splines+1):no.bspl] <- temp
            }
          }
  } ## end if transform
  if(!covariance){
       for(k in 1:K) covmat[k,,] <- 100*cov2cor(covmat[k,,])
  } 
  covmat1<-covmat
  #else 
  if(covariance){
         covmat<-covmat-min(covmat)
         # covmat<-100*max(covmat)*covmat
  }
  range.all<-range(round(covmat))
  range.all[1]<-range.all[1]-1
  range.all[2]<-range.all[2]+1
  tal<-diff(range.all)
  if(!covariates) no.bspl <- d.splines
  xval<-yval<-seq(0,1,len=no.bspl)
  yval<-rev(yval)
  farg<-wtcolors(tal)
  #par(mfrow=c(2,2))
  if(length(dim(covmat))==3){
       for(i in 1:K){ 
           minG<-min(round(covmat[i,no.bspl:1,]))
           maxG<-max(round(covmat[i,no.bspl:1,]))
           farger<-farg[match(minG:maxG,range.all[1]:range.all[2])]
           image(t(covmat[i,no.bspl:1,]),col=farger,axes=F)
           axis(1,at=seq(0,1,len=no.bspl),1:no.bspl)
           axis(2,at=seq(0,1,len=no.bspl),no.bspl:1,las=1)
           if(covariance) title(paste('Covariance of cluster',i)) else title(paste('Correlation of cluster',i))
           for(j in 1:no.bspl){
               for(k in 1:no.bspl)
                   text(xval[k],yval[j],round(covmat1[i,k,j],2),cex=0.8)
           }
        }
  } else{
        par(mfrow=c(1,1))
        image(t(covmat[no.bspl:1,]),col=wtcolors(tal),axes=F)
        axis(3,at=seq(0,1,len=no.bspl),1:no.bspl)
        axis(2,at=seq(0,1,len=no.bspl),no.bspl:1,las=1)
        xval<-yval<-seq(0,1,len=no.bspl)
        yval<-rev(yval)
        for(j in 1:no.bspl){
           for(k in 1:no.bspl)
            text(xval[k],yval[j],round(covmat[k,j]))
        }
      }
} ## end covariance.illustration()

####################################################################
## function to plot all cluster mean curves in one plot, 
## also to print loglik,  est. variances, probs...
####################################################################

illustrations.plot <- function(x,select =NULL,lwd=2,ylab="",ylim=NULL,ncolors=NULL){ # sekvens=NULL
 ## this is a default plot for plot.mocca() method
 ## that produces cluster means curves in one plot...
  par<-x$pars
  var<-x$vars
  data<-x$data
  design<-x$design
  iter.hist<-x$score.hist
  iter<- x$iter-1
  if (iter==0) stop("EM algorithm terminated at the 1st iteration, no plots produced")
  K <- x$initials$K
  if(is.null(select)) 
          select <- 1:K
  if(is.null(ncolors)) 
         ncolors <- K
    
  farg <- wtcolors(ncolors)
  ## logl<-round(iter.hist[iter+1,'ll'])
  ll <- round(iter.hist[iter+1,'ll'],5) ## max loglik
  ll.old<-iter.hist[iter,'ll']
  
  cluster.mean <- matrix(0,K,dim(design$FullS)[1])
  for(k in 1:K) 
          cluster.mean[k,] <- design$FullS %*% (par$lambda.zero + par$Lambda %*% par$alpha[k,  ])
  if (is.null(ylim)) ylim <- c(min(cluster.mean)-.1*min(cluster.mean), max(cluster.mean)+.1*max(cluster.mean))     

  if(x$initials$moc) { ## model with covariates
        v1<-table(apply(var$tag.piigivej,1,order)[K,])
        frekvens<-round(par$probs[select],3)
       ## Sigmas<-paste(round(iter.hist[iter+1,1],3),'and',round(iter.hist[iter+1,2],3))
       ## sigma.old<-iter.hist[iter,1]+iter.hist[iter,2]
       ## sigma.ny<-iter.hist[iter+1,1]+iter.hist[iter+1,2]
        sig2 <- round(iter.hist[iter+1,1],3) ## est. sig^2
        sig2x <- round(iter.hist[iter+1,2],3) ## est. sig^2 of covariates, sig2_x
        plot(data$grid,data$grid,type='n',ylab="", ylim=ylim,
            xlab=bquote("ll" ==.(ll) ~ ", " ~ hat(sigma)^2 ==.(sig2) ~ "and" ~ hat(sigma)[x]^2 ==.(sig2x)) )
  }
  else { ## model without covariates
        v1<-table(apply(var$piigivej,1,order)[K,])
        frekvens<-round(par$probs[select],3)
      ##  Sigmas<-round(iter.hist[iter+1,1],3)
      ##  sigma.old<-iter.hist[iter,1]
      ##  sigma.ny<-iter.hist[iter+1,1]
        sig2 <- round(iter.hist[iter+1,1],3) ## current sig^2
        plot(data$grid,data$grid,type='n',ylab="",ylim=ylim,xlab=bquote("ll" ==.(ll) ~ "and" ~ hat(sigma)^2  ==.(sig2)) ) ## xlab=paste('ll =',logl,'\n sigma(s) =',Sigmas),ylab="Cluster Mean"
  }
  
  ## title(paste('EM-iteration no:',iter+1));
  title(paste('Cluster means'))
  grid()
  legend('bottomleft',legend=paste('q = ',length(par$lambda.zero),'; #est. pars = ',numOfEst_pars(x), sep=""))
  legend('topright',inset=0, legend=paste('cluster',select, ': pr',' = ',frekvens,', nk = ',v1[select], sep=""),lty=1,col=farg[select]) # #objects = ',v1[select], 
  ##legend('bottomleft',legend=paste('loglik rel diff =',round(abs(ll.old - ll)/(.1 + abs(ll)),5)))
  lines(data$grid,design$FullS%*%par$lambda.zero,lty=2,lwd=lwd+1) ## overall mean...
  #i<-0
  for(k in select){ ## plotting cluster means..
      # i<-i+1
       lines(data$grid,cluster.mean[k,], col = farg[k], lwd = lwd)
  }
} ## end illustrations.plot()

#########################################################################
## fuction to produce separate plots of the cluster mean curve for each cluster...
##########################################################################

trendplot <- function(x,probs=FALSE,pts=FALSE,size=50,select=NULL,years=NULL,years.names=NULL,ylim=NULL){ # j=0,years==c(-4486:1901) 
## the function produces the trend of the frequencies of the different clusters, together with mean probabilites (if probs=T)
## It also gives the mean value of the included covariates within each cluster (not the model estimated covariate values).
## If pts=T then it plots points of the frequency trend.
## size' is the bin size (here 50 years) of how many of those 50 years belong to a specific cluster.
## select' can be given if you want to plot the clusters in a different select or just a few of them.
## j' can be used to tell which cluster it starts with when you plot a few of the clusters,
## ex., if you have 6 clusters and you want to plot cluster 2:4 then you set select=2:4 and j=1.
## stand.cov' tells whether standardized covariates were used in the modelling.
## OrigData' is needed for calculating the number of years within a cluster bin (Used in the function nysummering).
## years' gives what years you want to calculate frequencies in. 
   parameters <- x$pars
   vars <- x$var
   data <- x$data
   design <- x$design
   if (!is.null(years.names))
              years.names <- years.names
   else if (!is.null(x$initials$years.names))
              years.names <- x$initials$years.names
   else years.names <- as.character(years)
  
   nytabell <- nysummering(vars=vars, years=years, size=size, years.names=years.names) ##, orig.data=orig.data)
   #summary(nytabell)
   svar <- plotClusterMean(x)
   farger<-rev(c('darkgreen','green','lightgreen','orange','indianred','brown','darkred','white','white','white'))
   #antal<-apply(nytabell$tabell[,2:dim(svar)[2]],2,sum)
   medelvarde<-apply(nytabell$tabell[,2:dim(svar)[2]],2,mean)
   K <- dim(svar)[2]-2
   ylim.up<-max(pretty(nytabell$tabell[,2:(K+1)])) + .1*max(pretty(nytabell$tabell[,2:(K+1)]))
   pic.low <- 0.67*ylim.up
   frac <- 1.06 ## 6804/6387
   max.years<-round(frac*(max(years)))
   min.years<-round(frac*(min(years)))
   omf <- max.years-min.years
   prop.low <- 0.4 ## 1476/6804
   prop.up <- 0.6 ## 3176/6804
   svar1 <- omskalare(svar[,2:dim(svar)[2]],pic.low,ylim.up)
   left <- min.years+prop.low*omf # -3200
   right <- min.years+prop.up*omf# -1500
   svar2 <- omskalare(svar[,1],left,right)
   svar3 <- cbind(svar2,svar1)
#if(K<4) par(mfrow=c(2,2))
#if(K==4 | K==5) par(mfrow=c(2,3))
#if(K==6 | K==7| K==8) par(mfrow=c(3,3))
#if(K>8) par(mfrow=c(3,4))
   ## if(stand.cov) data$covariates<- data$covariates
   if(!is.null(parameters$vk)){  ## model with covariates...
      sannolikheter.hard<-apply(vars$tag.piigivej,1,order)[K,]
      parameters$vk <- NULL
      cnames<-colnames(x$pars$vk) ## names of the covariates used
      covariates <- data$covariates[,cnames]
      for(k in 1:K) parameters$vk <-rbind(parameters$vk,apply(as.matrix(covariates[sannolikheter.hard==k,]),2,mean))
      ce <- round(t(parameters$vk),3)
      
      springs<-NULL
      if ("landmark" %in% colnames(covariates)==1) {## 'landmark' covariate exists so, starting date of Springs can be estimated.... 
         landmark.index <- which(colnames(covariates)=="landmark")
         no.of.days<-365.25*parameters$vk[,landmark.index]
         for(i in 1:K) {
             what.day<-round(no.of.days[i]-31-28.25-31-30)
             if(what.day>0) spring<-paste('Spring starts May',what.day)
              else{
                 what.day<-round(no.of.days[i]-31-28.25-31)
                 if(what.day<0) {
                     what.day<-round(no.of.days[i]-31-28.25)
                     spring<-paste('Spring starts Mars',what.day)
                 } else {
                       if(what.day==0) what.day<-1
                       spring<-paste('Spring starts April',what.day)
                   }
              }
             springs<-c(springs,spring)
       }
     }
   } else{ ## model without covariates
          sannolikheter.hard<-apply(vars$piigivej,1,order)[K,]
         
        ##   for(k in 1:K) parameters$vk <-rbind(parameters$vk,apply(TestCov[sannolikheter.hard==k,],2,mean)) ## ??? TestCov
        ##   dimnames(parameters$vk)[[2]]<-dimnames(TestCov)[[2]]
        ##   x$pars$vk <- parameters$vk
        ##   ce<-round(t(parameters$vk),3)
        ##   no.of.days<-round(365*parameters$vk[,4])
        ##  springs<-NULL
        ##  for(i in 1:K) {
        ##       what.day<-round(no.of.days[i]-31-28.25-31-30)
        ##       if(what.day>0) spring<-paste('Spring starts May',what.day)
        ##         else{
        ##            what.day<-round(no.of.days[i]-31-28.25-31)
        ##            if(what.day<0) {
        ##               what.day<-round(no.of.days[i]-31-28.25)
        ##               spring<-paste('Spring starts Mars',what.day)
        ##            } else {
        ##                if(what.day==0) what.day<-1
        ##                spring<-paste('Spring starts April',what.day)
        ##              }
        ##       }
        ##       springs<-c(springs,spring)
        ##   }
       }
   
   farg <- wtcolors(10)
   antal<- table(sannolikheter.hard)
   
   if (is.null(ylim)) ylim <- c(0,ylim.up)

   if(is.null(select)) select <- 1:K
   m<-0
   if (x$initials$moc){ ## model with covariates...
     for(i in select){
        m<-m+1
        plot(nytabell$tabell[,1],nytabell$tabell[,i+1],type='l',col=1,ylim=ylim,lwd=1,lty=1,xlab=paste('NoOfObs =',antal[i],springs[i]),ylab=paste('Bin size =',size))
        grid()
        if(pts) points(nytabell$tabell[,1],nytabell$tabell[,i+1],col=farger[nytabell$fargning[,i+1]],pch=15,cex=1)
        #abline(h=medelvarde[i])
        lines(svar3[,1],svar3[,2],lwd=2,lty=2)
        lines(svar3[,1],svar3[,i+2],col=farg[m],lwd=2)
        lines(c(left,left),c(pic.low-0.5,ylim.up+0.5))
        lines(c(right,right),c(pic.low-0.5,ylim.up+0.5))
        lines(c(left,right),c(pic.low-0.5,pic.low-0.5))
        lines(c(left,right),c(ylim.up+0.5,ylim.up+0.5))
        if(probs) 
           legend("topleft",c('0.3<p<0.4','0.4<p<0.5','0.5<p<0.6','0.6<p<0.7','0.7<p<0.8','0.8<p<0.9','0.9<p<1.0'), pch=15, col=farger[4:10], bty="n")
        legend('topright',inset=0, legend=paste(t(c(paste(dimnames(as.matrix(x$pars$vk))[[2]],'='))),t(ce[,i])),xjust=0, bty="n")
        title(paste('Cluster',i))
      }
   } else{ ## model without covariates...
          for(i in select){
             m<-m+1
             plot(nytabell$tabell[,1],nytabell$tabell[,i+1],type='l',col=1,ylim=c(0,ylim.up),lwd=1,lty=1,xlab=paste('NoOfObs =',antal[i]),ylab=paste('Bin size =',size))
             grid()
             if(pts) points(nytabell$tabell[,1],nytabell$tabell[,i+1],col=farger[nytabell$fargning[,i+1]],pch=15,cex=1)
             #abline(h=medelvarde[i])
             lines(svar3[,1],svar3[,2],lwd=2,lty=2)
             lines(svar3[,1],svar3[,i+2],col=farg[m],lwd=2)
             lines(c(left,left),c(pic.low-0.5,ylim.up+0.5))
             lines(c(right,right),c(pic.low-0.5,ylim.up+0.5))
             lines(c(left,right),c(pic.low-0.5,pic.low-0.5))
             lines(c(left,right),c(ylim.up+0.5,ylim.up+0.5))
             if(probs) 
                legend("topleft",c('0.3<p<0.4','0.4<p<0.5','0.5<p<0.6','0.6<p<0.7','0.7<p<0.8','0.8<p<0.9','0.9<p<1.0'), pch=15, col=farger[4:10], bty="n")
             title(paste('Cluster',i))
          }
     }  

} ## end trendplot()


#######################################
## several supporting functions...
#######################################
wtcolors <- function(ncolors){
    upside <- rainbow(ncolors, start = 0, end = 0.7)
    down <- upside
    for (i in 1:ncolors) {
        down[i] <- upside[ncolors + 1 - i]
    }
    down
}


omskalare<-function (x = cbind(1:5,6:10), nedre = 10, ovre = 20){
#
# Sw                En
# omskalare     re-scale(r)
# ovre (övre)       upper
# nedre             lower
#
## This rescale something and it use it for plotting in plots. Here it is used in trendplot.        
    new.var <- x/max(x, na.rm = T)
    new.var <- new.var - mean(new.var, na.rm = T)
    rnv <- range(new.var, na.rm = T)
    lim <- ovre - nedre
    new.var <- (ovre - (rnv * lim/(rnv[2] - rnv[1])))[2] + new.var * 
        lim/(rnv[2] - rnv[1])
    new.var
}

nysummering <- function(vars,years, size=50, years.names){ ## ,orig.data
## This function gives the hard coded cluster belongings in bins of 'size'.
## It also gives the mean probabilities within the bins and the color pre-selected for that.
## Output is a list.The result is used in the function trendplot.
   if(!is.null(dim(vars$tag.piigivej))) { ## model with covariates...
      dimension <- dim(vars$tag.piigivej)
      lista <- apply(vars$tag.piigivej,1,order)[dimension[2],]
      sht <- apply(vars$tag.piigivej,1,max)
   } else{
        dimension <- dim(vars$piigivej)
        lista<-apply(vars$piigivej,1,order)[dimension[2],]
        sht<-apply(vars$piigivej,1,max)
     }
   
  ## nl2<-dimnames(orig.data)[[1]]
  ## nl2 <- dimnames(varve$covariates)[[1]]
   nl2 <- years.names

   names(lista) <- nl2
   nl2 <- years[is.na(match(years,nl2))]
   lista2 <- rep(max(lista)+1,length(nl2))
   sht2 <- rep(NA,length(nl2))
   names(sht2)<-nl2
   names(lista2)<-nl2
   lista3<-c(lista,lista2)
   lista2<-lista3[order(as.numeric(names(lista3)))]
   sht3<-c(sht,sht2)
   sht2<-sht3[order(as.numeric(names(lista3)))]
   names(sht2)<-years
   sht2<-cbind(sht2,lista2)
   bin.size <- size
   uppdelning <- seq(length(years),1,by=-bin.size)
   tabell <- NULL
   for(i in (length(uppdelning)-1):1){ 
      tabell <- rbind(tabell,tabulate(lista2[(uppdelning[i]):(uppdelning[i+1]+1)],nbins=max(lista2)))
   }
   artal <- rev(c(years)[uppdelning])
   tabell <- cbind(c(artal-(bin.size/2))[-1],tabell)
   k <- 0
   sum.prb <- fargning <- tabell
   for(i in (length(uppdelning)-1):1){
      temp <- sht2[(uppdelning[i]):(uppdelning[i+1]+1),]
      k <- k+1
      for(j in 2:(dimension[2]+1)) {
         sum.prb[k,j] <- sum(temp[temp[,2]==(j-1),1],na.rm=T)
      }
   }
   for(j in 2:(dimension[2]+1)) sum.prb[,j] <- sum.prb[,j]/tabell[,j]
   sum.prb <- ifelse(is.nan(sum.prb),0,sum.prb)
   # farger<-c('darkgreen','lightgreen','orange','indianred','brown','darkred')
   prb.mean <- sum.prb
   for(j in 2:(dimension[2]+1)) {
      fargning[,j]<-as.integer(round(10*sum.prb[,j]))
   }
   fargning<-ifelse(fargning==0,1,fargning)
   #minv<-min(fargning[,2:(dimension[2]+1)],na.rm=T)
   #for(j in 2:(dimension[2]+1)) {
   #   fargning[,j]<-fargning[,j]-minv+1
   #}
   list(tabell=tabell,mean.sannolikheter=sum.prb,sannolikheter=sht2,fargning=fargning)
} ## end nysummering()


plotClusterMean <- function(model,Plot=F, lty=2,lwd=3){
## plot the cluster mean in one plot. If Plot=FALSE then the non-plotted values are saved together with the overall mean  
   parameters <- model$pars
   FullS <- model$design$FullS
   K <- model$initials$K
   grid <- model$data$grid
   cluster.mean <- matrix(0,K,dim(FullS)[1])
   for(k in 1:K)
      cluster.mean[k,] <- FullS %*% (parameters$lambda.zero + parameters$Lambda %*% parameters$alpha[k,  ])
   if (Plot){
       print(dim(cluster.mean))
       plot(grid,grid,ylim=range(cluster.mean),type='n',ylab="Cluster Mean")
       grid()
       lines(grid,FullS%*%parameters$lambda.zero,lty=lty,lwd=lwd) 
       for(k in 1:K)
          lines(grid,cluster.mean[k,], col = k+1, lwd = lwd)
   } else
       cbind(grid=grid, overall.mean=FullS%*%parameters$lambda.zero,t(cluster.mean))
} ## end plotClusterMean()


#######################################
## main plot function.....
#######################################

plot.mocca <- function(x,type=1, select =NULL,transform=FALSE,covariance=TRUE, covariates =FALSE, lwd=2,ylab="",xlab="",main="", ylim=NULL,ncolors=NULL,
  probs=FALSE,pts=FALSE,size=50,years=NULL, years.names=NULL, ...)  ## years=c(-4486:1901)           
## creates plots for the mocca model ...
## x' is a mocca object
## type' - determines which type of plots to print. For 'type=1' (default) cluster mean curves are shown in one plot on one page together with the overall mean curve; 'type=2' cluster means are shown on separate plots;  'type=3' illustrates the covariance (or correlation) structure within each cluster.
## select' - the order of the cluster means to be printed  with 'type=1'. If NULL (default), the cluster mean curves are in {1,2,...,K} order, where K is the number of clusters.
## covariance' - whether covariance or correlation matrices should be plotted when printing covariance (type=3),
## covariates' - whether covariates should be included when printing covariance (type=3),
## transform' - whether svd back-transformation of the spline model matrix should be applied when printing covariance (type=3)
## if "type=2":
## the function produces the trend of the frequencies of the different clusters, together with mean probabilites (if probs=T)
## It also gives the mean value of the included covariates within each cluster (not the model estimated covariate values).
## If pts=T then it plots points of the frequency trend.
## size' is the bin size (here 50 years) of how many of those 50 years belong to a specific cluster.
## select' can be given if you want to plot the clusters in a different select or just a few of them.
## j' can be used to tell which cluster it starts with when you plot a few of the clusters,
## ex., if you have 6 clusters and you want to plot cluster 2:4 then you set select=2:4 and j=1.
## stand.cov' tells whether standardized covariates were used in the modelling.
## years' gives what years you want to calculate frequencies in. 
{ 
  oldpar <- par(no.readonly = TRUE)    
  on.exit(par(oldpar))  
  if (type==1) { ## plotting all cluster means in one plot...
      # old.par<-par(mfrow=c(1,1))
      par(mfrow=c(1,1))
      illustrations.plot(x=x,select=select,lwd=lwd,ylab=ylab,ylim=ylim, ncolors=ncolors)
   } else if (type==2) { ## plotting cluster means on separate plots...
         if (is.null(years)) stop("'years' should be supplied. type=2 plot works with yearly data.")
         else{
           if (is.null(select)){
              r <- ceiling(sqrt(x$initials$K))
              if (x$initials$K==2)
                   par(mfrow=c(1,2))  # old.par<- par(mfrow=c(1,2)) 
               else if(x$initials$K==5 | x$initials$K==6) 
                   par(mfrow=c(2,3)) # old.par<- par(mfrow=c(2,3))
                 else par(mfrow=c(r,r)) # old.par<- par(mfrow=c(r,r))
           } else{
                  r <- ceiling(sqrt(length(select)))
                  if (length(select)==2)
                      par(mfrow=c(2,1)) # old.par<- par(mfrow=c(2,1))
                  else par(mfrow=c(r,r)) # old.par<- par(mfrow=c(r,r))
             } 
           trendplot(x=x,probs=probs,pts=pts,size=size,select=select,years=years,ylim=ylim, years.names=years.names)
          }
      } else if (type==3) { ## plotting covariance structures...
           r <- ceiling(sqrt(x$initials$K))
           if (x$initials$K==2)
                  par(mfrow=c(1,2))  # old.par<- par(mfrow=c(1,2)) 
             else if(x$initials$K==5 | x$initials$K==6) 
                   par(mfrow=c(2,3)) # old.par<- par(mfrow=c(2,3))
             else par(mfrow=c(r,r)) # old.par<- par(mfrow=c(r,r))
           covariance.illustration(x=x,transform=transform,covariance=covariance,covariates=covariates)  
        } else stop("not sure what to plot. 'type' must be either 1, 2, or 3.")

     # par(old.par)
} ## end plot.mocca




