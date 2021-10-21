\name{mocca}
% \Rdversion{0.1-0}
\alias{mocca}
\title{Model-based clustering for functional data with covariates}

\description{This function fits a functional clustering model to observed independent functional subjects, where a functional subject consists of a function and possibly a set of covariates. Here, each curve is projected onto a finite dimensional basis and clustering is done on the resulting basis coefficients. However, rather than treating basis coefficients as parameters, mixed effect modelling is used for the coefficients. In the model-based functional clustering approach the functional subjects (i.e. the spline/basis coefficients and the potential covariates) are assumed to follow a multivariate Gaussian mixture model, where the number of distributions in the mixture model corresponds to the number of (predefined) clusters, \eqn{K}. Given that a functional subject belongs to a cluster \eqn{k}, the basis coefficients and covariate values are normally distributed with a cluster-specific mean and covariance structure. 

An EM-style algorithm based on James and Sugar (2003) is implemented to fit the Gaussian mixture model for a prespecified number of clusters \code{K}. The model allows for different cluster covariance structure for the spline coefficients and model coefficients for the covariates. See Arnqvist and \enc{Sjöstedt}{Sjostedt} de Luna (2019) for details about differences to the clustering model and its implementation. 

The routine calls \code{estimate.mocca} for the model fitting. 


}
\usage{
mocca(data=stop("No data supplied"), K = 5, q = 6, h = 2,
     use.covariates=FALSE,stand.cov=TRUE,index.cov=NULL,
     random=TRUE, B=NULL,svd=TRUE,lambda=1.4e-4, EM.maxit=50, 
     EMstep.tol=1e-6,Mstep.maxit=20,Mstep.tol=1e-4,EMplot=TRUE,
     trace=FALSE,n.cores=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a list containing at least three objects (vectors) named as \code{x}, \code{time}, \code{curve}, and optional \code{timeindex}, \code{grid} and \code{covariates}:

       i) suppose we observe \eqn{N} independent subjects, each consisting of a curve and potentially a set of scalar covariates, where the \eqn{i^{th}} curve has been observed at \eqn{n_i} different time points, \eqn{i=1,...,N}. \code{x} is a vector of length \eqn{\sum_{i=1}^N n_i} with the first \eqn{n_1} elements representing the observations of the first curve, followed by \eqn{n_2} observations of the second curve, etc;

      ii)  \code{time} is a \eqn{\sum_i n_i} vector of the concatenated time points for each curve \eqn{(t_{ij}, j=1,...,n_i, i=1,...,N)}, with the first \eqn{n_1} elements being the time points at which the first curve is observed, etc. Often, the time points within each curve are scaled to \eqn{[0,1]}. 
%- of time points (usually time points are scaled to \eqn{[0,1]} for each subject);

      iii) \code{timeindex} is a \eqn{\sum_i n_i} vector of time indices from \eqn{T} possible from \code{grid}. So each observation has a corresponding location (time index) within \eqn{[0,1]} uniquely specified time points. If not supplied, obtained from \code{time} and \code{grid};

      iv)  \code{curve} is a \eqn{\sum_i n_i} vector of integers from \eqn{1,..., N}, specifying the subject number for each observation in \code{x};

       v)  \code{grid} is a \eqn{T} vector of all unique time points (values within \eqn{[0,1]} interval) for all \eqn{N} subjects, needed for estimation of the B-spline coefficients in \code{fda::eval.basis()}. \code{timeindex} and \code{grid} together give the timepoint for each subject (curve). If not supplied, obtained from \code{time}.
       
       vi) if supplied, \code{covariates} is an \eqn{N \times r} matrix (or data frame) of scalar covariates (finite-dimensional covariates).

   } 

  \item{K}{number of clusters (default: \code{K=3}).}

   \item{q}{number of B-splines for the individual curves. Evenly spaced knots are used (default: \code{q=6}).} 

  \item{h}{a positive integer, parameter vector dimension in the low-dimensionality representation of the curves (spline coefficients).
          \eqn{h} should be smaller than the number of clusters \eqn{K} (default: \code{h=2}).}

  \item{use.covariates}{\code{TRUE/FALSE}, whether covariates should be used when modelling (default: \code{FALSE}).} 

  \item{stand.cov}{\code{TRUE/FALSE}, whether covariates should be standardized when modelling (default: \code{TRUE}).}

  \item{index.cov}{a vector of indices indicating which covariates should be used when modelling. If \code{NULL} (default) all present covariates are included.}

  \item{random}{\code{TRUE/FALSE},  if \code{TRUE} the initial cluster belongings is given by uniform distribution, otherwise \code{k-}means is used to initialize cluster belongings (default: \code{TRUE}).}

   \item{B}{ an \eqn{N \times q} matrix of spline coefficients,  the spline approximation of the yearly curves based on \eqn{p} number of splines. If \code{B=NULL} (default), the coefficients are estimated using \code{fda:: create.bspline.basis}.} 

  \item{svd}{\code{TRUE/FALSE}, whether SVD decomposition should be used for the matrix of spline coefficients (default: \code{TRUE}).}

  \item{lambda}{a positive real number, smoothing parameter value to be used when estimating B-spline coefficients.}

   \item{EM.maxit}{a positive integer which gives the maximum number of iterations for a EM algorithm (default:
 \code{EM.maxit=50}).} 

  \item{EMstep.tol}{the tolerance to use within iterative procedure of the EM algorithm (default: \code{EMstep.tol=1e-8}).}

   \item{Mstep.maxit}{a positive scalar which gives the maximum number of iterations for an inner loop of the parameter estimation in M step (default: \code{Mstep.maxit=20}).} 

  \item{Mstep.tol}{the tolerance to use within iterative procedure to estimate model parameters (default: Mstep.tol=1e-4).}

   \item{EMplot}{\code{TRUE/FALSE}, whether plots of cluster means with some summary information should be produced at each iteration of the EM algorithm (default: \code{FALSE}).} 

  \item{trace}{\code{TRUE/FALSE}, whether to print the current values of \eqn{\sigma^2} and \eqn{\sigma^2_x} for the covariates at each iteration of \code{M step}  (default: \code{FALSE}).}

  \item{n.cores}{number of cores to be used with parallel computing. If \code{NULL} (default) \code{n.cores} is set to the  numbers of available cores - 1 (\code{n.cores= detectCores()-1}). }

}

\details{
A model-based clustering with covariates (mocca) for the functional subjects (curves and potentially covariates) is a gaussian mixture model with \eqn{K} components. 
Let \eqn{g_i(t)} be the true function (curve) of the \eqn{i^{th}} subject, for a set of \eqn{N} independent subjects. Assume that for each subject we have a vector of observed values of the function \eqn{g_i(t)} at times \eqn{t_{i1},...,t_{in_i},} obtained with some measurement errors. We are interested in clustering the subjects into \eqn{K} (homogenous) groups. Let \eqn{y_{ij}} be the observed value of the \eqn{i}th curve at time point \eqn{t_{ij}}. Then  
\deqn{
y_{ij} = g_i(t_{ij})+ \epsilon_{ij},  i=1,...,N, j=1,...,n_i,}
where \eqn{\epsilon_{ij}} are assumed to be independent and normally distributed measurement errors
with mean \eqn{0} and variance \eqn{\sigma^2}. Let \eqn{\mathbf{y}_i}, \eqn{\mathbf{g}_i,} and \eqn{\boldsymbol{\epsilon}_i} be the \eqn{n_i}-dimensional vectors for subject \eqn{i}, corresponding to the observed values, true values and measurement errors, respectively. Then, in matrix notation, the above could be written as 
\deqn{
 \mathbf{y}_i=\mathbf{g}_i+\boldsymbol{\epsilon}_i, ~~~~i=1,\ldots, N,
}
where \eqn{\boldsymbol{\epsilon}_i ~\sim ~ N_{n_i}(\mathbf{0},\sigma^2 \mathbf{I}_{n_i}).} 
We further assume that the smooth function \eqn{g_i(t)} can be expressed as
\deqn{ 
g_i(t) = \boldsymbol{\phi}^T(t) \boldsymbol{\eta}_i,
}
where \eqn{\boldsymbol{\phi}(t)=\left(\phi_{1}(t),\ldots,\phi_{p}(t)\right)^T} 
is a \eqn{p}-dimensional vector of known basis functions evaluated at time t, e.g. B-splines, and \eqn{\boldsymbol{\eta}_i} a \eqn{p}-dimensional vector of unknown (random) coefficients. The \eqn{\boldsymbol{\eta}_i}'s are modelled as
\deqn{
\boldsymbol{\eta}_i = \boldsymbol{\mu}_{z_i} + \boldsymbol{\gamma}_i, ~~~ \boldsymbol{\eta}_i ~ \sim ~ N_p(\boldsymbol{\mu}_{z_i},\bm{\Gamma}_{z_i}),
}
where \eqn{\boldsymbol{\mu}_{z_i}} is a vector of expected spline coefficients for a cluster \eqn{k} and \eqn{z_i} denotes the unknown cluster membership, with \eqn{P(z_i=k)=\pi_k},  \eqn{k=1,\ldots,K}. The random vector \eqn{\boldsymbol{\gamma}_i} corresponds to subject-specific within-cluster variability.
 Note that this variability is allowed to be different in different clusters, due to \eqn{\bm\Gamma_{z_i}}. If desirable, given that subject \eqn{i} belongs to cluster \eqn{z_i=k}, a further parametrization of \eqn{\boldsymbol{\mu}_{k},~~ k=1,\ldots,K,} may prove useful, for producing low-dimensional representations of the curves as suggested by James and Sugar (2003):
\deqn{
\bm\mu_k = \bm\lambda_0+ \bm\Lambda \bm\alpha_k,
}
where \eqn{\bm\lambda_0} and \eqn{\bm\alpha_k} are \eqn{p}- and \eqn{h}-dimensional vectors respectively, and \eqn{\bm\Lambda} is a \eqn{p \times h} matrix with \eqn{h \leq K-1}. Choosing \eqn{h<K-1} may be valuable, especially for sparse data. In order to ensure identifiability, some restrictions need to be put on the parameters. Imposing the restriction that \eqn{\sum_{k=1}^K \bm\alpha_k=\mathbf{0}} implies that \eqn{\bm\phi^T(t)\bm\lambda_0} can be viewed as the overall mean curve. Depending on the choice of \eqn{h,p} and \eqn{K}, further restrictions may need to be imposed in order to have identifiability of the parameters (\eqn{\bm\lambda_0, \bm\Gamma} and \eqn{\bm\alpha_k} are confounded if no restrictions are imposed).
In vector-notation we thus have
\deqn{
\mathbf{y}_i = \mathbf{B}_i(\bm\lambda_0 +  \bm\Lambda\bm\alpha_{z_i}+\bm\gamma_i)+\bm\epsilon_i,~~ i=1,...,N,
}
where \eqn{\mathbf{B}_i} is an \eqn{n_i \times p} matrix with \eqn{\bm\phi^T(t_{ij})} on the \eqn{j^\textrm{th}} row, \eqn{j=1,\ldots,n_i.}
We will also assume that the \eqn{\bm\gamma_i}'s, \eqn{\bm\epsilon_i}'s and the \eqn{z_i}'s are independent. Hence, given that subject \eqn{i} belongs to cluster \eqn{z_i=k} we have
\deqn{
\mathbf{y}_i | z_i=k ~~\sim ~~ N_{n_i}\left(\mathbf{B}_i(\bm\lambda_0 +  \bm\Lambda \bm\alpha_k), ~~\mathbf{B}_i \bm\Gamma_k \mathbf{B}_i^T+ \sigma^2\mathbf{I}_{n_i}\right).
}
Based on the observed data \eqn{\mathbf{y}_1,\ldots,\mathbf{y}_N}, the parameters \eqn{\bm\theta}
of the model can be estimated by maximizing the observed likelihood
\deqn{
L_o(\bm\theta|\mathbf{y}_1,\ldots,\mathbf{y}_N)=\prod_{i=1}^N \sum_{k=1}^G \pi_k f_k(\mathbf{y}_i,\bm\theta),
}
where
%\eqn{\bm\theta = \left\{\bm\lambda_0,\bm\Lambda,\bm\alpha_1,\ldots,\bm\alpha_K,\pi_1,\ldots,\pi_K,\sigma^2_c,\rho_1,\ldots,\rho_K,\sigma^2,\bm\Gamma_1,\ldots,\bm\Gamma_K\right\},}
\eqn{\bm\theta = \left\{\bm\lambda_0,\bm\Lambda,\bm\alpha_1,\ldots,\bm\alpha_K,\pi_1,\ldots,\pi_K,\sigma^2,\bm\Gamma_1,\ldots,\bm\Gamma_K\right\},}
and \eqn{f_k(\mathbf{y}_i,\bm\theta)} is the normal density given above. Note that here \eqn{\bm\theta} will denote all scalar, vectors and matrices of parameters to be estimated. An EM-type algorithm is used to maximize the likelihood above.


If additional covariates  have been observed for each subject besides the curves, they can also be included in the model when clustering the subjects. Given that the subject \eqn{i} belongs to cluster \eqn{k, (z_{i}=k)} the \eqn{r} covariates \eqn{\boldsymbol{x}_i \in \mathbf{R}^r} are assumed to have mean value \eqn{\boldsymbol{\upsilon}_k} and moreover \eqn{\boldsymbol{x}_{i} = \boldsymbol{\upsilon}_{k} + \boldsymbol{\delta}_{i} + \boldsymbol{e}_i,}
where we assume that \eqn{\boldsymbol{\delta}_{i}|z_{i}=k \sim N_r(\boldsymbol{0}, \mathbf{D}_k)} is the random deviation within cluster and \eqn{\boldsymbol{e}_i \sim N_r(\boldsymbol{0},\sigma_x^2 \mathbf{I}_r)} independent remaining unexplained variability. 
Note that this model also incorporates the dependence between covariates and the random curves via the random basis coefficients. See Arnqvist and \enc{Sjöstedt}{Sjostedt} de Luna (2019) for further details.
EM-algorithm is implemented to maximize the mixture likelihood.


 The method is applied to annually varved lake sediment data from the lake \enc{Kassjön}{Kassjon} in Northern Sweden. See an example and also \code{\link{varve}} for the data description.
}

\value{The function returns an object of class \code{"mocca"} with the following elements:

 \item{loglik}{the maximized log likelihood value.} 
 
 \item{sig2}{estimated residual variance for the functional data (for the model without covariates), or a vector of the estimated residual variances for the functional data and for the covariates (for the model with covariates).}

  \item{conv}{indicates why the EM algorithm terminated:
            
             0: indicates successful completion.

             1: indicates that the iteration limit \code{EM.maxit} has been reached.
         } 
 
  \item{iter}{number of iterations of the EM algorithm taken to get convergence.}

  \item{nobs}{number of subjects/curves.}

  \item{score.hist}{a matrix of the succesive values of the scores, residual variances and log likelihood, up until convergence.}

 \item{pars}{a list containing all the estimated parameters: \eqn{\bm\lambda_0}, \eqn{\bm\Lambda}, \eqn{\bm\alpha_k}, \eqn{\bm\Gamma_k} (or \eqn{\bm\Delta_k} in presence of the covariates), \eqn{\pi_k} (probabilities of cluster belongnings), \eqn{\sigma^2}, \eqn{\sigma^2_x} (residual variance for the covariates if present), \eqn{\mathbf{v}_k} (mean values of the covariates for each cluster).} 

  \item{vars}{a list containing results from the E step of the algorithm: the posterior probabilities for each subject \eqn{\pi_{k|i}}'s, the expected values of the \eqn{\bm\gamma_i}'s, \eqn{\bm\gamma_i\bm\gamma_i^T}, and the covariance matrix of \eqn{\bm\gamma_i} given cluster membership and the observed values of the curve. See Arnqvist and \enc{Sjöstedt}{Sjostedt} de Luna (2019) that explains these values.}

   \item{data}{a list containing all the original data plus re-arranged functional data and covariates (if supplied).} 

  \item{design}{a list of spline basis matrices with and without covariates: 
   \code{FullS.bmat} is the spline basis matrix \eqn{\mathbf{S}} computed on the grid of uniquily specified time points; 
 \code{FullS} is the spline basis matrix \code{FullS.bmat} or \eqn{\mathbf U} matrix from the svd of \code{FullS} (if applied);
\eqn{\mathbf{S}} is the spline basis matrix computed on \code{timeindex}, a vector of time indices from \eqn{T} possible from \code{grid};
the inverse \eqn{(\mathbf{S}^T\mathbf{S})^{-1}};
\code{tag.S} is the matrix \eqn{\mathbf{S}} with covariates; 
\code{tag.FullS} is the matrix \code{FullS} with covariates. See Arnqvist and Sjöstedt de Luna (2019) for further details.
   }

  \item{initials}{a list of initial settings: \eqn{q} is the spline basis dimension, \eqn{N} is the number of subjects/curves, \eqn{Q} is the number of basis dimension plus the number of covariates (if present), \eqn{random} is whether k-means was used to initialize cluster belonings, \eqn{h} is the vector dimension in low-dimensionality representation of the curves, \eqn{K} is the number of clusters,  \eqn{r} is the number of scalar covariates, \eqn{moc} TRUE/FALSE signaling if the model includes covariates.}
  
}

\references{
Arnqvist, P., Bigler, C., Renberg, I., \enc{Sjöstedt}{Sjostedt} de Luna, S. (2016). Functional clustering of varved lake sediment to reconstruct past seasonal climate. \emph{Environmental and Ecological Statistics}, \bold{23}(4), 513--529.

Abramowicz, K., Arnqvist, P., Secchi, P., \enc{Sjöstedt}{Sjostedt} de Luna, S., Vantini, S., Vitelli, V. (2017). Clustering misaligned dependent curves applied to varved lake sediment for climate reconstruction. \emph{Stochastic Environmental Research and Risk Assessment}. Volume \bold{31.1}, 71--85.

Arnqvist, P., and \enc{Sjöstedt}{Sjostedt} de Luna, S. (2019). Model based functional clustering of varved lake sediments. \emph{arXiv preprint arXiv:1904.10265}.

James, G.M., Sugar, C.A. (2003). Clustering for sparsely sampled functional data.  \emph{Journal of the American Statistical Association}, 98.462, 397--408.
}
\author{
Per Arnqvist, Natalya Pya Arnqvist, Sara \enc{Sjöstedt}{Sjostedt} de Luna
}
\seealso{\code{\link{fdaMocca-package}}%, \code{\link{plot.mocca}}, \code{\link{summary.mocca}}
}

\examples{
 \donttest{
## example with lake sediment data from lake Kassjön...
library(fdaMocca)
data(varve) ## reduced data set

## run without covariates...
m <- mocca(data=varve,K=3,n.cores=2)
m
## some summary information...
summary(m)
criteria.mocca(m)
AIC(m)
BIC(m)
## various plots...
plot(m)
plot(m,select=2)
plot(m,type=2,years=c(-750:750)) 
plot(m,type=2,probs=TRUE,pts=TRUE,years=c(-750:750)) 
plot(m,type=2,pts=TRUE,select=c(3,1),years=c(-750:750))
plot(m,type=3)
plot(m,type=3,covariance=FALSE)


## model with two covariates...
## note, it takes some time to analyze the data...
m1 <- mocca(data=varve, use.covariates=TRUE,index.cov=c(2,3), K=3,n.cores=2)
m1
## summary information...
summary(m1)
criteria.mocca(m1)
## various plots...
plot(m1)
plot(m1,type=2,pts=TRUE,years=c(-750:750)) 
plot(m1,type=3)
plot(m1,type=3,covariance=FALSE)
plot(m1,type=3,covariates=TRUE)

## simple simulated data...
data(simdata)
set.seed(2)
m2 <- mocca(data=simdata,K=2,q=8,h=1,lambda=1e-10,n.cores=2,EMstep.tol=1e-3)
summary(m2)
criteria.mocca(m2)
plot(m2)
plot(m2,select=2)
}

## even simpler simulated data
##(reduced from 'simdata', EMstep.tol set high, q lower to allow automatic testing)...
library(fdaMocca)
data(simdata0)
set.seed(2)
m3 <- mocca(data=simdata0,K=2,q=5,h=1,lambda=1e-10,n.cores=2,EMstep.tol=.5,
      EMplot=FALSE,B=simdata0$B)
summary(m3)
#plot(m3)
#plot(m3,select=2))
 
}

\keyword{models} \keyword{functiona data analysis} \keyword{cluster analysis} \keyword{random effect}%-- one or more ..
\concept{EM algorithm}
\concept{B-spline}


