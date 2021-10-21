\name{estimate.mocca}
\alias{estimate.mocca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Model parameter estimation}
\description{
Function to estimate model parameters by maximizing the observed log likelihood via an EM algorithm. 
The estimation procedure is based on an algorithm proposed by James and Sugar (2003).

The function is not normally called directly, but rather service routines for \code{\link{mocca}}. See the description of the \code{\link{mocca}} function for more detailed information of arguments.
}
\usage{
estimate.mocca(data,K=5,q=6,h=2,random=TRUE,B=NULL,svd=TRUE,
       use.covariates=FALSE,stand.cov=TRUE,index.cov=NULL,
       lambda=1.4e-4,EM.maxit=50, EMstep.tol=1e-8,Mstep.maxit=10,
       Mstep.tol=1e-4, EMplot=TRUE,trace=TRUE,n.cores=NULL)

%- maybe also 'usage' for other objects documented here.
}
\arguments{
  \item{data}{a list containing at least five objects (vectors) named as \code{x}, \code{time}, \code{timeindex}, \code{curve}, \code{grid}, \code{covariates} (optional). See \code{\link{mocca}} for the detailed explanation of each object.}

  \item{K}{number of clusters (default: \code{K=3}).}

   \item{q}{number of B-splines used to describe the individual curves. Evenly spaced knots are used (default: \code{q=6}). (currently only B-splines are implemented, however, it is possible  to use other basis functions such as, e.g. Fourier basis functions)} 

  \item{h}{a positive integer, parameter vector dimension in low-dimensionality representation of the curves (spline coefficients).
          \eqn{h} should be less or equal to the number of clusters \eqn{K} (default: \code{h=2}).}

  \item{random}{\code{TRUE/FALSE},  if \code{TRUE} each subject is randomly assigned to one of the \code{K} clusters initially, otherwise \code{k-}means is used to initialize cluster belongings (default: TRUE).}

   \item{B}{ an \eqn{N x q} matrix of spline coefficients,  the spline approximation of the yearly curves based on \eqn{p} number of splines. If \code{B=NULL} (default), the coefficients are estimated using \code{fda:: create.bspline.basis}.} 

  \item{svd}{\code{TRUE/FALSE}, whether SVD decomposition should be used for the matrix of spline coefficients (default: TRUE).}

   \item{use.covariates}{\code{TRUE/FALSE}, whether covariates should be included when modelling (default: FALSE).} 

   \item{stand.cov}{\code{TRUE/FALSE}, whether covariates should be standardized when modelling (default: \code{TRUE}).}

  \item{index.cov}{a vector of indices indicating which covariates should be used when modelling. If \code{NULL} (default) all present covariates are included.}

  \item{lambda}{a positive real number, smoothing parameter value to be used when estimating B-spline coefficients.}

   \item{EM.maxit}{a positive integer which gives the maximum number of iterations for a EM algorithm (default: EM.maxit=50).} 

  \item{EMstep.tol}{the tolerance to use within iterative procedure of the EM algorithm (default: EMstep.tol=1e-8).}

   \item{Mstep.maxit}{a positive scalar which gives the maximum number of iterations for an inner loop of the parameter estimation in M step (default: Mstep.maxit=20).} 

  \item{Mstep.tol}{the tolerance to use within iterative procedure to estimate model parameters (default: Mstep.tol=1e-4).}

   \item{EMplot}{\code{TRUE/FALSE}, whether plots of cluster means with some summary information should be produced at each iteration of the EM algorithm (default: FALSE).} 

  \item{trace}{\code{TRUE/FALSE}, whether to print the current values of \eqn{\sigma^2} and \eqn{\sigma^2_x} of the covariates at each iteration of \code{M step}  (default: TRUE).}

  \item{n.cores}{number of cores to be used with parallel computing.}
 
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{A list is returned with the following items:

  \item{loglik}{the maximized log likelihood value.} 
 
 \item{sig2}{estimated residual variance for the spline coefficients (for the model without covariates), or a vector of the estimated residual variances for the spline coefficients and for the covariates (for the model with covariates).}

  \item{conv}{indicates why the EM algorithm terminated:
            
             0: indicates successful completion.

             1: indicates that the iteration limit \code{EM.maxit} has been reached.
         } 
 
  \item{iter}{number of iterations of the EM algorithm taken to get convergence.}

  \item{score.hist}{a matrix of the succesive values of the scores: residual variances and log likelihood, up until convergence.}

 \item{parameters}{a list containing all the estimated parameters: \eqn{\bm\lambda_0}, \eqn{\bm\Lambda}, \eqn{\bm\alpha_k}, \eqn{\bm\Gamma_k} (or \eqn{\bm\Delta_k} in presence of covariates), \eqn{\pi_k} (probabilities of cluster belongnings), \eqn{\sigma^2_x} (residual variance for the covariates if present), \eqn{\mathbf{v}_k} (mean values of the covariates for each cluster, in presence of covariates), \eqn{k=1,..., K}, where \eqn{K} is the number of clusters.} 

  \item{vars}{a list containing results from the E step of the algorithm: the posterior probabilities for each subject \eqn{\pi_{k|i}}'s, the expected values of the \eqn{\bm\gamma_i}'s, \eqn{\bm\gamma_i\bm\gamma_i^T}, and the covariance matrix of \eqn{\bm\gamma_i} given cluster membership and the observed values of the curve. See Arnqvist and \enc{Sjöstedt}{Sjostedt} de Luna (2019) that explains these values.}

   \item{data}{a list containing all the original data plus re-arranged functional data and covariates (if supplied) needed for EM-steps.} 

  \item{design}{a list of spline basis matrices with and without covariates: 
   \code{FullS.bmat} is the spline basis matrix \eqn{\mathbf{S}} computed on the grid of uniquily specified time points; 
 \code{FullS} is the spline basis matrix \code{FullS.bmat} or \eqn{\mathbf U} matrix from the svd of \code{FullS} (if applied);
\eqn{\mathbf{S}} is the spline basis matrix computed on \code{timeindex}, a vector of time indices from \eqn{T} possible from \code{grid};
the inverse \eqn{(\mathbf{S}^T\mathbf{S})^{-1}};
\code{tag.S} is the matrix \eqn{\mathbf{S}} with covariates; 
\code{tag.FullS} is the matrix \code{FullS} with covariates. 
   }

  \item{initials}{a list of initial settings: \eqn{q} is the spline basis dimension, \eqn{N} is the number of objects/curves, \eqn{Q} is the number of basis dimension plus the number of covariates (if present), \eqn{random} is whether k-means was used to initialize cluster belonings, \eqn{h} is the vector dimension in low-dimensionality representation of the curves, \eqn{K} is the number of clusters,  \eqn{r} is the number of scalar covariates.}
  
}
\references{
James, G.M., Sugar, C.A. (2003). Clustering for sparsely sampled functional data.  \emph{Journal of the American Statistical Association}, 98.462, 397--408.

Arnqvist, P., and \enc{Sjöstedt}{Sjostedt} de Luna, S. (2019). Model based functional clustering of varved lake sediments. \emph{arXiv preprint arXiv:1904.10265}.
}
\author{
Per Arnqvist, Natalya Pya Arnqvist, Sara \enc{Sjöstedt}{Sjostedt} de Luna
}
\seealso{\code{\link{fdaMocca-package}}, \code{\link{mocca}}
}
\keyword{fda} 
\keyword{covariates} 
\keyword{EM algorithm}



