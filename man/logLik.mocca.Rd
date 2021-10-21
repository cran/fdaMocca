\name{logLik.mocca}
\alias{logLik.mocca}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Log-likelihood for a functional clustering model
}
\description{Function to extract the log-likelihood for a fitted functional clustering \code{mocca}
model (fitted by mixture likelihood maximization). 

Note: \code{estimate.mocca} uses \code{loglik.EMmocca()} for calculating the log likelihood at each iterative step.

}
\usage{
  \method{logLik}{mocca}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{fitted model objects of class \code{mocca} as produced by \code{mocca()}.}
\item{...}{unused in this case} 

}

%\details{ 
%}

\value{The log-likehood value as \code{logLik} object. 
}
\references{
Arnqvist, P., and \enc{Sjöstedt}{Sjostedt} de Luna, S. (2019). Model based functional clustering of varved lake sediments. \emph{arXiv preprint arXiv:1904.10265}.

}

\author{
Per Arnqvist
}

\seealso{  \code{\link{estimate.mocca}}, \code{\link{mocca}} 
}

\keyword{clustering} \keyword{fda} \keyword{covariates}%-- one or more ...






