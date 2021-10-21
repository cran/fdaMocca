\name{criteria.mocca}
\alias{criteria.mocca}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{ AIC, BIC, entropy for a functional clustering model
}
\description{Function to extract the information criteria AIC and BIC, as well as the average Shannon entropy over all functional objects, for a fitted functional clustering \code{mocca}. The Shannon entropy is computed over the posterior probability distribution of belonging to a specific cluster given the functional object (see Arnqvist and Sjöstedt de Luna, 2019, for further details).  

}
\usage{
  criteria.mocca(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x}{fitted model objects of class \code{mocca} as produced by \code{mocca()}.}

}

%\details{ 
%}

\value{A table with the AIC, BIC and Shannon entropy values of the fitted model. 
}
\references{
Arnqvist, P., and \enc{Sjöstedt}{Sjostedt} de Luna, S. (2019). Model based functional clustering of varved lake sediments. \emph{arXiv preprint arXiv:1904.10265}.

}

\author{
Per Arnqvist
}

\seealso{  \code{\link{logLik.mocca}}, \code{\link{mocca}} 
}
\examples{
## see examples in mocca() 
}


\keyword{clustering} \keyword{fda} \keyword{covariates}%-- one or more ...






