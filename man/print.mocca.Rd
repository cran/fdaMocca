\name{print.mocca}
\alias{print.mocca}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Print a mocca object}
\description{ The default print method for a \code{mocca} object.
}
\usage{
\method{print}{mocca}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x, ...}{ fitted model objects of class \code{mocca} as produced by \code{mocca()}.}
}
\details{ Prints out whether the model is fitted with or without covariates, the number of clusters, the estimated residual variance for the functional data and for the scalar covariates (if present), the number of covariates (if present), the maximized value of the log likelihood, and the number of subjects/curves.

}
\value{ No return value, the function prints out several fitted results.

}

\author{
Per Arnqvist, Natalya Pya Arnqvist, Sara \enc{Sjöstedt}{Sjostedt} de Luna
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%    
%}

%\examples{

%}











