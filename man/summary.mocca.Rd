\name{summary.mocca}
%\Rdversion{0.1-0}
\alias{summary.mocca}
\alias{print.summary.mocca}

\title{Summary for a mocca fit
}
\description{Takes a \code{mocca} object produced by \code{mocca()} and produces various useful
summaries from it. 
}
\usage{
\method{summary}{mocca}(object,...)

\method{print}{summary.mocca}(x,digits = max(3, getOption("digits") - 3),...)
}

\arguments{
\item{object}{a fitted \code{mocca} object as produced by \code{mocca()}.}

\item{x}{a \code{summary.mocca} object produced by \code{summary.mocca()}.} 

\item{digits}{controls the number of digits printed in the output.}

\item{...}{ other arguments.}

}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }Filter


\value{
\code{summary.mocca} produces the following list of summary information for a \code{mocca} object. 


\item{N}{number of observations}
\item{K}{number of clusters}
\item{r}{number of scalar covariates if model with covariates}
\item{sig2}{residual variance estimate for the functional data and for the scalar covariates (if the model is with covariates)}
\item{p}{total number of the estimated parameters in the model}
\item{tab_numOfCurves_cluster}{number of objects/curves in each cluster as a table. Here 'hard' clustering is applied, where each object/curve belongs to a cluster with the highest posterior probability.}
\item{covariates_est}{mean value estimates for scalar covariates given cluster belongings (if the model is with covariates)}
\item{t.probs}{estimated probabilities of belonging to each cluster}
\item{crita}{a table with the maximized log likelihood, AIC, BIC and Shannon entropy values of the fitted model}
}


\author{
Per Arnqvist, Sara \enc{Sjöstedt}{Sjostedt} de Luna,  Natalya Pya Arnqvist
}
%\references{

%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{mocca}} 
}

\examples{
## see ?mocca help files
}

\keyword{fda} \keyword{clustering} 


